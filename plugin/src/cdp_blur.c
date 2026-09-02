/* CDP Blur — spectral time-smear as a CLAP effect.
 *
 * Port of CDP "blur blur" mode (dev/blur/ap_blur.c): each new analysis frame
 * becomes the target of a linear interpolation ramp over 'span' frames —
 * amp prev->next, phase prev->next (unwrapped), matching CDP's do_the_bltr.
 * CDP runs on .pvx analysis files; here analysis/resynthesis is done in-plugin
 * via radix-2 FFT with Hann window, HOP = N/4, exact COLA normalization.
 */

#include <clap/clap.h>
#include <clap/ext/thread-pool.h>

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define FFT_N 1024
#define HOP   (FFT_N / 4)
#define BINS  (FFT_N / 2 + 1)
#define MAX_CH 2

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

enum { P_SPAN = 0, P_MODE = 1, P_AVGW = 2, P_SUPN = 3, P_NOISE = 4,
       P_CHOR_AMPR = 5, P_CHOR_FRQR = 6, P_CHOR_DIR = 7,
       P_SCAT_KEEP = 8, P_SCAT_BLOK = 9, NUM_PARAMS };

enum { MODE_BLUR = 0, MODE_AVRG = 1, MODE_SUPPRESS = 2, MODE_NOISE = 3,
       MODE_CHORUS = 4, MODE_SCATTER = 5 };

/* ---------------- FFT (in-place radix-2, forward when dir<0) ---------------- */

static float s_tw_re[FFT_N / 2], s_tw_im[FFT_N / 2];
static int   s_rev[FFT_N];
static int   s_fft_ready = 0;

static void fft_init(void) {
    if (s_fft_ready) return;
    for (int i = 0; i < FFT_N / 2; ++i) {
        s_tw_re[i] = cosf(2.0f * (float)M_PI * i / FFT_N);
        s_tw_im[i] = -sinf(2.0f * (float)M_PI * i / FFT_N); /* forward */
    }
    for (int i = 0; i < FFT_N; ++i) {
        int r = i, b = 0;
        for (int k = 0; k < 10; ++k) { b = (b << 1) | (r & 1); r >>= 1; }
        s_rev[i] = b;
    }
    s_fft_ready = 1;
}

static void fft(float *re, float *im, int inverse) {
    for (int i = 0; i < FFT_N; ++i)
        if (s_rev[i] > i) {
            float tr = re[i]; re[i] = re[s_rev[i]]; re[s_rev[i]] = tr;
            float ti = im[i]; im[i] = im[s_rev[i]]; im[s_rev[i]] = ti;
        }
    for (int len = 2; len <= FFT_N; len <<= 1) {
        int half = len >> 1;
        int step = FFT_N / len;
        for (int i = 0; i < FFT_N; i += len)
            for (int j = 0; j < half; ++j) {
                float wr = s_tw_re[j * step];
                float wi = inverse ? -s_tw_im[j * step] : s_tw_im[j * step];
                float xr = re[i + j + half] * wr - im[i + j + half] * wi;
                float xi = re[i + j + half] * wi + im[i + j + half] * wr;
                re[i + j + half] = re[i + j] - xr;
                im[i + j + half] = im[i + j] - xi;
                re[i + j] += xr;
                im[i + j] += xi;
            }
    }
    if (inverse) {
        float s = 1.0f / FFT_N;
        for (int i = 0; i < FFT_N; ++i) { re[i] *= s; im[i] *= s; }
    }
}

/* ---------------- per-channel STFT state ---------------- */

/* per-channel scratch + state (so channels can be processed in parallel) */
typedef struct {
    float  win[FFT_N];            /* time-domain window accumulator */
    float  hop[HOP];              /* sample accumulator until hop full */
    float  ola[FFT_N];            /* overlap-add accumulator */
    float  re[FFT_N], im[FFT_N];  /* FFT work buffers */
    float  mag_prev[BINS], ph_prev[BINS];
    float  mag_tb[BINS], ph_tb[BINS]; /* per-frame polar scratch */
    int    bloks[BINS];           /* scatter block-index scratch */
    float  re_out[HOP * 2];       /* drained output fifo */
    int    qr, qf;                /* fifo read index, fill count */
    int    blend;                 /* blur ramp counter 0..span-1 */
    int    have_prev;
    uint32_t rng;                 /* per-channel xorshift state */
} blur_ch_t;

typedef struct {
    clap_plugin_t     plugin;
    const clap_host_t *host;
    const clap_host_thread_pool_t *thread_pool;

    blur_ch_t ch[MAX_CH];
    float  window[FFT_N];         /* Hann (read-only during process) */
    float  norm[HOP];             /* COLA normalization per phase */
    double span;                  /* blur span in analysis frames (mode 0) */
    double mode;                  /* 0=blur 1=avrg 2=suppress 3=noise 4=chorus 5=scatter */
    double avgw;                  /* neighbor-average width in bins (mode 1) */
    double supn;                  /* loudest partials to zero (mode 2) */
    double noise;                 /* noise blend 0..1 (mode 3) */
    double chor_ampr;             /* chorus max amp ratio >= 1 (mode 4) */
    double chor_frqr;             /* chorus freq jitter in bins (mode 4) */
    double chor_dir;              /* 0=both 1=up 2=down (mode 4) */
    double scat_keep;             /* blocks kept per window (mode 5) */
    double scat_blok;             /* scatter block size in bins (mode 5) */
    int    fill;                  /* samples accumulated toward hop */
} cdp_blur_t;

static void blur_clear(cdp_blur_t *p) {
    p->fill = 0;
    for (int c = 0; c < MAX_CH; ++c) {
        blur_ch_t *b = &p->ch[c];
        memset(b->win, 0, sizeof(b->win));
        memset(b->hop, 0, sizeof(b->hop));
        memset(b->ola, 0, sizeof(b->ola));
        memset(b->re_out, 0, sizeof(b->re_out));
        memset(b->mag_prev, 0, sizeof(b->mag_prev));
        memset(b->ph_prev, 0, sizeof(b->ph_prev));
        b->qr = 0; b->qf = 0;
        b->blend = 0; b->have_prev = 0;
        b->rng = 0x9E3779B9u * (uint32_t)(c + 1);
    }
}

/* per-channel RNG (rand() shared state would race when channels run in parallel) */
static uint32_t ch_rand(blur_ch_t *b) {
    uint32_t x = b->rng;
    x ^= x << 13; x ^= x >> 17; x ^= x << 5;
    b->rng = x ? x : 0x9E3779B9u;
    return b->rng;
}
static float rnd_uni(blur_ch_t *b)     { return (ch_rand(b) >> 8) * (1.0f / 16777216.0f); }  /* [0,1) */
static float rnd_bipolar(blur_ch_t *b) { return 2.0f * rnd_uni(b) - 1.0f; }                  /* [-1,1) */

/* one hop per channel: window -> fft -> make next target frame -> resynth
 * blended frame -> ola -> drain HOP samples into fifo */
static void blur_frame(cdp_blur_t *p, int ci) {
    blur_ch_t *b = &p->ch[ci];
    int span = (int)p->span;

    /* slide window, append new hop */
    memmove(b->win, b->win + HOP, (FFT_N - HOP) * sizeof(float));
    memcpy(b->win + FFT_N - HOP, b->hop, HOP * sizeof(float));

    for (int i = 0; i < FFT_N; ++i) {
        b->re[i] = b->win[i] * p->window[i];
        b->im[i] = 0.0f;
    }
    fft(b->re, b->im, 0);

    if (p->mode == MODE_AVRG) {
        /* polar -> neighbor-averaged magnitudes, phase preserved per bin */
        for (int k = 0; k < BINS; ++k) {
            b->mag_tb[k] = hypotf(b->re[k], b->im[k]);
            b->ph_tb[k]  = atan2f(b->im[k], b->re[k]);
        }
        int w = ((int)p->avgw) / 2;            /* half-width */
        for (int k = 0; k < BINS; ++k) {
            int lo = k - w; if (lo < 0) lo = 0;
            int hi = k + w; if (hi >= BINS) hi = BINS - 1;
            float sum = 0; int cnt = 0;
            for (int j = lo; j <= hi; ++j) { sum += b->mag_tb[j]; cnt++; }
            float mg = sum / cnt;
            b->re[k] = mg * cosf(b->ph_tb[k]);
            b->im[k] = mg * sinf(b->ph_tb[k]);
        }
    } else if (p->mode == MODE_SUPPRESS) {
        /* zero the N loudest partials (CDP blur suppress) */
        int n = (int)p->supn;
        for (int k = 0; k < BINS; ++k)
            b->mag_tb[k] = hypotf(b->re[k], b->im[k]);
        for (int i = 0; i < n; ++i) {
            int best = 0; float bm = 0.0f;
            for (int k = 0; k < BINS; ++k)
                if (b->mag_tb[k] > bm) { bm = b->mag_tb[k]; best = k; }
            if (bm == 0.0f) break;
            b->mag_tb[best] = 0.0f;
            b->re[best] = 0.0f; b->im[best] = 0.0f;
        }
    } else if (p->mode == MODE_NOISE) {
        /* blend white noise into spectrum (CDP blur noise: 0=dry, 1=all noise) */
        float t = (float)p->noise;
        for (int k = 0; k < BINS; ++k) {
            float mg = hypotf(b->re[k], b->im[k]);
            float ph = atan2f(b->im[k], b->re[k]);
            mg = mg * (1.0f - t) + t * rnd_uni(b) * mg;
            b->re[k] = mg * cosf(ph);
            b->im[k] = mg * sinf(ph);
        }
    } else if (p->mode == MODE_CHORUS) {
        /* randomize partial amps and/or freqs by a spread (CDP blur chorus) */
        float ampr = (float)p->chor_ampr;
        float frqr = (float)p->chor_frqr;
        int   dir  = (int)p->chor_dir;
        for (int k = 0; k < BINS; ++k) {
            b->mag_tb[k] = hypotf(b->re[k], b->im[k]);
            b->ph_tb[k]  = atan2f(b->im[k], b->re[k]);
        }
        if (frqr > 0.0f) {                    /* freq scatter: move energy between bins */
            for (int k = 0; k < BINS; ++k) { b->re[k] = 0.0f; b->im[k] = 0.0f; }
            for (int k = 0; k < BINS; ++k) {
                float v = dir == 1 ? rnd_uni(b) : dir == 2 ? -rnd_uni(b) : rnd_bipolar(b);
                int k2 = k + (int)lrintf(v * frqr);
                if (k2 < 0) k2 = 0; if (k2 >= BINS) k2 = BINS - 1;
                b->re[k2] += b->mag_tb[k] * cosf(b->ph_tb[k]);
                b->im[k2] += b->mag_tb[k] * sinf(b->ph_tb[k]);
            }
        } else {
            for (int k = 0; k < BINS; ++k) {
                b->re[k] = b->mag_tb[k] * cosf(b->ph_tb[k]);
                b->im[k] = b->mag_tb[k] * sinf(b->ph_tb[k]);
            }
        }
        if (ampr > 1.0f) {                    /* amp scatter: mag *= ampr^u, u in [-1,1] */
            for (int k = 0; k < BINS; ++k) {
                float m = powf(ampr, rnd_bipolar(b));
                b->re[k] *= m; b->im[k] *= m;
            }
        }
    } else if (p->mode == MODE_SCATTER) {
        /* randomly thin the spectrum (CDP blur scatter; no output normalizing) */
        int blok   = (int)p->scat_blok;
        int nbloks = BINS / blok;
        int keep   = (int)p->scat_keep;
        if (blok < 1) blok = 1;
        if (BINS % blok) nbloks++;
        if (keep > nbloks) keep = nbloks;
        if (keep < nbloks) {
            for (int i = 0; i < nbloks; ++i) b->bloks[i] = i;
            for (int i = 0; i < keep; ++i) {          /* partial Fisher-Yates */
                int j = i + (int)(ch_rand(b) % (uint32_t)(nbloks - i));
                int t = b->bloks[i]; b->bloks[i] = b->bloks[j]; b->bloks[j] = t;
            }
            for (int k = 0; k < BINS; ++k) {
                int bi = k / blok, kept = 0;
                for (int i = 0; i < keep; ++i)
                    if (b->bloks[i] == bi) { kept = 1; break; }
                if (!kept) { b->re[k] = 0.0f; b->im[k] = 0.0f; }
            }
        }
    } else {
        /* blur: ramp prev->next over 'span' frames (CDP do_the_bltr) */
        if (b->blend == 0 || !b->have_prev) {
            for (int k = 0; k < BINS; ++k) {
                float mg = hypotf(b->re[k], b->im[k]);
                float ph = atan2f(b->im[k], b->re[k]);
                if (!b->have_prev) {
                    b->mag_prev[k] = mg;
                    b->ph_prev[k]  = ph;
                }
                b->re[k] = mg;   /* "next" frame polar coords */
                b->im[k] = ph;
            }
            b->have_prev = 1;
        }

        int span = (int)p->span;
        float t = span > 1 ? (float)(b->blend + 1) / (float)span : 1.0f;

        for (int k = 0; k < BINS; ++k) {
            float mg_next = b->re[k], ph_next = b->im[k];
            float d = ph_next - b->ph_prev[k];
            while (d > (float)M_PI)  d -= 2.0f * (float)M_PI;
            while (d < -(float)M_PI) d += 2.0f * (float)M_PI;
            float mg = b->mag_prev[k] + t * (mg_next - b->mag_prev[k]);
            float ph = b->ph_prev[k] + t * d;
            b->re[k] = mg * cosf(ph);
            b->im[k] = mg * sinf(ph);
        }
        b->blend++;
        if (b->blend >= span) {
            b->blend = 0;
            for (int k = 0; k < BINS; ++k) {
                b->mag_prev[k] = hypotf(b->re[k], b->im[k]);
                b->ph_prev[k]  = atan2f(b->im[k], b->re[k]);
            }
        }
    }

    /* conjugate symmetric iFFT: only BINS stored, mirror */
    for (int k = 0; k < BINS; ++k) { /* bins 0..N/2 */
        float r = b->re[k], i = b->im[k];
        b->re[k] = r; b->im[k] = i;
    }
    for (int k = 1; k < FFT_N / 2; ++k) {
        b->re[FFT_N - k] = b->re[k];
        b->im[FFT_N - k] = -b->im[k];
    }
    fft(b->re, b->im, 1);

    /* window again + overlap-add; norm[] folds analysis/synthesis compensation */
    for (int i = 0; i < FFT_N; ++i)
        b->ola[i] += b->re[i] * p->window[i];

    /* drain HOP from ola into fifo, shift ola */
    for (int i = 0; i < HOP; ++i) {
        if (b->qf < HOP * 2) {
            b->re_out[(b->qr + b->qf) & (HOP * 2 - 1)] =
                b->ola[i] / p->norm[i];
            b->qf++;
        }
    }
    memmove(b->ola, b->ola + HOP, (FFT_N - HOP) * sizeof(float));
    memset(b->ola + FFT_N - HOP, 0, HOP * sizeof(float));
}

/* ---------------- CLAP boilerplate ---------------- */

const clap_plugin_descriptor_t s_blur_desc = {
    .clap_version = CLAP_VERSION_INIT,
    .id           = "com.composersdesktop.cdp.blur",
    .name         = "CDP Blur",
    .vendor       = "Composers Desktop Project",
    .url          = "https://composersdesktop.com",
    .manual_url   = "",
    .support_url  = "",
    .version      = "0.2.0",
    .description  = "Spectral processing (CDP blur: blur, avrg, suppress, noise, chorus, scatter).",
    .features =
        (const char *[]){CLAP_PLUGIN_FEATURE_AUDIO_EFFECT, CLAP_PLUGIN_FEATURE_STEREO, NULL},
};

static uint32_t blur_ports_count(const clap_plugin_t *plugin, bool is_input) { return 1; }
static bool blur_ports_get(const clap_plugin_t *plugin, uint32_t index, bool is_input,
                           clap_audio_port_info_t *info) {
    if (index > 0) return false;
    info->id = 0; info->channel_count = 2; info->flags = CLAP_AUDIO_PORT_IS_MAIN;
    info->port_type = CLAP_PORT_STEREO; info->in_place_pair = CLAP_INVALID_ID;
    snprintf(info->name, sizeof(info->name), "Main"); return true;
}
static const clap_plugin_audio_ports_t s_blur_ports = { .count = blur_ports_count, .get = blur_ports_get };

static uint32_t blur_params_count(const clap_plugin_t *plugin) { return NUM_PARAMS; }
static bool blur_params_get_info(const clap_plugin_t *plugin, uint32_t idx, clap_param_info_t *info) {
    memset(info, 0, sizeof(*info));
    if (idx == P_SPAN) {
        info->id = P_SPAN;
        snprintf(info->name, sizeof(info->name), "Blur span (frames)");
        info->min_value = 1.0; info->max_value = 64.0; info->default_value = 8.0;
        info->flags = CLAP_PARAM_IS_AUTOMATABLE;
        return true;
    }
    if (idx == P_MODE) {
        info->id = P_MODE;
        snprintf(info->name, sizeof(info->name), "Mode");
        info->min_value = 0.0; info->max_value = 5.0; info->default_value = 0.0;
        info->flags = CLAP_PARAM_IS_STEPPED | CLAP_PARAM_IS_AUTOMATABLE;
        return true;
    }
    if (idx == P_AVGW) {
        info->id = P_AVGW;
        snprintf(info->name, sizeof(info->name), "Averge width (bins)");
        info->min_value = 1.0; info->max_value = 31.0; info->default_value = 3.0;
        info->flags = CLAP_PARAM_IS_AUTOMATABLE;
        return true;
    }
    if (idx == P_SUPN) {
        info->id = P_SUPN;
        snprintf(info->name, sizeof(info->name), "Suppress count (partials)");
        info->min_value = 0.0; info->max_value = 64.0; info->default_value = 4.0;
        info->flags = CLAP_PARAM_IS_AUTOMATABLE;
        return true;
    }
    if (idx == P_NOISE) {
        info->id = P_NOISE;
        snprintf(info->name, sizeof(info->name), "Noise level");
        info->min_value = 0.0; info->max_value = 1.0; info->default_value = 0.5;
        info->flags = CLAP_PARAM_IS_AUTOMATABLE;
        return true;
    }
    if (idx == P_CHOR_AMPR) {
        info->id = P_CHOR_AMPR;
        snprintf(info->name, sizeof(info->name), "Chorus amp spread");
        info->min_value = 1.0; info->max_value = 16.0; info->default_value = 2.0;
        info->flags = CLAP_PARAM_IS_AUTOMATABLE;
        return true;
    }
    if (idx == P_CHOR_FRQR) {
        info->id = P_CHOR_FRQR;
        snprintf(info->name, sizeof(info->name), "Chorus freq jitter (bins)");
        info->min_value = 0.0; info->max_value = 32.0; info->default_value = 4.0;
        info->flags = CLAP_PARAM_IS_AUTOMATABLE;
        return true;
    }
    if (idx == P_CHOR_DIR) {
        info->id = P_CHOR_DIR;
        snprintf(info->name, sizeof(info->name), "Chorus direction");
        info->min_value = 0.0; info->max_value = 2.0; info->default_value = 0.0;
        info->flags = CLAP_PARAM_IS_STEPPED | CLAP_PARAM_IS_AUTOMATABLE;
        return true;
    }
    if (idx == P_SCAT_KEEP) {
        info->id = P_SCAT_KEEP;
        snprintf(info->name, sizeof(info->name), "Scatter keep (blocks)");
        info->min_value = 1.0; info->max_value = 64.0; info->default_value = 4.0;
        info->flags = CLAP_PARAM_IS_AUTOMATABLE;
        return true;
    }
    if (idx == P_SCAT_BLOK) {
        info->id = P_SCAT_BLOK;
        snprintf(info->name, sizeof(info->name), "Scatter block size (bins)");
        info->min_value = 1.0; info->max_value = 64.0; info->default_value = 8.0;
        info->flags = CLAP_PARAM_IS_AUTOMATABLE;
        return true;
    }
    return false;
}
static bool blur_params_get_value(const clap_plugin_t *plugin, clap_id id, double *v) {
    const cdp_blur_t *p = plugin->plugin_data;
    if (id == P_SPAN) { *v = p->span; return true; }
    if (id == P_MODE) { *v = p->mode; return true; }
    if (id == P_AVGW) { *v = p->avgw; return true; }
    if (id == P_SUPN) { *v = p->supn; return true; }
    if (id == P_NOISE) { *v = p->noise; return true; }
    if (id == P_CHOR_AMPR) { *v = p->chor_ampr; return true; }
    if (id == P_CHOR_FRQR) { *v = p->chor_frqr; return true; }
    if (id == P_CHOR_DIR) { *v = p->chor_dir; return true; }
    if (id == P_SCAT_KEEP) { *v = p->scat_keep; return true; }
    if (id == P_SCAT_BLOK) { *v = p->scat_blok; return true; }
    return false;
}
static bool blur_params_value_to_text(const clap_plugin_t *plugin, clap_id id, double v,
                                      char *display, uint32_t size) {
    if (id == P_MODE) {
        const char *names[] = {"Blur", "Avrg", "Suppress", "Noise", "Chorus", "Scatter"};
        int m = (int)v; if (m < 0) m = 0; if (m > 5) m = 5;
        snprintf(display, size, "%s", names[m]);
        return true;
    }
    if (id == P_CHOR_DIR) {
        const char *names[] = {"Both", "Up", "Down"};
        int m = (int)v; if (m < 0) m = 0; if (m > 2) m = 2;
        snprintf(display, size, "%s", names[m]);
        return true;
    }
    if (id == P_NOISE || id == P_CHOR_AMPR) {
        snprintf(display, size, "%.2f", v);
        return true;
    }
    if (id == P_SPAN || id == P_AVGW || id == P_SUPN ||
        id == P_CHOR_FRQR || id == P_SCAT_KEEP || id == P_SCAT_BLOK) {
        snprintf(display, size, "%.0f", v);
        return true;
    }
    return false;
}
static bool blur_params_text_to_value(const clap_plugin_t *plugin, clap_id id, const char *display, double *v) {
    if (id == P_MODE) {
        if (!strcmp(display, "Avrg"))     { *v = MODE_AVRG; return true; }
        if (!strcmp(display, "Blur"))     { *v = MODE_BLUR; return true; }
        if (!strcmp(display, "Suppress")) { *v = MODE_SUPPRESS; return true; }
        if (!strcmp(display, "Noise"))    { *v = MODE_NOISE; return true; }
        if (!strcmp(display, "Chorus"))   { *v = MODE_CHORUS; return true; }
        if (!strcmp(display, "Scatter"))  { *v = MODE_SCATTER; return true; }
        return false;
    }
    if (id == P_CHOR_DIR) {
        if (!strcmp(display, "Both")) { *v = 0.0; return true; }
        if (!strcmp(display, "Up"))   { *v = 1.0; return true; }
        if (!strcmp(display, "Down")) { *v = 2.0; return true; }
        return false;
    }
    char *end = NULL; double t = strtod(display, &end);
    if (end == display) return false;
    if (id == P_SPAN) {
        if (t < 1.0) t = 1.0; if (t > 64.0) t = 64.0;
        *v = t; return true;
    }
    if (id == P_AVGW) {
        if (t < 1.0) t = 1.0; if (t > 31.0) t = 31.0;
        *v = t; return true;
    }
    if (id == P_SUPN) {
        if (t < 0.0) t = 0.0; if (t > 64.0) t = 64.0;
        *v = t; return true;
    }
    if (id == P_NOISE) {
        if (t < 0.0) t = 0.0; if (t > 1.0) t = 1.0;
        *v = t; return true;
    }
    if (id == P_CHOR_AMPR) {
        if (t < 1.0) t = 1.0; if (t > 16.0) t = 16.0;
        *v = t; return true;
    }
    if (id == P_CHOR_FRQR || id == P_SCAT_KEEP || id == P_SCAT_BLOK) {
        if (t < 0.0) t = 0.0; if (t > 64.0) t = 64.0;
        *v = t; return true;
    }
    return false;
}
static void blur_apply_event(cdp_blur_t *p, const clap_event_header_t *hdr) {
    if (hdr->space_id == CLAP_CORE_EVENT_SPACE_ID && hdr->type == CLAP_EVENT_PARAM_VALUE) {
        const clap_event_param_value_t *ev = (const clap_event_param_value_t *)hdr;
        if (ev->param_id == P_SPAN) {
            double t = ev->value; if (t < 1.0) t = 1.0; if (t > 64.0) t = 64.0;
            p->span = t;
        }
        if (ev->param_id == P_MODE) {
            int m = (int)ev->value; if (m < 0) m = 0; if (m > 5) m = 5;
            if ((int)p->mode != m) {
                for (int c = 0; c < MAX_CH; ++c) { p->ch[c].blend = 0; p->ch[c].have_prev = 0; }
            }
            p->mode = (double)m;
        }
        if (ev->param_id == P_AVGW) {
            double t = ev->value; if (t < 1.0) t = 1.0; if (t > 31.0) t = 31.0;
            p->avgw = t;
        }
        if (ev->param_id == P_SUPN) {
            double t = ev->value; if (t < 0.0) t = 0.0; if (t > 64.0) t = 64.0;
            p->supn = t;
        }
        if (ev->param_id == P_NOISE) {
            double t = ev->value; if (t < 0.0) t = 0.0; if (t > 1.0) t = 1.0;
            p->noise = t;
        }
        if (ev->param_id == P_CHOR_AMPR) {
            double t = ev->value; if (t < 1.0) t = 1.0; if (t > 16.0) t = 16.0;
            p->chor_ampr = t;
        }
        if (ev->param_id == P_CHOR_FRQR) {
            double t = ev->value; if (t < 0.0) t = 0.0; if (t > 32.0) t = 32.0;
            p->chor_frqr = t;
        }
        if (ev->param_id == P_CHOR_DIR) {
            int m = (int)ev->value; if (m < 0) m = 0; if (m > 2) m = 2;
            p->chor_dir = (double)m;
        }
        if (ev->param_id == P_SCAT_KEEP) {
            double t = ev->value; if (t < 1.0) t = 1.0; if (t > 64.0) t = 64.0;
            p->scat_keep = t;
        }
        if (ev->param_id == P_SCAT_BLOK) {
            double t = ev->value; if (t < 1.0) t = 1.0; if (t > 64.0) t = 64.0;
            p->scat_blok = t;
        }
    }
}
static void blur_params_flush(const clap_plugin_t *plugin, const clap_input_events_t *in,
                              const clap_output_events_t *out) {
    cdp_blur_t *p = plugin->plugin_data;
    const uint32_t nev = in->size(in);
    for (uint32_t i = 0; i < nev; ++i) blur_apply_event(p, in->get(in, i));
}
static const clap_plugin_params_t s_blur_params = {
    .count = blur_params_count, .get_info = blur_params_get_info, .get_value = blur_params_get_value,
    .value_to_text = blur_params_value_to_text, .text_to_value = blur_params_text_to_value,
    .flush = blur_params_flush,
};

static bool blur_state_save(const clap_plugin_t *plugin, const clap_ostream_t *stream) {
    const cdp_blur_t *p = plugin->plugin_data;
    double st[NUM_PARAMS] = {p->span, p->mode, p->avgw, p->supn, p->noise,
                             p->chor_ampr, p->chor_frqr, p->chor_dir,
                             p->scat_keep, p->scat_blok};
    return stream->write(stream, st, sizeof(st)) == sizeof(st);
}
static bool blur_state_load(const clap_plugin_t *plugin, const clap_istream_t *stream) {
    cdp_blur_t *p = plugin->plugin_data; double st[NUM_PARAMS];
    if (stream->read(stream, st, sizeof(st)) != sizeof(st)) return false;
    if (st[0] < 1.0) st[0] = 1.0; if (st[0] > 64.0) st[0] = 64.0;
    if (st[1] < 0.0) st[1] = 0.0; if (st[1] > 5.0) st[1] = 5.0;
    if (st[2] < 1.0) st[2] = 1.0; if (st[2] > 31.0) st[2] = 31.0;
    if (st[3] < 0.0) st[3] = 0.0; if (st[3] > 64.0) st[3] = 64.0;
    if (st[4] < 0.0) st[4] = 0.0; if (st[4] > 1.0) st[4] = 1.0;
    if (st[5] < 1.0) st[5] = 1.0; if (st[5] > 16.0) st[5] = 16.0;
    if (st[6] < 0.0) st[6] = 0.0; if (st[6] > 32.0) st[6] = 32.0;
    if (st[7] < 0.0) st[7] = 0.0; if (st[7] > 2.0) st[7] = 2.0;
    if (st[8] < 1.0) st[8] = 1.0; if (st[8] > 64.0) st[8] = 64.0;
    if (st[9] < 1.0) st[9] = 1.0; if (st[9] > 64.0) st[9] = 64.0;
    p->span = st[0]; p->mode = st[1]; p->avgw = st[2]; p->supn = st[3]; p->noise = st[4];
    p->chor_ampr = st[5]; p->chor_frqr = st[6]; p->chor_dir = st[7];
    p->scat_keep = st[8]; p->scat_blok = st[9];
    return true;
}
static const clap_plugin_state_t s_blur_state = { .save = blur_state_save, .load = blur_state_load };

static bool blur_init(const clap_plugin_t *plugin) {
    cdp_blur_t *p = plugin->plugin_data;
    p->thread_pool = (const clap_host_thread_pool_t *)
        p->host->get_extension(p->host, CLAP_EXT_THREAD_POOL);
    return true;
}
static void blur_destroy(const clap_plugin_t *plugin) { free(plugin->plugin_data); }
static bool blur_activate(const clap_plugin_t *plugin, double sr, uint32_t minf, uint32_t maxf) {
    cdp_blur_t *p = plugin->plugin_data;
    fft_init();
    for (int i = 0; i < FFT_N; ++i)
        p->window[i] = 0.5f - 0.5f * cosf(2.0f * (float)M_PI * i / FFT_N);
    /* exact per-phase COLA normalization: sum of w^2 over the 4 overlaps */
    for (int ph = 0; ph < HOP; ++ph) {
        float s = 0;
        for (int k = 0; k < 4; ++k) {
            float w = p->window[ph + k * HOP];
            s += w * w;
        }
        p->norm[ph] = s;
    }
    blur_clear(p);
    return true;
}
static void blur_deactivate(const clap_plugin_t *plugin) {}
static bool blur_start_processing(const clap_plugin_t *plugin) { return true; }
static void blur_stop_processing(const clap_plugin_t *plugin) {}
static void blur_reset(const clap_plugin_t *plugin) {
    cdp_blur_t *p = plugin->plugin_data;
    blur_clear(p);
}

static float read_samp(const clap_audio_buffer_t *b, uint32_t ch, uint32_t i) {
    return b->data32 ? b->data32[ch][i] : (float)b->data64[ch][i];
}
static void write_samp(clap_audio_buffer_t *b, uint32_t ch, uint32_t i, float v) {
    if (b->data32) b->data32[ch][i] = v; else b->data64[ch][i] = v;
}

static clap_process_status blur_process(const clap_plugin_t *plugin, const clap_process_t *process) {
    cdp_blur_t *p = plugin->plugin_data;
    const clap_input_events_t *evs = process->in_events;
    uint32_t nev = evs ? evs->size(evs) : 0;
    for (uint32_t i = 0; i < nev; ++i) blur_apply_event(p, evs->get(evs, i));

    if (process->audio_outputs_count == 0) return CLAP_PROCESS_CONTINUE;
    const clap_audio_buffer_t *inb = process->audio_inputs_count > 0 ? &process->audio_inputs[0] : NULL;
    clap_audio_buffer_t *outb = &process->audio_outputs[0];
    uint32_t nch = outb->channel_count;
    if (inb && inb->channel_count < nch) nch = inb->channel_count;
    if (nch > MAX_CH) nch = MAX_CH;

    for (uint32_t i = 0; i < process->frames_count; ++i) {
        for (uint32_t c = 0; c < nch; ++c) {
            blur_ch_t *b = &p->ch[c];
            b->hop[p->fill] = inb ? read_samp(inb, c, i) : 0.0f;
            float out;
            if (b->qf > 0) {
                out = b->re_out[b->qr];
                b->qr = (b->qr + 1) & (HOP * 2 - 1);
                b->qf--;
            } else
                out = 0.0f;   /* initial latency HOP */
            write_samp(outb, c, i, out);
        }
        if (++p->fill == HOP) {
            p->fill = 0;
            /* per-channel frames are independent → host thread pool if offered */
            int did_parallel = 0;
            if (nch > 1 && p->thread_pool && p->thread_pool->request_exec)
                did_parallel = p->thread_pool->request_exec(p->host, nch);
            if (!did_parallel)
                for (uint32_t c = 0; c < nch; ++c)
                    blur_frame(p, c);
        }
    }
    return CLAP_PROCESS_CONTINUE;
}

static void blur_thread_pool_exec(const clap_plugin_t *plugin, uint32_t task_index) {
    cdp_blur_t *p = plugin->plugin_data;
    if (task_index < MAX_CH)
        blur_frame(p, (int)task_index);
}
static const clap_plugin_thread_pool_t s_blur_thread_pool = { .exec = blur_thread_pool_exec };

static const void *blur_get_extension(const clap_plugin_t *plugin, const char *id) {
    if (!strcmp(id, CLAP_EXT_AUDIO_PORTS)) return &s_blur_ports;
    if (!strcmp(id, CLAP_EXT_PARAMS))      return &s_blur_params;
    if (!strcmp(id, CLAP_EXT_STATE))       return &s_blur_state;
    if (!strcmp(id, CLAP_EXT_THREAD_POOL)) return &s_blur_thread_pool;
    return NULL;
}
static void blur_on_main_thread(const clap_plugin_t *plugin) {}

clap_plugin_t *blur_create(const clap_host_t *host) {
    cdp_blur_t *p = calloc(1, sizeof(*p));
    p->host = host; p->span = 8.0; p->mode = MODE_BLUR; p->avgw = 3.0; p->supn = 4.0; p->noise = 0.5;
    p->chor_ampr = 2.0; p->chor_frqr = 4.0; p->chor_dir = 0.0;
    p->scat_keep = 4.0; p->scat_blok = 8.0;
    p->plugin.desc = &s_blur_desc; p->plugin.plugin_data = p;
    p->plugin.init = blur_init; p->plugin.destroy = blur_destroy;
    p->plugin.activate = blur_activate; p->plugin.deactivate = blur_deactivate;
    p->plugin.start_processing = blur_start_processing;
    p->plugin.stop_processing = blur_stop_processing;
    p->plugin.reset = blur_reset; p->plugin.process = blur_process;
    p->plugin.get_extension = blur_get_extension;
    p->plugin.on_main_thread = blur_on_main_thread;
    return &p->plugin;
}
