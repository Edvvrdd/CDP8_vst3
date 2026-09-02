/* CDP Distort — wavecycle and sample distortion as a CLAP effect.
 *
 * Realtime ports from CDP 'distort' (dev/distort/):
 *
 * Mode 0: Overload (Noise)   - Clips peaks and replaces with noise (overload 1)
 * Mode 1: Overload (Cosine)  - Clips peaks and replaces with cosine modulation (overload 2)
 * Mode 2: Reform Square      - Converts each wavecycle to square shape (reform 2)
 * Mode 3: Reform Triangle    - Converts each wavecycle to triangular ramp (reform 4)
 * Mode 4: Reform Invert      - Inverts the waveform profile of each half-cycle (reform 5)
 * Mode 5: Reform Sine        - Replaces half-cycles with synthesized sine arcs (reform 7)
 * Mode 6: Exaggerate         - Non-linear wave contour expansion (reform 8)
 * Mode 7: Multiply           - Multiplies wavecycle frequency (multiply)
 *
 * All modes work in the sample and wavecycle domain without external files.
 */

#include <clap/clap.h>

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MAX_CH 2
#define MAX_CYCLE_LEN 4096

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

enum {
    P_CLIP = 0,
    P_DEPTH = 1,
    P_FREQ = 2,
    P_EXAGG = 3,
    P_MULT = 4,
    P_ON_0 = 5,
    NUM_PARAMS = P_ON_0 + 8
};

static const char *s_mode_names[8] = {
    "Overload Noise",
    "Overload Sine",
    "Reform Square",
    "Reform Triangle",
    "Reform Invert",
    "Reform Sine",
    "Reform Exaggerate",
    "Multiply"
};

typedef struct {
    uint32_t rng;
    int phase;                    /* 1 for positive, -1 for negative */
    float cycle_buf[MAX_CYCLE_LEN];
    float temp[MAX_CYCLE_LEN];    /* per-mode output scratch (off the stack) */
    int cycle_pos;
    float cycle_peak;
    float cos_phase;
} dist_ch_t;

typedef struct {
    clap_plugin_t plugin;
    const clap_host_t *host;
    dist_ch_t ch[MAX_CH];
    double clip;
    double depth;
    double freq;
    double exagg;
    double mult;
    double on[8];
    double sample_rate;
} cdp_dist_t;

static uint32_t ch_rand(dist_ch_t *c) {
    uint32_t x = c->rng;
    x ^= x << 13; x ^= x >> 17; x ^= x << 5;
    c->rng = x ? x : 0x9E3779B9u;
    return c->rng;
}

static float rnd_uni(dist_ch_t *c) {
    return (ch_rand(c) >> 8) * (1.0f / 16777216.0f);
}

/* Overload 1: Noise */
static float proc_overload_noise(cdp_dist_t *p, dist_ch_t *c, float s) {
    float clip = (float)p->clip;
    if (s >= clip) {
        s = clip * (1.0f - rnd_uni(c) * (float)p->depth);
    } else if (s <= -clip) {
        s = -clip * (1.0f - rnd_uni(c) * (float)p->depth);
    }
    return s / clip;
}

/* Overload 2: Cosine Modulation */
static float proc_overload_sine(cdp_dist_t *p, dist_ch_t *c, float s) {
    float clip = (float)p->clip;
    float inc = (float)(2.0 * M_PI * p->freq / p->sample_rate);
    c->cos_phase = fmodf(c->cos_phase + inc, (float)(2.0 * M_PI));
    float mod = (cosf(c->cos_phase) - 1.0f) * 0.5f; /* -1 to 0 */

    if (s >= clip) {
        s = clip * (1.0f + mod * (float)p->depth);
    } else if (s <= -clip) {
        s = -clip * (1.0f + mod * (float)p->depth);
    }
    return s / clip;
}

/* Process completed half-cycle for wavecycle modes */
static void process_half_cycle(cdp_dist_t *p, dist_ch_t *c, float *out_buf, int len, int mode) {
    if (len <= 0) return;
    float peak = c->cycle_peak;

    switch (mode) {
    case 2: /* Reform Square */
        for (int i = 0; i < len; ++i) {
            out_buf[i] = peak;
        }
        break;

    case 3: /* Reform Triangle */
        if (len <= 2) {
            for (int i = 0; i < len; ++i) out_buf[i] = peak;
        } else {
            int mid = len / 2;
            for (int i = 0; i < len; ++i) {
                if (i <= mid) {
                    out_buf[i] = (mid > 0) ? peak * ((float)i / (float)mid) : peak;
                } else {
                    out_buf[i] = peak * ((float)(len - i) / (float)(len - mid));
                }
            }
        }
        break;

    case 4: /* Reform Invert */
        for (int i = 0; i < len; ++i) {
            out_buf[i] = peak - c->cycle_buf[i];
        }
        break;

    case 5: /* Reform Sine */
        for (int i = 0; i < len; ++i) {
            float ph = (float)M_PI * ((float)i / (float)len);
            out_buf[i] = (c->phase > 0) ? peak * sinf(ph) : -peak * sinf(ph);
        }
        break;

    case 6: /* Exaggerate contour */
        for (int i = 0; i < len; ++i) {
            float val = c->cycle_buf[i];
            float norm = (peak > 0.00001f || peak < -0.00001f) ? val / peak : 0.0f;
            if (norm < 0.0f) norm = 0.0f;
            if (norm > 1.0f) norm = 1.0f;
            out_buf[i] = powf(norm, (float)p->exagg) * peak;
        }
        break;

    case 7: /* Multiply */
        {
            int mult = (int)p->mult;
            if (mult < 1) mult = 1;
            for (int i = 0; i < len; ++i) {
                int src_idx = (i * mult) % len;
                out_buf[i] = c->cycle_buf[src_idx];
            }
        }
        break;

    default:
        for (int i = 0; i < len; ++i) out_buf[i] = c->cycle_buf[i];
        break;
    }
}

/* buffer access: hosts may deliver either data32 or data64 */
static float read_samp(const clap_audio_buffer_t *b, uint32_t ch, uint32_t i) {
    return b->data32 ? b->data32[ch][i] : (float)b->data64[ch][i];
}
static void write_samp(clap_audio_buffer_t *b, uint32_t ch, uint32_t i, float v) {
    if (b->data32) b->data32[ch][i] = v; else b->data64[ch][i] = v;
}

/* Single channel frame processing with wavecycle tracker */
static void process_channel_samples(cdp_dist_t *p, int ch_idx, const clap_audio_buffer_t *inb,
                                    clap_audio_buffer_t *outb, uint32_t frames) {
    dist_ch_t *c = &p->ch[ch_idx];
    uint32_t ch = (uint32_t)ch_idx;
    int n_on = 0;
    for (int m = 0; m < 8; ++m) {
        if (p->on[m] >= 0.5) n_on++;
    }

    if (n_on == 0) {
        /* Passthrough */
        for (uint32_t i = 0; i < frames; ++i)
            write_samp(outb, ch, i, inb ? read_samp(inb, ch, i) : 0.0f);
        return;
    }

    /* Process sample-by-sample */
    for (uint32_t i = 0; i < frames; ++i) {
        float in_samp = inb ? read_samp(inb, ch, i) : 0.0f;
        float acc = 0.0f;

        /* Mode 0: Overload Noise */
        if (p->on[0] >= 0.5) {
            acc += proc_overload_noise(p, c, in_samp);
        }

        /* Mode 1: Overload Sine */
        if (p->on[1] >= 0.5) {
            acc += proc_overload_sine(p, c, in_samp);
        }

        /* Zero-crossing detection for wavecycle modes (2..7) */
        int cur_phase = (in_samp >= 0.0f) ? 1 : -1;
        if (cur_phase != c->phase && c->cycle_pos > 0) {
            /* Zero-crossing event! Trigger active wavecycle modes */
            for (int m = 2; m < 8; ++m) {
                if (p->on[m] >= 0.5) {
                    process_half_cycle(p, c, c->temp, c->cycle_pos, m);
                    /* Replace the most recent processed slice in out buffer if safe */
                    int start = (int)i - c->cycle_pos;
                    if (start >= 0) {
                        for (int k = 0; k < c->cycle_pos; ++k) {
                            float prev = read_samp(outb, ch, (uint32_t)(start + k));
                            write_samp(outb, ch, (uint32_t)(start + k),
                                       (prev * (n_on - 1) + c->temp[k]) / (float)n_on);
                        }
                    }
                }
            }
            c->phase = cur_phase;
            c->cycle_pos = 0;
            c->cycle_peak = 0.0f;
        }

        /* Buffer sample into wavecycle tracker */
        if (c->cycle_pos < MAX_CYCLE_LEN) {
            c->cycle_buf[c->cycle_pos++] = in_samp;
            if (fabsf(in_samp) > fabsf(c->cycle_peak)) {
                c->cycle_peak = in_samp;
            }
        }

        /* For sample modes, write accumulator scaled by number of active modes */
        float out_val = (p->on[0] >= 0.5 || p->on[1] >= 0.5) ? acc / (float)n_on : in_samp;
        write_samp(outb, ch, i, out_val);
    }
}

/* ---------- Descriptor ---------- */

const clap_plugin_descriptor_t s_distort_desc = {
    .clap_version = CLAP_VERSION_INIT,
    .id           = "com.composersdesktop.cdp.distort",
    .name         = "CDP Distort",
    .vendor       = "Composers Desktop Project",
    .url          = "https://composersdesktop.com",
    .manual_url   = "",
    .support_url  = "",
    .version      = "0.2.0",
    .description  = "CDP Distort collection (Overload, Reform, Multiply).",
    .features =
        (const char *[]){CLAP_PLUGIN_FEATURE_AUDIO_EFFECT, CLAP_PLUGIN_FEATURE_STEREO, CLAP_PLUGIN_FEATURE_DISTORTION, NULL},
};

/* ---------- Ports ---------- */

static uint32_t dist_ports_count(const clap_plugin_t *plugin, bool is_input) { return 1; }
static bool dist_ports_get(const clap_plugin_t *plugin, uint32_t index, bool is_input,
                           clap_audio_port_info_t *info) {
    if (index > 0) return false;
    info->id = 0; info->channel_count = 2; info->flags = CLAP_AUDIO_PORT_IS_MAIN;
    info->port_type = CLAP_PORT_STEREO; info->in_place_pair = CLAP_INVALID_ID;
    snprintf(info->name, sizeof(info->name), "Main"); return true;
}
static const clap_plugin_audio_ports_t s_dist_ports = { .count = dist_ports_count, .get = dist_ports_get };

/* ---------- Params ---------- */

static uint32_t dist_params_count(const clap_plugin_t *plugin) { return NUM_PARAMS; }
static bool dist_params_get_info(const clap_plugin_t *plugin, uint32_t idx, clap_param_info_t *info) {
    memset(info, 0, sizeof(*info));
    if (idx == P_CLIP) {
        info->id = P_CLIP;
        snprintf(info->name, sizeof(info->name), "Clip Level");
        info->min_value = 0.01; info->max_value = 1.0; info->default_value = 0.5;
        info->flags = CLAP_PARAM_IS_AUTOMATABLE;
        return true;
    }
    if (idx == P_DEPTH) {
        info->id = P_DEPTH;
        snprintf(info->name, sizeof(info->name), "Depth");
        info->min_value = 0.0; info->max_value = 1.0; info->default_value = 0.5;
        info->flags = CLAP_PARAM_IS_AUTOMATABLE;
        return true;
    }
    if (idx == P_FREQ) {
        info->id = P_FREQ;
        snprintf(info->name, sizeof(info->name), "Sine Mod Freq (Hz)");
        info->min_value = 1.0; info->max_value = 2000.0; info->default_value = 440.0;
        info->flags = CLAP_PARAM_IS_AUTOMATABLE;
        return true;
    }
    if (idx == P_EXAGG) {
        info->id = P_EXAGG;
        snprintf(info->name, sizeof(info->name), "Exaggeration Power");
        info->min_value = 0.1; info->max_value = 8.0; info->default_value = 2.0;
        info->flags = CLAP_PARAM_IS_AUTOMATABLE;
        return true;
    }
    if (idx == P_MULT) {
        info->id = P_MULT;
        snprintf(info->name, sizeof(info->name), "Frequency Multiplier");
        info->min_value = 2.0; info->max_value = 16.0; info->default_value = 2.0;
        info->flags = CLAP_PARAM_IS_AUTOMATABLE;
        return true;
    }
    if (idx >= P_ON_0 && idx < NUM_PARAMS) {
        info->id = (clap_id)idx;
        snprintf(info->name, sizeof(info->name), "%s on/off", s_mode_names[idx - P_ON_0]);
        info->min_value = 0.0; info->max_value = 1.0; info->default_value = (idx == P_ON_0) ? 1.0 : 0.0;
        info->flags = CLAP_PARAM_IS_STEPPED | CLAP_PARAM_IS_AUTOMATABLE;
        return true;
    }
    return false;
}

static bool dist_params_get_value(const clap_plugin_t *plugin, clap_id id, double *v) {
    const cdp_dist_t *p = plugin->plugin_data;
    if (id == P_CLIP)  { *v = p->clip; return true; }
    if (id == P_DEPTH) { *v = p->depth; return true; }
    if (id == P_FREQ)  { *v = p->freq; return true; }
    if (id == P_EXAGG) { *v = p->exagg; return true; }
    if (id == P_MULT)  { *v = p->mult; return true; }
    if (id >= P_ON_0 && id < NUM_PARAMS) { *v = p->on[id - P_ON_0]; return true; }
    return false;
}

static bool dist_params_value_to_text(const clap_plugin_t *plugin, clap_id id, double v,
                                      char *display, uint32_t size) {
    if (id >= P_ON_0 && id < NUM_PARAMS) {
        snprintf(display, size, "%s", v >= 0.5 ? "On" : "Off");
        return true;
    }
    if (id == P_CLIP || id == P_DEPTH) {
        snprintf(display, size, "%.2f", v);
        return true;
    }
    if (id == P_FREQ) {
        snprintf(display, size, "%.1f Hz", v);
        return true;
    }
    if (id == P_EXAGG || id == P_MULT) {
        snprintf(display, size, "%.1f", v);
        return true;
    }
    return false;
}

static bool dist_params_text_to_value(const clap_plugin_t *plugin, clap_id id, const char *display, double *v) {
    if (id >= P_ON_0 && id < NUM_PARAMS) {
        if (!strcmp(display, "On"))  { *v = 1.0; return true; }
        if (!strcmp(display, "Off")) { *v = 0.0; return true; }
        return false;
    }
    char *end = NULL; double t = strtod(display, &end);
    if (end == display) return false;
    if (id == P_CLIP) {
        if (t < 0.01) t = 0.01; if (t > 1.0) t = 1.0;
        *v = t; return true;
    }
    if (id == P_DEPTH) {
        if (t < 0.0) t = 0.0; if (t > 1.0) t = 1.0;
        *v = t; return true;
    }
    if (id == P_FREQ) {
        if (t < 1.0) t = 1.0; if (t > 2000.0) t = 2000.0;
        *v = t; return true;
    }
    if (id == P_EXAGG) {
        if (t < 0.1) t = 0.1; if (t > 8.0) t = 8.0;
        *v = t; return true;
    }
    if (id == P_MULT) {
        if (t < 2.0) t = 2.0; if (t > 16.0) t = 16.0;
        *v = t; return true;
    }
    return false;
}

static void dist_apply_event(cdp_dist_t *p, const clap_event_header_t *hdr) {
    if (hdr->space_id == CLAP_CORE_EVENT_SPACE_ID && hdr->type == CLAP_EVENT_PARAM_VALUE) {
        const clap_event_param_value_t *ev = (const clap_event_param_value_t *)hdr;
        if (ev->param_id == P_CLIP)  p->clip = ev->value;
        if (ev->param_id == P_DEPTH) p->depth = ev->value;
        if (ev->param_id == P_FREQ)  p->freq = ev->value;
        if (ev->param_id == P_EXAGG) p->exagg = ev->value;
        if (ev->param_id == P_MULT)  p->mult = ev->value;
        if (ev->param_id >= P_ON_0 && ev->param_id < NUM_PARAMS) {
            p->on[ev->param_id - P_ON_0] = ev->value >= 0.5 ? 1.0 : 0.0;
        }
    }
}

static void dist_params_flush(const clap_plugin_t *plugin, const clap_input_events_t *in,
                              const clap_output_events_t *out) {
    cdp_dist_t *p = plugin->plugin_data;
    const uint32_t nev = in->size(in);
    for (uint32_t i = 0; i < nev; ++i) dist_apply_event(p, in->get(in, i));
}

static const clap_plugin_params_t s_dist_params = {
    .count = dist_params_count, .get_info = dist_params_get_info, .get_value = dist_params_get_value,
    .value_to_text = dist_params_value_to_text, .text_to_value = dist_params_text_to_value,
    .flush = dist_params_flush,
};

/* ---------- State ---------- */

static bool dist_state_save(const clap_plugin_t *plugin, const clap_ostream_t *stream) {
    const cdp_dist_t *p = plugin->plugin_data;
    double st[NUM_PARAMS] = {
        p->clip, p->depth, p->freq, p->exagg, p->mult,
        p->on[0], p->on[1], p->on[2], p->on[3], p->on[4], p->on[5], p->on[6], p->on[7]
    };
    return stream->write(stream, st, sizeof(st)) == sizeof(st);
}

static bool dist_state_load(const clap_plugin_t *plugin, const clap_istream_t *stream) {
    cdp_dist_t *p = plugin->plugin_data; double st[NUM_PARAMS];
    if (stream->read(stream, st, sizeof(st)) != sizeof(st)) return false;
    p->clip  = st[0];
    p->depth = st[1];
    p->freq  = st[2];
    p->exagg = st[3];
    p->mult  = st[4];
    for (int m = 0; m < 8; ++m) {
        p->on[m] = st[P_ON_0 + m] >= 0.5 ? 1.0 : 0.0;
    }
    return true;
}

static const clap_plugin_state_t s_dist_state = { .save = dist_state_save, .load = dist_state_load };

/* ---------- Lifecycle ---------- */

static bool dist_init(const clap_plugin_t *plugin) { return true; }
static void dist_destroy(const clap_plugin_t *plugin) { free(plugin->plugin_data); }

static bool dist_activate(const clap_plugin_t *plugin, double sr, uint32_t minf, uint32_t maxf) {
    cdp_dist_t *p = plugin->plugin_data;
    p->sample_rate = sr > 1000.0 ? sr : 44100.0;
    for (int c = 0; c < MAX_CH; ++c) {
        p->ch[c].rng = 0x9E3779B9u * (uint32_t)(c + 3);
        p->ch[c].phase = 1;
        p->ch[c].cycle_pos = 0;
        p->ch[c].cycle_peak = 0.0f;
        p->ch[c].cos_phase = 0.0f;
    }
    return true;
}

static void dist_deactivate(const clap_plugin_t *plugin) {}
static bool dist_start_processing(const clap_plugin_t *plugin) { return true; }
static void dist_stop_processing(const clap_plugin_t *plugin) {}
static void dist_reset(const clap_plugin_t *plugin) {
    cdp_dist_t *p = plugin->plugin_data;
    for (int c = 0; c < MAX_CH; ++c) {
        p->ch[c].cycle_pos = 0;
        p->ch[c].cycle_peak = 0.0f;
    }
}

static clap_process_status dist_process(const clap_plugin_t *plugin, const clap_process_t *process) {
    cdp_dist_t *p = plugin->plugin_data;
    const clap_input_events_t *evs = process->in_events;
    uint32_t nev = evs ? evs->size(evs) : 0;
    for (uint32_t i = 0; i < nev; ++i) dist_apply_event(p, evs->get(evs, i));

    if (process->audio_outputs_count == 0) return CLAP_PROCESS_CONTINUE;
    const clap_audio_buffer_t *inb = process->audio_inputs_count > 0 ? &process->audio_inputs[0] : NULL;
    clap_audio_buffer_t *outb = &process->audio_outputs[0];
    uint32_t nch = outb->channel_count;
    if (inb && inb->channel_count < nch) nch = inb->channel_count;
    if (nch > MAX_CH) nch = MAX_CH;

    for (uint32_t c = 0; c < nch; ++c)
        process_channel_samples(p, (int)c, inb, outb, process->frames_count);

    return CLAP_PROCESS_CONTINUE;
}

static const void *dist_get_extension(const clap_plugin_t *plugin, const char *id) {
    if (!strcmp(id, CLAP_EXT_AUDIO_PORTS)) return &s_dist_ports;
    if (!strcmp(id, CLAP_EXT_PARAMS))      return &s_dist_params;
    if (!strcmp(id, CLAP_EXT_STATE))       return &s_dist_state;
    return NULL;
}

static void dist_on_main_thread(const clap_plugin_t *plugin) {}

clap_plugin_t *distort_create(const clap_host_t *host) {
    cdp_dist_t *p = calloc(1, sizeof(*p));
    p->host = host;
    p->clip = 0.5;
    p->depth = 0.5;
    p->freq = 440.0;
    p->exagg = 2.0;
    p->mult = 2.0;
    p->on[0] = 1.0; /* Overload Noise enabled by default */
    p->sample_rate = 44100.0;
    p->plugin.desc = &s_distort_desc;
    p->plugin.plugin_data = p;
    p->plugin.init = dist_init;
    p->plugin.destroy = dist_destroy;
    p->plugin.activate = dist_activate;
    p->plugin.deactivate = dist_deactivate;
    p->plugin.start_processing = dist_start_processing;
    p->plugin.stop_processing = dist_stop_processing;
    p->plugin.reset = dist_reset;
    p->plugin.process = dist_process;
    p->plugin.get_extension = dist_get_extension;
    p->plugin.on_main_thread = dist_on_main_thread;
    return &p->plugin;
}
