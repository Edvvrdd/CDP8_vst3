/* scan_check: load the built .clap, enumerate the factory, and smoke-test
 * every plugin: instantiate, activate, and process audio with BOTH data32
 * and data64 buffers (hosts may deliver either; a NULL deref here = host
 * crash/blacklist in a real DAW).
 *
 * usage: scan_check "path/to/CDP FX Collection.clap"
 */

#include <clap/clap.h>
#include <windows.h>
#include <stdio.h>
#include <string.h>

#define FRAMES 128
#define NCH 2

static const void *host_get_extension(const clap_host_t *host, const char *id) { return NULL; }
static void host_request_restart(const clap_host_t *host) {}
static void host_request_process(const clap_host_t *host) {}
static void host_request_callback(const clap_host_t *host) {}

int main(int argc, char **argv) {
    if (argc < 2) { printf("usage: scan_check <path-to.clap>\n"); return 2; }

    HMODULE dll = LoadLibraryA(argv[1]);
    if (!dll) { printf("FAIL: LoadLibrary error %lu\n", GetLastError()); return 1; }

    const clap_plugin_entry_t *entry = (const clap_plugin_entry_t *)GetProcAddress(dll, "clap_entry");
    if (!entry) { printf("FAIL: no clap_entry\n"); return 1; }
    if (!entry->init(argv[1])) { printf("FAIL: entry init\n"); return 1; }

    const clap_plugin_factory_t *f =
        (const clap_plugin_factory_t *)entry->get_factory(CLAP_PLUGIN_FACTORY_ID);
    if (!f) { printf("FAIL: no plugin factory\n"); return 1; }

    clap_host_t host = {
        .clap_version = CLAP_VERSION_INIT,
        .host_data = NULL,
        .name = "scan-check",
        .vendor = "CDP",
        .url = "https://localhost",
        .version = "0.0.1",
        .get_extension = host_get_extension,
        .request_restart = host_request_restart,
        .request_process = host_request_process,
        .request_callback = host_request_callback,
    };

    float  in32[NCH][FRAMES] = {{0}}, out32[NCH][FRAMES] = {{0}};
    double in64[NCH][FRAMES] = {{0}}, out64[NCH][FRAMES] = {{0}};
    for (int c = 0; c < NCH; ++c)
        for (int i = 0; i < FRAMES; ++i) {
            in32[c][i] = 0.5f * (i & 1 ? 1.0f : -1.0f) * ((i % 16) / 16.0f);
            in64[c][i] = in32[c][i];
        }

    uint32_t n = f->get_plugin_count(f);
    printf("%u plugins:\n", n);

    for (uint32_t i = 0; i < n; ++i) {
        const clap_plugin_descriptor_t *d = f->get_plugin_descriptor(f, i);
        printf("  [%u] %-40s %s\n", i, d->id, d->name);

        const clap_plugin_t *plg = f->create_plugin(f, &host, d->id);
        if (!plg) { printf("      FAIL: create\n"); continue; }
        if (!plg->init(plg)) { printf("      FAIL: init\n"); plg->destroy(plg); continue; }
        if (!plg->activate(plg, 44100.0, 32, 512)) { printf("      FAIL: activate\n"); plg->destroy(plg); continue; }

        for (int pass = 0; pass < 2; ++pass) {   /* pass 0: data32, pass 1: data64 */
            clap_audio_buffer_t inb = {0}, outb = {0};
            float  *p32in[NCH]  = {in32[0],  in32[1]};
            float  *p32out[NCH] = {out32[0], out32[1]};
            double *p64in[NCH]  = {in64[0],  in64[1]};
            double *p64out[NCH] = {out64[0], out64[1]};
            if (pass == 0) { inb.data32 = p32in; outb.data32 = p32out; }
            else           { inb.data64 = p64in; outb.data64 = p64out; }
            inb.channel_count = outb.channel_count = NCH;
            inb.latency = outb.latency = 0;

            clap_process_t proc = {0};
            proc.steady_time = -1;
            proc.frames_count = FRAMES;
            proc.transport = NULL;
            proc.audio_inputs = &inb;
            proc.audio_outputs = &outb;
            proc.audio_inputs_count = 1;
            proc.audio_outputs_count = 1;
            proc.in_events = NULL;
            proc.out_events = NULL;

            clap_process_status st = plg->process(plg, &proc);
            if (st == CLAP_PROCESS_ERROR) {
                printf("      FAIL: process (%s)\n", pass ? "data64" : "data32");
                break;
            }
            /* NaN/inf check on output */
            int bad = 0;
            for (int c = 0; c < NCH && !bad; ++c)
                for (int k = 0; k < FRAMES; ++k) {
                    double v = pass ? out64[c][k] : out32[c][k];
                    if (v != v || v > 1e30 || v < -1e30) { bad = 1; break; }
                }
            if (bad) { printf("      FAIL: NaN/inf output (%s)\n", pass ? "data64" : "data32"); break; }
        }

        plg->deactivate(plg);
        plg->destroy(plg);
        printf("      OK\n");
    }

    entry->deinit();
    FreeLibrary(dll);
    return 0;
}
