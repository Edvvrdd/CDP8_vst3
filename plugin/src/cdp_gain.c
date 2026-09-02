/* CDP FX Collection - CLAP plugin skeleton.
 *
 * One binary, many plugins: each CDP port adds a descriptor + create function
 * to s_plugins[] at the bottom of this file (or its own file later).
 * Model: Airwindows Consolidated (github.com/baconpaul/airwin2rack).
 *
 * This first plugin is a minimal gain effect to prove the build and to serve
 * as the copy-me template: audio ports, one automatable param, state save.
 * Debug GUI: a win32 window with "hello world" text and a dummy scrollbar,
 * to verify CLAP_EXT_GUI embedding in a DAW.
 */

#include <clap/clap.h>
#include <clap/ext/gui.h>

#ifdef _WIN32
#include <windows.h>
#endif

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

enum { PARAM_GAIN = 0, NUM_PARAMS };

typedef struct {
   clap_plugin_t     plugin;
   const clap_host_t *host;
   double            sample_rate;
   double            gain_db;
#ifdef _WIN32
   HWND            gui_hwnd;
#endif
} cdp_gain_t;

static const clap_plugin_descriptor_t s_gain_desc = {
   .clap_version = CLAP_VERSION_INIT,
   .id           = "com.composersdesktop.cdp.gain",
   .name         = "CDP Gain",
   .vendor       = "Composers Desktop Project",
   .url          = "https://composersdesktop.com",
   .manual_url   = "",
   .support_url  = "",
   .version      = "0.1.0",
   .description  = "Skeleton plugin for the CDP CLAP collection.",
   .features =
      (const char *[]){CLAP_PLUGIN_FEATURE_AUDIO_EFFECT, CLAP_PLUGIN_FEATURE_STEREO, NULL},
};

/* ---------- audio ports: one stereo in, one stereo out ---------- */

static uint32_t gain_audio_ports_count(const clap_plugin_t *plugin, bool is_input) { return 1; }

static bool gain_audio_ports_get(const clap_plugin_t    *plugin,
                                 uint32_t                index,
                                 bool                    is_input,
                                 clap_audio_port_info_t *info) {
   if (index > 0)
      return false;
   info->id            = 0;
   snprintf(info->name, sizeof(info->name), "Main");
   info->channel_count = 2;
   info->flags         = CLAP_AUDIO_PORT_IS_MAIN;
   info->port_type     = CLAP_PORT_STEREO;
   /* allow host to run us in-place: output reuses input buffers */
   info->in_place_pair = is_input ? CLAP_INVALID_ID : 0;
   return true;
}

static const clap_plugin_audio_ports_t s_gain_audio_ports = {
   .count = gain_audio_ports_count,
   .get   = gain_audio_ports_get,
};

/* ---------- params ---------- */

static uint32_t gain_params_count(const clap_plugin_t *plugin) { return NUM_PARAMS; }

static bool gain_params_get_info(const clap_plugin_t *plugin,
                                 uint32_t             param_index,
                                 clap_param_info_t   *param_info) {
   if (param_index != PARAM_GAIN)
      return false;
   memset(param_info, 0, sizeof(*param_info));
   param_info->id            = PARAM_GAIN;
   snprintf(param_info->name, sizeof(param_info->name), "Gain");
   param_info->min_value     = -60.0;
   param_info->max_value     = 24.0;
   param_info->default_value = 0.0;
   param_info->flags         = CLAP_PARAM_IS_AUTOMATABLE;
   return true;
}

static bool gain_params_get_value(const clap_plugin_t *plugin, clap_id param_id, double *value) {
   const cdp_gain_t *p = plugin->plugin_data;
   if (param_id != PARAM_GAIN)
      return false;
   *value = p->gain_db;
   return true;
}

static bool gain_params_value_to_text(const clap_plugin_t *plugin,
                                      clap_id              param_id,
                                      double               value,
                                      char                *display,
                                      uint32_t             size) {
   if (param_id != PARAM_GAIN)
      return false;
   snprintf(display, size, "%.1f dB", value);
   return true;
}

static bool gain_params_text_to_value(const clap_plugin_t *plugin,
                                      clap_id              param_id,
                                      const char          *display,
                                      double              *value) {
   char *end = NULL;
   double v;
   if (param_id != PARAM_GAIN)
      return false;
   v = strtod(display, &end);
   if (end == display)
      return false;
   if (v < -60.0)
      v = -60.0;
   if (v > 24.0)
      v = 24.0;
   *value = v;
   return true;
}

static void gain_apply_event(cdp_gain_t *p, const clap_event_header_t *hdr) {
   if (hdr->space_id == CLAP_CORE_EVENT_SPACE_ID && hdr->type == CLAP_EVENT_PARAM_VALUE) {
      const clap_event_param_value_t *ev = (const clap_event_param_value_t *)hdr;
      if (ev->param_id == PARAM_GAIN)
         p->gain_db = ev->value;
   }
}

static void gain_params_flush(const clap_plugin_t        *plugin,
                              const clap_input_events_t  *in,
                              const clap_output_events_t *out) {
   cdp_gain_t *p = plugin->plugin_data;
   const uint32_t nev = in->size(in);
   for (uint32_t i = 0; i < nev; ++i)
      gain_apply_event(p, in->get(in, i));
}

static const clap_plugin_params_t s_gain_params = {
   .count         = gain_params_count,
   .get_info      = gain_params_get_info,
   .get_value     = gain_params_get_value,
   .value_to_text = gain_params_value_to_text,
   .text_to_value = gain_params_text_to_value,
   .flush         = gain_params_flush,
};

/* ---------- state: one double ----------
 * ponytail: raw double, not portable across endianness; use a text/agnostic
 * format if presets must move between ARM/x86 machines. */

static bool gain_state_save(const clap_plugin_t *plugin, const clap_ostream_t *stream) {
   const cdp_gain_t *p = plugin->plugin_data;
   return stream->write(stream, &p->gain_db, sizeof(p->gain_db)) == sizeof(p->gain_db);
}

static bool gain_state_load(const clap_plugin_t *plugin, const clap_istream_t *stream) {
   cdp_gain_t *p = plugin->plugin_data;
   double v;
   if (stream->read(stream, &v, sizeof(v)) != sizeof(v))
      return false;
   if (v < -60.0)
      v = -60.0;
   if (v > 24.0)
      v = 24.0;
   p->gain_db = v;
   return true;
}

static const clap_plugin_state_t s_gain_state = {
   .save = gain_state_save,
   .load = gain_state_load,
};

/* ---------- plugin ---------- */

static bool gain_init(const clap_plugin_t *plugin) { return true; }

static void gain_destroy(const clap_plugin_t *plugin) { free(plugin->plugin_data); }

static bool gain_activate(const clap_plugin_t *plugin,
                          double               sample_rate,
                          uint32_t             min_frames_count,
                          uint32_t             max_frames_count) {
   cdp_gain_t *p = plugin->plugin_data;
   p->sample_rate = sample_rate;
   return true;
}

static void gain_deactivate(const clap_plugin_t *plugin) {}

static bool gain_start_processing(const clap_plugin_t *plugin) { return true; }

static void gain_stop_processing(const clap_plugin_t *plugin) {}

static void gain_reset(const clap_plugin_t *plugin) {}

static float read_sample(const clap_audio_buffer_t *buf, uint32_t ch, uint32_t frame) {
   return buf->data32 ? buf->data32[ch][frame] : (float)buf->data64[ch][frame];
}

static void write_sample(clap_audio_buffer_t *buf, uint32_t ch, uint32_t frame, float v) {
   if (buf->data32)
      buf->data32[ch][frame] = v;
   else
      buf->data64[ch][frame] = v;
}

/* ponytail: params applied at block start, not sample-accurate; add per-event
 * sub-slicing and gain smoothing if zipper noise becomes audible. */
static clap_process_status gain_process(const clap_plugin_t *plugin, const clap_process_t *process) {
   cdp_gain_t *p = plugin->plugin_data;
   const clap_input_events_t *evs = process->in_events;
   const uint32_t nev = evs ? evs->size(evs) : 0;
   const clap_audio_buffer_t *in;
   clap_audio_buffer_t *out;
   uint32_t nch;
   double g;

   for (uint32_t i = 0; i < nev; ++i)
      gain_apply_event(p, evs->get(evs, i));

   if (process->audio_outputs_count == 0)
      return CLAP_PROCESS_CONTINUE;

   in  = process->audio_inputs_count > 0 ? &process->audio_inputs[0] : NULL;
   out = &process->audio_outputs[0];
   nch = out->channel_count;
   if (in && in->channel_count < nch)
      nch = in->channel_count;
   g = pow(10.0, p->gain_db / 20.0);

   for (uint32_t c = 0; c < nch; ++c)
      for (uint32_t i = 0; i < process->frames_count; ++i) {
         const float s = in ? read_sample(in, c, i) : 0.0f;
         write_sample(out, c, i, (float)(s * g));
      }

   return CLAP_PROCESS_CONTINUE;
}

/* ---------- debug GUI: win32 window, hello world + dummy scrollbar ---------- */

#ifdef _WIN32

#define GUI_W 320
#define GUI_H 120

static bool gui_is_api_supported(const clap_plugin_t *plugin, const char *api, bool is_floating) {
   return !strcmp(api, CLAP_WINDOW_API_WIN32);
}

static bool gui_get_preferred_api(const clap_plugin_t *plugin, const char **api, bool *is_floating) {
   *api        = CLAP_WINDOW_API_WIN32;
   *is_floating = false;
   return true;
}

static bool gui_create(const clap_plugin_t *plugin, const char *api, bool is_floating) {
   cdp_gain_t *p     = plugin->plugin_data;
   HINSTANCE   inst  = GetModuleHandle(NULL);
   WNDCLASSW   wc    = {0};
   const wchar_t *cls = L"CDPDebugGui";
   DWORD style;

   if (api && *api && strcmp(api, CLAP_WINDOW_API_WIN32))
      return false;

   if (!GetClassInfoW(inst, cls, &wc)) {
      wc.lpfnWndProc   = DefWindowProcW;
      wc.hInstance     = inst;
      wc.lpszClassName = cls;
      wc.hbrBackground = (HBRUSH)(COLOR_WINDOW + 1);
      if (!RegisterClassW(&wc))
         return false;
   }

   style = is_floating ? (WS_POPUPWINDOW | WS_CAPTION | WS_SYSMENU) : (WS_CHILD | WS_VISIBLE);
   /* set_parent() comes later; create detached */
   p->gui_hwnd = CreateWindowW(cls, L"CDP FX debug", style, 0, 0, GUI_W, GUI_H,
                               NULL, NULL, inst, NULL);
   if (!p->gui_hwnd)
      return false;

   CreateWindowW(L"STATIC", L"hello world", WS_CHILD | WS_VISIBLE | SS_CENTER,
                 0, 0, GUI_W, 40, p->gui_hwnd, NULL, inst, NULL);

   HWND sb = CreateWindowW(L"SCROLLBAR", NULL, WS_CHILD | WS_VISIBLE | SBS_HORZ,
                           0, 50, GUI_W, 20, p->gui_hwnd, NULL, inst, NULL);
   SCROLLINFO si = {sizeof(si), SIF_ALL, 0, 100, 1, 50, 0};
   SetScrollInfo(sb, SB_CTL, &si, TRUE);

   return true;
}

static void gui_destroy(const clap_plugin_t *plugin) {
   cdp_gain_t *p = plugin->plugin_data;
   if (p->gui_hwnd)
      DestroyWindow(p->gui_hwnd);
   p->gui_hwnd = NULL;
}

static bool gui_set_scale(const clap_plugin_t *plugin, double scale) { return true; }

static bool gui_get_size(const clap_plugin_t *plugin, uint32_t *width, uint32_t *height) {
   *width  = GUI_W;
   *height = GUI_H;
   return true;
}

static bool gui_can_resize(const clap_plugin_t *plugin) { return false; }

static bool gui_get_resize_hints(const clap_plugin_t *plugin, clap_gui_resize_hints_t *hints) {
   return false;
}

static bool gui_adjust_size(const clap_plugin_t *plugin, uint32_t *width, uint32_t *height) {
   *width  = GUI_W;
   *height = GUI_H;
   return true;
}

static bool gui_set_size(const clap_plugin_t *plugin, uint32_t width, uint32_t height) {
   return width == GUI_W && height == GUI_H;
}

static bool gui_set_parent(const clap_plugin_t *plugin, const clap_window_t *window) {
   cdp_gain_t *p = plugin->plugin_data;
   DWORD style;

   if (strcmp(window->api, CLAP_WINDOW_API_WIN32) || !window->win32 || !p->gui_hwnd)
      return false;

   style = GetWindowLong(p->gui_hwnd, GWL_STYLE) | WS_CHILD;
   SetWindowLong(p->gui_hwnd, GWL_STYLE, style);
   SetParent(p->gui_hwnd, window->win32);
   return true;
}

static bool gui_set_transient(const clap_plugin_t *plugin, const clap_window_t *window) {
   return false;
}

static void gui_suggest_title(const clap_plugin_t *plugin, const char *title) {}

static bool gui_show(const clap_plugin_t *plugin) {
   cdp_gain_t *p = plugin->plugin_data;
   if (!p->gui_hwnd)
      return false;
   ShowWindow(p->gui_hwnd, SW_SHOW);
   return true;
}

static bool gui_hide(const clap_plugin_t *plugin) {
   cdp_gain_t *p = plugin->plugin_data;
   if (!p->gui_hwnd)
      return false;
   ShowWindow(p->gui_hwnd, SW_HIDE);
   return true;
}

static const clap_plugin_gui_t s_gain_gui = {
   .is_api_supported  = gui_is_api_supported,
   .get_preferred_api = gui_get_preferred_api,
   .create            = gui_create,
   .destroy           = gui_destroy,
   .set_scale         = gui_set_scale,
   .get_size          = gui_get_size,
   .can_resize        = gui_can_resize,
   .get_resize_hints  = gui_get_resize_hints,
   .adjust_size       = gui_adjust_size,
   .set_size          = gui_set_size,
   .set_parent        = gui_set_parent,
   .set_transient     = gui_set_transient,
   .suggest_title     = gui_suggest_title,
   .show              = gui_show,
   .hide              = gui_hide,
};

#endif /* _WIN32 */

static const void *gain_get_extension(const clap_plugin_t *plugin, const char *id) {
   if (!strcmp(id, CLAP_EXT_AUDIO_PORTS))
      return &s_gain_audio_ports;
   if (!strcmp(id, CLAP_EXT_PARAMS))
      return &s_gain_params;
   if (!strcmp(id, CLAP_EXT_STATE))
      return &s_gain_state;
#ifdef _WIN32
   if (!strcmp(id, CLAP_EXT_GUI))
      return &s_gain_gui;
#endif
   return NULL;
}

static void gain_on_main_thread(const clap_plugin_t *plugin) {}

static clap_plugin_t *gain_create(const clap_host_t *host) {
   cdp_gain_t *p = calloc(1, sizeof(*p));
   p->host                  = host;
   p->gain_db               = 0.0;
   p->plugin.desc           = &s_gain_desc;
   p->plugin.plugin_data    = p;
   p->plugin.init           = gain_init;
   p->plugin.destroy        = gain_destroy;
   p->plugin.activate       = gain_activate;
   p->plugin.deactivate     = gain_deactivate;
   p->plugin.start_processing = gain_start_processing;
   p->plugin.stop_processing  = gain_stop_processing;
   p->plugin.reset          = gain_reset;
   p->plugin.process        = gain_process;
   p->plugin.get_extension  = gain_get_extension;
   p->plugin.on_main_thread = gain_on_main_thread;
   return &p->plugin;
}

/* ---------- factory: add one entry per ported CDP process ---------- */

extern const clap_plugin_descriptor_t s_blur_desc;
extern clap_plugin_t *blur_create(const clap_host_t *host);

static struct {
   const clap_plugin_descriptor_t *desc;
   clap_plugin_t *(CLAP_ABI *create)(const clap_host_t *host);
} s_plugins[] = {
   {.desc = &s_gain_desc, .create = gain_create},
   {.desc = &s_blur_desc, .create = blur_create},
};

static uint32_t plugin_factory_get_plugin_count(const clap_plugin_factory_t *factory) {
   return sizeof(s_plugins) / sizeof(s_plugins[0]);
}

static const clap_plugin_descriptor_t *
plugin_factory_get_plugin_descriptor(const clap_plugin_factory_t *factory, uint32_t index) {
   return s_plugins[index].desc;
}

static const clap_plugin_t *plugin_factory_create_plugin(const clap_plugin_factory_t *factory,
                                                         const clap_host_t          *host,
                                                         const char                 *plugin_id) {
   const int n = sizeof(s_plugins) / sizeof(s_plugins[0]);
   if (!clap_version_is_compatible(host->clap_version))
      return NULL;
   for (int i = 0; i < n; ++i)
      if (!strcmp(plugin_id, s_plugins[i].desc->id))
         return s_plugins[i].create(host);
   return NULL;
}

static const clap_plugin_factory_t s_plugin_factory = {
   .get_plugin_count      = plugin_factory_get_plugin_count,
   .get_plugin_descriptor = plugin_factory_get_plugin_descriptor,
   .create_plugin         = plugin_factory_create_plugin,
};

/* ---------- entry ---------- */

static int s_entry_init_counter = 0;

static bool entry_init(const char *plugin_path) {
   ++s_entry_init_counter;
   return true;
}

static void entry_deinit(void) { --s_entry_init_counter; }

static const void *entry_get_factory(const char *factory_id) {
   if (!strcmp(factory_id, CLAP_PLUGIN_FACTORY_ID))
      return &s_plugin_factory;
   return NULL;
}

CLAP_EXPORT const clap_plugin_entry_t clap_entry = {
   .clap_version = CLAP_VERSION_INIT,
   .init         = entry_init,
   .deinit       = entry_deinit,
   .get_factory  = entry_get_factory,
};
