# CDP FX Collection

A simple project to bring the [Composers Desktop Project](https://composersdesktop.com) sound processors into your DAW as a **CLAP plugin collection** (VST3 may follow). One binary, many effects — modelled after the [Airwindows Consolidated](https://github.com/baconpaul/airwin2rack) approach.

CDP is ~200 command-line sound-processing tools by Trevor Wishart, in continuous development since the 1980s. They're beloved by sound designers, but running them means command lines and file juggling. This project wraps them as realtime effects with automatable parameters, so you can just drop them on a track.

## Status

Work in progress. Currently ported:

- **CDP Blur** — spectral processing (STFT) with 6 modes from the CDP `blur` tool: **Blur** (time-smear), **Avrg** (bin averaging), **Suppress** (kill loudest partials), **Noise**, **Chorus**, **Scatter**. Each mode has its own on/off toggle — enable several at once and they run in parallel (summed per frame, plus channels split across the host thread pool).
- **CDP Distort** — realtime wavecycle and sample distortion from the CDP `distort` suite with parallel toggles:
  1. **Overload Noise** — peak clipping with pseudo-random grit texture.
  2. **Overload Sine** — peak clipping with sinusoidal frequency modulation.
  3. **Reform Square** — wavecycle-to-square conversion.
  4. **Reform Triangle** — wavecycle-to-triangle conversion.
  5. **Reform Invert** — inverted half-cycle contour.
  6. **Reform Sine** — sine-fitted wavecycle synthesis.
  7. **Reform Exaggerate** — non-linear contour expansion.
  8. **Multiply** — wavecycle frequency multiplier.
- **CDP Gain** — minimal template/plugin plumbing example.

Not every CDP process makes sense as a plugin — some need whole files or text data files. See [docs/effects-catalog.md](docs/effects-catalog.md) for the full inventory and suitability notes.

## Install

1. Build (below) or grab `CDP FX Collection.clap`.
2. Copy it into your CLAP plugins folder:
   - Windows: `C:\Program Files\Common Files\CLAP\`
3. Rescan plugins in your DAW (Reaper, Bitwig, FL Studio, etc.).

## Build

Windows, MSVC:

```
cmake -G "Visual Studio 17 2022" -A x64 -B build
cmake --build build --target cdp_fx --config Release
```

Output lands in `NewRelease/Release/CDP FX Collection.clap`.

To build the full CDP command-line suite instead, just build the whole project (`cmake --build build --config Release`) — executables land in `NewRelease/`.

## Credits & license

- **CDP System Software, Release 8** — Composers Desktop Project / Trevor Wishart, LGPL-2.1+ (see LICENSE). Upstream: [github.com/ComposersDesktop/CDP8](https://github.com/ComposersDesktop/CDP8)
- **CLAP** headers — [free-audio/clap](https://github.com/free-audio/clap), MIT (vendored in `plugin/clap`)
- Plugin wrapper code in `plugin/` — LGPL-2.1+

The original CDP release notes are in [docs/](./docs), and `building.txt` has upstream build details for the CLI tools.
