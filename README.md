# VoxelRPG (D3D11, C++)

Windows voxel sandbox built with Direct3D 11 and CMake. Ships as a self-contained ZIP via CPack.

## Highlights (v0.2.0)
- Chunked voxel world (16×16×128) with greedy meshing.
- CPU-baked sun shadows at meshing time (3D DDA): per-vertex shadow stored in alpha, adaptive quad subdivision where variance is high, water excluded from darkening.
- Quality controls: presets from FASTEST → ACCURATE and live sliders (samples, distance, radius, subdivision/variance, hard/soft shadows, sun samples).
- Sky pass with gradient and sun disk.
- Smooth frame pacing: precise FPS limiter (QPC) and optional VSync; FPS overlay.
- Loading UX: startup seed input screen; blocking loading screen with progress bar; fast initial meshing path followed by HQ remesh.
- Freecam: press ` (backtick) to toggle; WASD + mouse fly; Shift boost; Space/Ctrl up/down; player model hidden in freecam.
- Terrain details: tan paths, shallow ponds (opaque water voxel).

## Build

Prerequisites
- Visual Studio 2019/2022 (Desktop development with C++)
- Windows 10/11 SDK
- CMake

Configure and build (Developer PowerShell)
```
cmake -S . -B build -DCMAKE_BUILD_TYPE=RelWithDebInfo
cmake --build build --config RelWithDebInfo
```

Run (from repo root)
```
& "build/RelWithDebInfo/voxel.exe"
```

Packaging (ZIP)
```
cmake --build build --config RelWithDebInfo --target PACKAGE
# Output: build/VoxelRPG-0.2.0-win64.zip
```

## Controls
- Mouse: look
- WASD: move
- Shift: sprint / boost (freecam)
- Space / Ctrl: ascend / descend (freecam)
- ` (backtick): toggle freecam
- Enter (on seed screen): confirm seed and start

## Settings & UI
- Overlay shows FPS and a settings menu:
	- Presets: FASTEST → ACCURATE
	- Sliders: shadow samples, max distance, penumbra radius, adaptive subdivision threshold, hard/soft classification, sun samples
	- VSync toggle and FPS cap
- Changes trigger live remesh of visible chunks.

## Notes / Limitations
- Water is opaque and unlit by baked shadow (intended). No transparency/refraction yet.
- Freecam has no collision or physics.
- Hardware shadow maps were removed in favor of stable CPU-baked shadows.

## License
MIT (see LICENSE).
