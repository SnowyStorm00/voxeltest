# Voxel RPG Starter (D3D11, C++)

Minimal voxel world rendering on Windows using Direct3D 11 and CMake. No external dependencies.

Features:
- Window + D3D11 device setup
- Simple camera and matrices (DirectXMath)
- Voxel world with chunked storage (16x64x16)
- Face-culling mesher
- Basic directional light shader

## Build

Prereqs: Visual Studio 2019/2022 with C++ and CMake, Windows 10 SDK.

### Configure and build
Use CMake Presets or from a Developer PowerShell:

```
cmake -S . -B build -DCMAKE_BUILD_TYPE=RelWithDebInfo
cmake --build build --config RelWithDebInfo
```

Run the executable from the build dir so shaders can be found:

```
./build/RelWithDebInfo/voxel.exe
```

## Next steps
- Add input handling (WASD + mouse)
- Greedy meshing for fewer triangles
- Chunk streaming and saving
- Block types and textures (atlas)
- Basic gameplay loop: player, NPCs, items

## Settings Menu:
- VSYNC: toggle vertical sync.
- FPS CAP: cycle caps; Enter to set uncapped.
- TAA: toggle Temporal AA with Halton jitter. When ON, the scene renders to an offscreen target and resolves with a simple history clamp.
