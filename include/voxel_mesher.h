#pragma once
#include "voxel_world.h"
#include <vector>

// Shadow / meshing quality knobs used by the software shadow baker.
struct MesherSettings {
    int   maxSteps = 96;          // DDA steps (shadow distance)
    int   vertexSamples = 4;      // per-vertex in-plane samples (1,2,4,8)
    float sampleRadius = 0.15f;   // in-plane sampling radius in voxel units
    float normalBias   = 0.02f;   // push-off along face normal to avoid self-hit
    float sunBias      = 0.02f;   // tiny push toward sun to avoid grazing hits
    float inset        = 0.08f;   // inset from vertex toward face center
    float minShadow    = 0.05f;   // clamp minimum light (avoid full black)
    float mergeLitThreshold = 0.60f; // center-sample lit threshold for greedy merge split
    // Subdivision to reduce interpolation artifacts across big quads
    bool  enableSubdivision = true; // split quads when vertex shadows differ a lot
    float subdivThreshold   = 0.20f; // if (max-min) > threshold, subdivide
    int   subdivMaxDepth    = 6;     // recursive depth cap (6 => up to ~64x64 split if needed)
    // Rendering style: hard vs. smooth baked shadows
    bool  hardShadows       = true;  // when true, binarize per-quad shadow with mergeLitThreshold
    // Soft shadow (area sun) approximation
    int   sunSamples        = 1;     // >1 enables multi-direction sampling for soft penumbra
    float sunAngularRadius  = 0.0092f; // radians (~0.53 degrees)
    // Initial-load fast path: skip all shadow sampling and subdivision
    bool  fastBakeOnly      = false;  // when true, emit quads with no shadow tests (alpha=1)
};

class VoxelMesher {
public:
    // Builds a chunk mesh and bakes per-vertex shadowing using CPU DDA.
    static ChunkMesh BuildChunkMesh(const VoxelWorld& world, int cx, int cz, const DirectX::XMFLOAT3& sunDir, const MesherSettings& settings);
};
