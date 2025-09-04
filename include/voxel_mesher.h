#pragma once
#include "voxel_world.h"
#include <vector>

class VoxelMesher {
public:
    // Simple greedy meshing could be added later; start with face culling meshing
    static ChunkMesh BuildChunkMesh(const VoxelWorld& world, int cx, int cz);
};
