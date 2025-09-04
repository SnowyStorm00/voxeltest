#include "voxel_mesher.h"
#include <array>
using namespace DirectX;

static const std::array<XMFLOAT3, 6> normals = {
    XMFLOAT3{-1,0,0}, XMFLOAT3{1,0,0}, XMFLOAT3{0,0,-1}, XMFLOAT3{0,0,1}, XMFLOAT3{0,1,0}, XMFLOAT3{0,-1,0}
};

static const std::array<std::array<XMFLOAT3,4>,6> faceVerts = {
    std::array<XMFLOAT3,4>{ XMFLOAT3{0,0,0}, XMFLOAT3{0,1,0}, XMFLOAT3{0,1,1}, XMFLOAT3{0,0,1} }, // -X
    std::array<XMFLOAT3,4>{ XMFLOAT3{1,0,1}, XMFLOAT3{1,1,1}, XMFLOAT3{1,1,0}, XMFLOAT3{1,0,0} }, // +X
    std::array<XMFLOAT3,4>{ XMFLOAT3{1,0,0}, XMFLOAT3{1,1,0}, XMFLOAT3{0,1,0}, XMFLOAT3{0,0,0} }, // -Z
    std::array<XMFLOAT3,4>{ XMFLOAT3{0,0,1}, XMFLOAT3{0,1,1}, XMFLOAT3{1,1,1}, XMFLOAT3{1,0,1} }, // +Z
    std::array<XMFLOAT3,4>{ XMFLOAT3{0,1,1}, XMFLOAT3{0,1,0}, XMFLOAT3{1,1,0}, XMFLOAT3{1,1,1} }, // +Y
    std::array<XMFLOAT3,4>{ XMFLOAT3{0,0,0}, XMFLOAT3{0,0,1}, XMFLOAT3{1,0,1}, XMFLOAT3{1,0,0} }, // -Y
};

static inline uint32_t packColor(uint8_t r, uint8_t g, uint8_t b) {
    return (0xFFu << 24) | (uint32_t(r) << 16) | (uint32_t(g) << 8) | uint32_t(b);
}

ChunkMesh VoxelMesher::BuildChunkMesh(const VoxelWorld& world, int cx, int cz) {
    ChunkMesh mesh;
    const int baseX = cx * CHUNK_SIZE;
    const int baseZ = cz * CHUNK_SIZE;

    for(int z = 0; z < CHUNK_SIZE; ++z) {
        for(int x = 0; x < CHUNK_SIZE; ++x) {
            for(int y = 0; y < CHUNK_HEIGHT; ++y) {
                uint8_t id = world.GetBlock(baseX + x, y, baseZ + z);
                if(id == 0) continue;

                // Pick a color from id
                uint32_t color = packColor(200,200,200);
                if(id == 2) color = packColor(80,200,120);      // grass
                else if(id == 3) color = packColor(155,118,83); // dirt
                else if(id == 4) color = packColor(150,150,160); // stone

                // For each face, if neighbor is empty, emit face
                const int nx[6] = {-1,1,0,0,0,0};
                const int ny[6] = {0,0,0,0,1,-1};
                const int nz[6] = {0,0,-1,1,0,0};
                for(int f=0; f<6; ++f) {
                    uint8_t nid = world.GetBlock(baseX + x + nx[f], y + ny[f], baseZ + z + nz[f]);
                    if(nid != 0) continue;
                    uint32_t start = (uint32_t)mesh.vertices.size();
                    for(int v=0; v<4; ++v) {
                        XMFLOAT3 p = faceVerts[f][v];
                        mesh.vertices.push_back({ XMFLOAT3{p.x + x + baseX, p.y + y, p.z + z + baseZ}, normals[f], color });
                    }
                    // two triangles
                    // Ensure counter-clockwise winding for front faces
                    mesh.indices.push_back(start + 0);
                    mesh.indices.push_back(start + 1);
                    mesh.indices.push_back(start + 2);
                    mesh.indices.push_back(start + 0);
                    mesh.indices.push_back(start + 2);
                    mesh.indices.push_back(start + 3);
                }
            }
        }
    }

    return mesh;
}
