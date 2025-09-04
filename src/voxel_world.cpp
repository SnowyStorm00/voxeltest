#include "voxel_world.h"
#include <algorithm>
#include <cmath>
#include <random>
#include <cmath>

static inline int idx(int x, int y, int z) {
    return x + CHUNK_SIZE*(z + CHUNK_SIZE*y);
}

VoxelWorld::VoxelWorld(int radius) : m_radius(radius) {
    int dim = 2*radius + 1;
    m_chunks.resize(dim * dim);
}

bool VoxelWorld::InBounds(int cx, int cz) const {
    int r = m_radius;
    return cx >= -r && cx <= r && cz >= -r && cz <= r;
}

const Chunk& VoxelWorld::GetChunk(int cx, int cz) const {
    int r = m_radius;
    int dim = 2*r + 1;
    return m_chunks[(cx + r) + dim * (cz + r)];
}

Chunk& VoxelWorld::GetChunk(int cx, int cz) {
    int r = m_radius;
    int dim = 2*r + 1;
    return m_chunks[(cx + r) + dim * (cz + r)];
}

uint8_t VoxelWorld::GetBlock(int x, int y, int z) const {
    int cx = floorf((float)x / CHUNK_SIZE);
    int cz = floorf((float)z / CHUNK_SIZE);
    if(!InBounds(cx, cz) || y < 0 || y >= CHUNK_HEIGHT) return 0;
    int lx = (x % CHUNK_SIZE + CHUNK_SIZE) % CHUNK_SIZE;
    int lz = (z % CHUNK_SIZE + CHUNK_SIZE) % CHUNK_SIZE;
    auto& c = m_chunks[(cx + m_radius) + (2*m_radius + 1) * (cz + m_radius)];
    return c.blocks[idx(lx, y, lz)];
}

void VoxelWorld::SetBlock(int x, int y, int z, uint8_t id) {
    int cx = floorf((float)x / CHUNK_SIZE);
    int cz = floorf((float)z / CHUNK_SIZE);
    if(!InBounds(cx, cz) || y < 0 || y >= CHUNK_HEIGHT) return;
    int lx = (x % CHUNK_SIZE + CHUNK_SIZE) % CHUNK_SIZE;
    int lz = (z % CHUNK_SIZE + CHUNK_SIZE) % CHUNK_SIZE;
    auto& c = m_chunks[(cx + m_radius) + (2*m_radius + 1) * (cz + m_radius)];
    c.blocks[idx(lx, y, lz)] = id;
}

void VoxelWorld::GenerateFlat() {
    int r = m_radius;
    for(int cz = -r; cz <= r; ++cz) {
        for(int cx = -r; cx <= r; ++cx) {
            auto& c = GetChunk(cx, cz);
            for(int z = 0; z < CHUNK_SIZE; ++z) {
                for(int x = 0; x < CHUNK_SIZE; ++x) {
                    for(int y = 0; y < CHUNK_HEIGHT; ++y) {
                        uint8_t id = y < 8 ? 2 : 0; // flat ground
                        c.blocks[idx(x,y,z)] = id;
                    }
                }
            }
        }
    }
}

// Simple 2D value noise with bilinear interpolation and multiple octaves
static float hash2i(int x, int y, uint32_t seed) {
    uint32_t h = (uint32_t)x * 374761393u + (uint32_t)y * 668265263u + seed * 2246822519u;
    h = (h ^ (h >> 13)) * 1274126177u;
    return (h ^ (h >> 16)) * (1.0f/4294967295.0f);
}

static float valueNoise2D(float x, float y, uint32_t seed) {
    int xi = (int)floorf(x);
    int yi = (int)floorf(y);
    float tx = x - xi;
    float ty = y - yi;
    auto lerp = [](float a, float b, float t){ return a + (b - a) * (t * t * (3 - 2*t)); };
    float v00 = hash2i(xi, yi, seed);
    float v10 = hash2i(xi+1, yi, seed);
    float v01 = hash2i(xi, yi+1, seed);
    float v11 = hash2i(xi+1, yi+1, seed);
    float vx0 = lerp(v00, v10, tx);
    float vx1 = lerp(v01, v11, tx);
    return lerp(vx0, vx1, ty);
}

static float fbm2D(float x, float y, uint32_t seed, int octaves, float lacunarity, float gain) {
    float amp = 1.0f;
    float freq = 1.0f;
    float sum = 0.0f;
    float norm = 0.0f;
    for(int i=0;i<octaves;++i){
        sum += valueNoise2D(x*freq, y*freq, seed + i*1013u) * amp;
        norm += amp;
        amp *= gain;
        freq *= lacunarity;
    }
    return sum / (norm > 0 ? norm : 1.0f);
}

void VoxelWorld::GenerateProcedural(unsigned seed) {
    int r = m_radius;
    float baseScale = 1.0f / 48.0f; // big features
    float detailScale = 1.0f / 12.0f; // smaller shapes
    for(int cz = -r; cz <= r; ++cz) {
        for(int cx = -r; cx <= r; ++cx) {
            auto& c = GetChunk(cx, cz);
            for(int z = 0; z < CHUNK_SIZE; ++z) {
                for(int x = 0; x < CHUNK_SIZE; ++x) {
                    int wx = cx*CHUNK_SIZE + x;
                    int wz = cz*CHUNK_SIZE + z;
                    // Height via FBM
                    float h0 = fbm2D(wx * baseScale, wz * baseScale, seed, 3, 2.0f, 0.5f);
                    float h1 = fbm2D(wx * detailScale, wz * detailScale, seed+1337u, 2, 2.0f, 0.5f);
                    float h = 0.4f*h0 + 0.6f*h1; // blend
                    int height = (int)floorf(8 + h * 18); // 0..~26
                    height = std::clamp(height, 0, CHUNK_HEIGHT-1);
                    for(int y = 0; y < CHUNK_HEIGHT; ++y) {
                        uint8_t id = 0;
                        if(y <= height) {
                            if(y == height) id = 2; // grass
                            else if(y >= height-3) id = 3; // dirt
                            else id = 4; // stone
                        }
                        c.blocks[idx(x,y,z)] = id;
                    }
                }
            }
        }
    }
}
