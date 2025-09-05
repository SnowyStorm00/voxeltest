#include "voxel_world.h"
#include <algorithm>
#include <cmath>
#include <random>

static inline int idx(int x, int y, int z) {
    return x + CHUNK_SIZE*(z + CHUNK_SIZE*y);
}

VoxelWorld::VoxelWorld(unsigned seed) : m_seed(seed) {}

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

void VoxelWorld::GenerateChunkData(Chunk& c, int cx, int cz) const {
    // Discrete biomes: 0=plains, 1=hills, 2=mountains, 3=desert
    const float biomeScale = 1.0f / 256.0f; // very low frequency for large regions
    for(int z = 0; z < CHUNK_SIZE; ++z) {
        for(int x = 0; x < CHUNK_SIZE; ++x) {
            int wx = cx*CHUNK_SIZE + x;
            int wz = cz*CHUNK_SIZE + z;

            // Biome control: bias toward lower values to increase plains frequency
            float b = fbm2D(wx * biomeScale, wz * biomeScale, m_seed + 9001u, 3, 2.0f, 0.5f);
            b = 0.5f * (b + 1.0f);
            b = powf(std::clamp(b, 0.0f, 1.0f), 1.2f);
            // Desert mask is computed later as a smooth factor for blending (not a hard boolean)
            // Smoothly blend heights in narrow bands around thresholds to avoid cliffs
            const float t0 = 0.65f, t1 = 0.90f; // thresholds
            const float w = 0.06f;              // blend half-width
            auto sstep = [](float e0, float e1, float x){ float t = (x - e0) / (e1 - e0); t = std::clamp(t, 0.0f, 1.0f); return t*t*(3.0f - 2.0f*t); };

            // Height sampler for a specific biome index (0=plains,1=hills,2=mountains,3=desert)
            auto sampleHeight = [&](int biomeIdx)->std::pair<float,uint8_t>{
                float baseHeight, amp, freq; int oct; float warpScale, warpAmp; float ridgedMix = 0.0f; uint8_t grassId;
                if(biomeIdx == 0){ // plains
                    baseHeight = 10.0f; amp = 16.0f; freq = 1.0f/72.0f; oct=3; warpScale=1.0f/200.0f; warpAmp=3.0f; grassId=2; ridgedMix=0.0f;
                } else if(biomeIdx == 1){ // hills
                    baseHeight = 14.0f; amp = 36.0f; freq = 1.0f/44.0f; oct=5; warpScale=1.0f/150.0f; warpAmp=10.0f; grassId=5; ridgedMix=0.35f;
                } else if(biomeIdx == 2){ // mountains
                    baseHeight = 18.0f; amp = 72.0f; freq = 1.0f/28.0f; oct=5; warpScale=1.0f/120.0f; warpAmp=22.0f; grassId=6; ridgedMix=0.85f;
                } else { // desert
                    // Flat like plains, with occasional larger dunes/hills
                    baseHeight = 9.0f; amp = 10.0f; freq = 1.0f/90.0f; oct=3; warpScale=1.0f/240.0f; warpAmp=2.0f; grassId=9; // 9 = sand
                }

                float wxWarp = fbm2D(wx * warpScale, wz * warpScale, m_seed + 1717u, 2, 2.0f, 0.5f) * warpAmp;
                float wzWarp = fbm2D(wx * warpScale, wz * warpScale, m_seed + 2727u, 2, 2.0f, 0.5f) * warpAmp;
                float n = fbm2D((wx + wxWarp) * freq, (wz + wzWarp) * freq, m_seed + 17u, oct, 2.0f, 0.5f);
                if(ridgedMix > 0.0f){
                    float n2 = fbm2D((wx - wxWarp) * (freq*1.25f), (wz - wzWarp) * (freq*1.25f), m_seed + 4242u, oct, 2.0f, 0.5f);
                    float ridged = 1.0f - fabsf(n2);
                    float ridgedSharp = powf(std::clamp(ridged, 0.0f, 1.0f), (biomeIdx==2 ? 0.55f : 0.85f));
                    n = (1.0f - ridgedMix)*n + ridgedMix*ridgedSharp;
                }
                if(biomeIdx == 2){
                    // Range modulation for mountains
                    const float deg = 25.0f * 3.14159265f / 180.0f;
                    float cs = cosf(deg), sn = sinf(deg);
                    float xr = wx * cs - wz * sn;
                    float zr = wx * sn + wz * cs;
                    float rangeScaleX = 1.0f / 420.0f;
                    float rangeScaleZ = 1.0f / 180.0f;
                    float rangeNoise = fbm2D(xr * rangeScaleX, zr * rangeScaleZ, m_seed + 5555u, 3, 2.0f, 0.5f);
                    float rangeRidged = 1.0f - fabsf(rangeNoise);
                    rangeRidged = powf(std::clamp(rangeRidged, 0.0f, 1.0f), 1.15f);
                    float ampBoost = 0.7f + 0.6f * rangeRidged;
                    n *= (0.8f + 0.2f * rangeRidged);
                    amp *= ampBoost;
                    float t = 0.5f*(n+1.0f);
                    t = powf(std::clamp(t, 0.0f, 1.0f), 0.7f);
                    n = t*2.0f - 1.0f;
                }
                float hf = baseHeight + amp * n;
                return {hf, grassId};
            };

            // Smooth desert factor for blending with surrounding biomes
            const float desertScale = 1.0f / 800.0f; // lower freq => larger deserts
            float dNoise = fbm2D(wx * desertScale, wz * desertScale, m_seed + 13337u, 2, 2.0f, 0.5f);
            dNoise = std::clamp(dNoise, 0.0f, 1.0f);
            // Macro mask to enforce minimum region size for deserts (prevents tiny specks)
            const float macroScale = 1.0f / 2200.0f;
            float macroN = fbm2D(wx * macroScale, wz * macroScale, m_seed + 202020u, 2, 2.0f, 0.5f);
            macroN = std::clamp(0.5f * (macroN + 1.0f), 0.0f, 1.0f);
            float macroBand = sstep(0.46f, 0.58f, macroN);
            // Prefer deserts in plains-like areas (low biome roughness)
            float plainsMask = 1.0f - sstep(0.60f, 0.88f, b);
            // Slightly relaxed inner band for a bit more coverage; multiplied by macro band for larger minimum size
            float dBand = sstep(0.66f, 0.80f, dNoise);
            float desertFactor = std::clamp(dBand * macroBand * plainsMask, 0.0f, 1.0f);

            // Base (non-desert) blended height via existing biomes
            float hBase = 0.0f; uint8_t baseSurf = 2;
            if(b < t0 - w){
                auto p = sampleHeight(0); hBase = p.first; baseSurf = p.second;
            } else if(b < t0 + w){
                auto p0 = sampleHeight(0); auto p1 = sampleHeight(1); float s = sstep(t0 - w, t0 + w, b);
                hBase = (1.0f - s) * p0.first + s * p1.first; baseSurf = (s >= 0.5f) ? p1.second : p0.second;
            } else if(b < t1 - w){
                auto p = sampleHeight(1); hBase = p.first; baseSurf = p.second;
            } else if(b < t1 + w){
                auto p0 = sampleHeight(1); auto p1 = sampleHeight(2); float s = sstep(t1 - w, t1 + w, b);
                hBase = (1.0f - s) * p0.first + s * p1.first; baseSurf = (s >= 0.5f) ? p1.second : p0.second;
            } else {
                auto p = sampleHeight(2); hBase = p.first; baseSurf = p.second;
            }
            // Desert height sample and final blend
            auto pDes = sampleHeight(3);
            float hDes = pDes.first;
            float hF = (1.0f - desertFactor) * hBase + desertFactor * hDes;
            uint8_t surfGrass = (desertFactor > 0.55f ? 9 : baseSurf);

            int height = (int)floorf(hF);
            height = std::clamp(height, 0, CHUNK_HEIGHT-1);

            // Small ponds: only in lower terrain; carve shallow bowls and fill with water (id 8)
            // Evaluate a low-frequency mask and a local radial/banded shape to decide ponds
            const float lowlandThreshold = 20.0f;
            bool allowPond = (hF < lowlandThreshold) && (desertFactor < 0.20f);
            uint8_t waterId = 8; // new voxel id for water
            int pondDepth = 0;
            if(allowPond){
                // Mask controls where ponds are allowed
                float mask = fbm2D(wx * (1.0f/180.0f), wz * (1.0f/180.0f), m_seed + 8181u, 3, 2.0f, 0.5f);
                // Center-ish when near 0.5
                float m = 1.0f - abs(mask - 0.5f) * 2.0f; // ~1 near 0.5, 0 at edges
                m = std::clamp(m, 0.0f, 1.0f);
                // Shape noise for bowl depth variation
                float shape = fbm2D(wx * (1.0f/32.0f), wz * (1.0f/32.0f), m_seed + 9191u, 2, 2.0f, 0.5f);
                shape = 0.5f * (shape + 1.0f);
                float pondChance = m * 0.35f; // keep rare
                if(pondChance > 0.32f){
                    // Carve a shallow basin of 1-3 voxels below height
                    pondDepth = 1 + (int)floorf(2.5f * shape);
                    pondDepth = std::clamp(pondDepth, 1, 3);
                }
            }

            // Meandering tan paths: banded FBM noise creates continuous contour-like lines across chunks (fade in desert)
            const float pathScale = 1.0f / 140.0f; // lower => larger features
            const int   pathBands = 5;             // number of repeating bands
            const float pathWidth = 0.035f;        // half-width in band space
            float vPath = fbm2D(wx * pathScale, wz * pathScale, m_seed + 7777u, 4, 2.0f, 0.5f);
            // Distance to band center (0.5) after tiling into bands
            float bandPhase = vPath * pathBands;
            float fracPhase = bandPhase - floorf(bandPhase);
            float dPathDist = fabsf(fracPhase - 0.5f);
            bool isPathHere = (dPathDist < pathWidth);
            // Slight widening by checking a tiny offset along two axes to avoid gaps on steep diagonals
            if(!isPathHere){
                float vP2 = fbm2D((wx+0.6f) * pathScale, (wz+0.2f) * pathScale, m_seed + 7777u, 4, 2.0f, 0.5f);
                float ph2 = vP2 * pathBands; float fr2 = ph2 - floorf(ph2);
                if(fabsf(fr2 - 0.5f) < pathWidth) isPathHere = true;
            }
            // Fade out paths progressively toward desert core
            if(desertFactor > 0.35f){
                float fadeN = fbm2D(wx * 0.05f, wz * 0.05f, m_seed + 2468u, 2, 2.0f, 0.5f);
                fadeN = 0.5f * (fadeN + 1.0f);
                float keep = 1.0f - sstep(0.35f, 0.70f, desertFactor);
                if(fadeN > keep) isPathHere = false;
            }

            for(int y = 0; y < CHUNK_HEIGHT; ++y) {
                uint8_t id = 0;
                if(y <= height) {
                    if(pondDepth > 0 && y >= height - pondDepth && y <= height){
                        // Make the bowl hollow: air up to waterline, then fill with water at or below height-1
                        // Set an approximate waterline one voxel below original surface to reduce shoreline stair-steps
                        int waterline = height - 1;
                        if(y <= waterline){
                            id = waterId; // water inside pond up to the waterline
                        } else {
                            id = 0; // hollowed air above waterline (including original surface)
                        }
                    } else if(y == height) {
                        // Surface material blending: sand in core, transition band mixes with base surface
                        if(desertFactor >= 0.60f){
                            id = 9; // sand
                        } else if(desertFactor <= 0.30f) {
                            id = (isPathHere ? 7 : surfGrass);
                        } else {
                            // Transition: stochastic mix to avoid sharp edge
                            uint32_t hh = (uint32_t)wx * 374761393u + (uint32_t)wz * 668265263u + (uint32_t)m_seed * 2246822519u;
                            hh ^= hh >> 13; hh *= 1274126177u; hh ^= hh >> 16;
                            float pr = (hh & 0xFFFF) / 65535.0f;
                            float pSand = sstep(0.30f, 0.60f, desertFactor);
                            id = (pr < pSand) ? 9 : (isPathHere ? 7 : surfGrass);
                        }
                    }
                    else if(y >= height-3) {
                        // In desert, use sand layers instead of dirt
                        id = (desertFactor > 0.55f ? 9 : 3);
                    }
                    else id = 4; // stone
                }
                c.blocks[idx(x,y,z)] = id;
            }

            // After setting terrain, add cacti across desert and transition, scaled by desert strength; very sparse overall
            if(desertFactor > 0.45f){
                // Use a repeatable hash to decide placements
                uint32_t h = (uint32_t)wx * 2166136261u ^ (uint32_t)wz * 16777619u ^ (uint32_t)m_seed;
                h ^= h >> 13; h *= 1274126177u; h ^= h >> 16;
                float chance = (h & 0xFFFF) / 65535.0f;
                float occ = sstep(0.45f, 0.95f, desertFactor);
                // Density grows toward core; square to taper fringe strongly
                if(chance < 0.0030f * occ * occ){
                    // Find surface height again (height variable is surface).
                    int baseY = height + 1; // place cactus starting above ground
                    int tall = 2 + (int)((h >> 16) % 4); // 2..5
                    for(int i=0;i<tall && baseY+i < CHUNK_HEIGHT; ++i){
                        // Only place if current is air
                        int ly = baseY + i;
                        if(c.blocks[idx(x,ly,z)] == 0) c.blocks[idx(x,ly,z)] = 10; // 10 = cactus
                    }
                    // Arms: horizontal out then vertical up (L-shape), classic saguaro style
                    auto rand01 = [&](uint32_t k){ uint32_t r=k; r ^= r>>13; r *= 1274126177u; r ^= r>>16; return (r & 0xFFFF) / 65535.0f; };
                    float armProb = 0.10f + 0.20f * occ; // ~10% at fringe up to ~30% at core per direction
                    auto tryArm = [&](int dirIdx, uint32_t salt, int dx, int dz){
                        if(rand01(h ^ salt) >= armProb) return;
                        int ay = baseY + 1 + (int)((h >> (20+dirIdx)) % std::max(1, tall-2));
                        ay = std::min(ay, CHUNK_HEIGHT-2);
                        // Fixed L-shape: go out 2, then up 1
                        const int hLen = 2;
                        const int vLen = 1;
                        int ex = x, ez = z;
                        // Horizontal segment
                        for(int i=1;i<=hLen;++i){
                            int nx = x + dx*i, nz = z + dz*i;
                            if(nx < 0 || nx >= CHUNK_SIZE || nz < 0 || nz >= CHUNK_SIZE) { break; }
                            if(c.blocks[idx(nx, ay, nz)] != 0) { break; }
                            c.blocks[idx(nx, ay, nz)] = 10;
                            ex = nx; ez = nz;
                        }
                        // Vertical segment up from tip
                        for(int j=1;j<=vLen && ay+j < CHUNK_HEIGHT; ++j){
                            if(c.blocks[idx(ex, ay+j, ez)] != 0) break;
                            c.blocks[idx(ex, ay+j, ez)] = 10;
                        }
                    };
                    // Try arms in 4 directions
                    tryArm(0, 0xA1B2C3u, +1,  0);
                    tryArm(1, 0xB2C3D4u, -1,  0);
                    tryArm(2, 0xC3D4E5u,  0, +1);
                    tryArm(3, 0xD4E5F6u,  0, -1);
                }
            }
        }
    }
}

long long VoxelWorld::Key(int cx, int cz){ return ( (long long)cx << 32 ) ^ (unsigned long long)(unsigned int)cz; }

void VoxelWorld::EnsureChunkGenerated(int cx, int cz) const {
    long long k = Key(cx, cz);
    {
        std::lock_guard<std::mutex> lk(m_mutex);
        auto it = m_chunks.find(k);
        if(it != m_chunks.end()) return;
    }
    // Generate outside the lock
    Chunk c{};
    const_cast<VoxelWorld*>(this)->GenerateChunkData(c, cx, cz);
    const_cast<VoxelWorld*>(this)->ApplyModsIfAny(c, k);
    // Insert under lock if still missing
    {
        std::lock_guard<std::mutex> lk(m_mutex);
        if(m_chunks.find(k) == m_chunks.end()){
            m_chunks.emplace(k, std::move(c));
            m_versions[k] = 0; // initial version
        }
    }
}

const Chunk& VoxelWorld::GetChunk(int cx, int cz) const {
    EnsureChunkGenerated(cx, cz);
    std::lock_guard<std::mutex> lk(m_mutex);
    return m_chunks.find(Key(cx,cz))->second;
}

Chunk& VoxelWorld::GetChunk(int cx, int cz) {
    EnsureChunkGenerated(cx, cz);
    std::lock_guard<std::mutex> lk(m_mutex);
    return m_chunks.find(Key(cx,cz))->second;
}

uint8_t VoxelWorld::GetBlock(int x, int y, int z) const {
    int cx = (int)floorf((float)x / CHUNK_SIZE);
    int cz = (int)floorf((float)z / CHUNK_SIZE);
    if(y < 0 || y >= CHUNK_HEIGHT) return 0;
    int lx = (x % CHUNK_SIZE + CHUNK_SIZE) % CHUNK_SIZE;
    int lz = (z % CHUNK_SIZE + CHUNK_SIZE) % CHUNK_SIZE;
    // Ensure exists, then read under lock to avoid invalid references
    EnsureChunkGenerated(cx, cz);
    std::lock_guard<std::mutex> lk(m_mutex);
    auto it = m_chunks.find(Key(cx, cz));
    if(it == m_chunks.end()) return 0; // should not happen
    return it->second.blocks[idx(lx, y, lz)];
}

void VoxelWorld::SetBlock(int x, int y, int z, uint8_t id) {
    int cx = (int)floorf((float)x / CHUNK_SIZE);
    int cz = (int)floorf((float)z / CHUNK_SIZE);
    if(y < 0 || y >= CHUNK_HEIGHT) return;
    int lx = (x % CHUNK_SIZE + CHUNK_SIZE) % CHUNK_SIZE;
    int lz = (z % CHUNK_SIZE + CHUNK_SIZE) % CHUNK_SIZE;
    // Ensure exists, then write under lock
    long long k = Key(cx, cz);
    EnsureChunkGenerated(cx, cz);
    std::lock_guard<std::mutex> lk(m_mutex);
    // Update mods overlay and in-memory chunk
    m_mods[k][idx(lx, y, lz)] = id;
    auto it = m_chunks.find(k);
    if(it != m_chunks.end()) it->second.blocks[idx(lx, y, lz)] = id;
    // bump version for this chunk
    m_versions[k] += 1;
}

void VoxelWorld::ApplyModsIfAny(Chunk& out, long long key) const {
    auto mit = m_mods.find(key);
    if(mit == m_mods.end()) return;
    for(const auto& kv : mit->second){
        int linear = kv.first; uint8_t id = kv.second;
        if(linear >= 0 && linear < (int)out.blocks.size()) out.blocks[linear] = id;
    }
}

unsigned VoxelWorld::GetChunkVersion(int cx, int cz) const {
    long long k = Key(cx, cz);
    std::lock_guard<std::mutex> lk(m_mutex);
    auto it = m_versions.find(k);
    return (it == m_versions.end()) ? 0u : it->second;
}

VoxelWorld::Biome VoxelWorld::GetBiomeAt(int wx, int wz) const {
    const float biomeScale = 1.0f / 256.0f;
    float b = fbm2D(wx * biomeScale, wz * biomeScale, m_seed + 9001u, 3, 2.0f, 0.5f);
    b = 0.5f * (b + 1.0f);
    b = powf(std::clamp(b, 0.0f, 1.0f), 1.2f);
    // Smooth desert factor matching terrain generation (with macro mask)
    const float desertScale = 1.0f / 800.0f;
    float dNoise = fbm2D(wx * desertScale, wz * desertScale, m_seed + 13337u, 2, 2.0f, 0.5f);
    dNoise = std::clamp(dNoise, 0.0f, 1.0f);
    auto sstep = [](float e0, float e1, float x){ float t = (x - e0) / (e1 - e0); t = std::clamp(t, 0.0f, 1.0f); return t*t*(3.0f - 2.0f*t); };
    const float macroScale = 1.0f / 2200.0f;
    float macroN = fbm2D(wx * macroScale, wz * macroScale, m_seed + 202020u, 2, 2.0f, 0.5f);
    macroN = std::clamp(0.5f * (macroN + 1.0f), 0.0f, 1.0f);
    float macroBand = sstep(0.46f, 0.58f, macroN);
    float plainsMask = 1.0f - sstep(0.60f, 0.88f, b);
    float dBand = sstep(0.66f, 0.80f, dNoise);
    float desertFactor = std::clamp(dBand * macroBand * plainsMask, 0.0f, 1.0f);
    if(desertFactor > 0.55f) return Biome::Desert;
    if(b < 0.65f) return Biome::Plains;
    if(b < 0.90f) return Biome::Hills;
    return Biome::Mountains;
}
