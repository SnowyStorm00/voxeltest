#include "voxel_mesher.h"
#include <array>
#include <vector>
#include <functional>
using namespace DirectX;

static inline uint32_t packColor(uint8_t r, uint8_t g, uint8_t b, uint8_t a) {
    // DXGI_FORMAT_R8G8B8A8_UNORM expects byte order in memory as RGBA. On little-endian CPUs,
    // the least-significant byte is stored first, so pack as 0xAABBGGRR in register yields [RR GG BB AA] in memory.
    // To avoid R/B swap, pack explicitly in little-endian RGBA order: [r][g][b][a] -> r | g<<8 | b<<16 | a<<24.
    return (uint32_t)r | (uint32_t(g) << 8) | (uint32_t(b) << 16) | (uint32_t(a) << 24);
}

static inline uint32_t colorForId(uint8_t id) {
    if(id == 2) return packColor(80,200,120, 255);      // grass
    if(id == 5) return packColor(60,160,95, 255);       // grass (hills) darker green
    if(id == 6) return packColor(45,120,75, 255);       // grass (mountains) darkest green
    if(id == 3) return packColor(155,118,83, 255);      // dirt
    if(id == 4) return packColor(150,150,160, 255);     // stone
    if(id == 7) return packColor(210, 180, 140, 255);   // path (tan)
    if(id == 8) return packColor(50, 110, 220, 255);    // water (blue)
    if(id == 9) return packColor(237, 201, 175, 255);   // sand (light desert sand)
    if(id == 10) return packColor(34, 120, 60, 255);    // cactus (green)
    return packColor(200,200,200, 255);
}

// Ray march through voxels along sun direction using 3D DDA; returns shadow factor in [0,1]
// Uses provided getBlock accessor to avoid locking the world per query.
static float VertexShadowDDA(const std::function<uint8_t(int,int,int)>& getBlock, const XMFLOAT3& p, const XMFLOAT3& sunDir, int maxSteps) {
    // Normalize direction and ensure it points from surface toward light (against light dir)
    XMVECTOR d = XMVector3Normalize(XMLoadFloat3(&sunDir));
    // March from the surface toward the sun direction; gSunDir is used for lighting as the sun-to-surface direction,
    // so use +sunDir here.
    XMFLOAT3 dir; XMStoreFloat3(&dir, d);
    // Small offset to avoid immediate self-hit; bias along normal approximated from face direction will be handled by caller via vertex offset if needed.
    const float eps = 0.001f;
    float x = p.x + dir.x * eps;
    float y = p.y + dir.y * eps;
    float z = p.z + dir.z * eps;

    // Early out: if direction is nearly zero
    float adx = fabsf(dir.x), ady = fabsf(dir.y), adz = fabsf(dir.z);
    if(adx < 1e-5f && ady < 1e-5f && adz < 1e-5f) return 1.0f;

    // Voxel coordinates
    int ix = (int)floorf(x);
    int iy = (int)floorf(y);
    int iz = (int)floorf(z);

    int stepX = (dir.x > 0) ? 1 : -1;
    int stepY = (dir.y > 0) ? 1 : -1;
    int stepZ = (dir.z > 0) ? 1 : -1;

    float nextVoxBoundaryX = (stepX > 0) ? (float)(ix + 1) : (float)ix;
    float nextVoxBoundaryY = (stepY > 0) ? (float)(iy + 1) : (float)iy;
    float nextVoxBoundaryZ = (stepZ > 0) ? (float)(iz + 1) : (float)iz;

    float txMax = (adx > 1e-6f) ? (nextVoxBoundaryX - x) / dir.x : FLT_MAX;
    float tyMax = (ady > 1e-6f) ? (nextVoxBoundaryY - y) / dir.y : FLT_MAX;
    float tzMax = (adz > 1e-6f) ? (nextVoxBoundaryZ - z) / dir.z : FLT_MAX;

    float txDelta = (adx > 1e-6f) ? (float)stepX / dir.x : FLT_MAX; // how far along t to cross a cell in X
    float tyDelta = (ady > 1e-6f) ? (float)stepY / dir.y : FLT_MAX;
    float tzDelta = (adz > 1e-6f) ? (float)stepZ / dir.z : FLT_MAX;

    // March up to N steps or until we exit plausible world bounds
    const int MaxSteps = (maxSteps > 0 ? maxSteps : 1);
    for(int i=0;i<MaxSteps;++i){
        // Move to next voxel boundary along the smallest t
        if(txMax < tyMax && txMax < tzMax){
            ix += stepX; txMax += txDelta;
        } else if(tyMax < tzMax){
            iy += stepY; tyMax += tyDelta;
        } else {
            iz += stepZ; tzMax += tzDelta;
        }
        // If below world, treat as blocked (terrain goes to -inf effectively)
        if(iy < 0) return 0.0f;
        if(iy >= CHUNK_HEIGHT) return 1.0f; // above world, unblocked
        // Check occupancy
    if(getBlock(ix, iy, iz) != 0){
            return 0.0f; // in shadow
        }
    }
    // No hit within range; consider lit
    return 1.0f;
}

// Greedy meshing across all 3 axes. Merges adjacent faces with same material and facing.
ChunkMesh VoxelMesher::BuildChunkMesh(const VoxelWorld& world, int cx, int cz, const DirectX::XMFLOAT3& sunDir, const MesherSettings& settings) {
    ChunkMesh mesh;
    const int baseX = cx * CHUNK_SIZE;
    const int baseZ = cz * CHUNK_SIZE;
    const int sizeX = CHUNK_SIZE;
    const int sizeY = CHUNK_HEIGHT;
    const int sizeZ = CHUNK_SIZE;

    struct MaskCell { uint8_t id; bool back; uint8_t lit; }; // back=false => +axis normal, back=true => -axis normal; lit: 1=lit,0=shadow

    // ---- Local 3x3 chunk cache to avoid world locks during meshing ----
    const int cminX = cx - 1, cmaxX = cx + 1;
    const int cminZ = cz - 1, cmaxZ = cz + 1;
    std::vector<std::vector<uint8_t>> cache(9); // [dx+1 + 3*(dz+1)]
    auto keyIdx = [&](int ccx, int ccz){ return (ccx - cminX) + 3*(ccz - cminZ); };
    for(int ccz = cminZ; ccz <= cmaxZ; ++ccz){
        for(int ccx = cminX; ccx <= cmaxX; ++ccx){
            const Chunk& ch = world.GetChunk(ccx, ccz);
            auto &v = cache[keyIdx(ccx, ccz)];
            v.assign(ch.blocks.begin(), ch.blocks.end());
        }
    }
    auto idx3 = [](int x,int y,int z){ return x + CHUNK_SIZE*(z + CHUNK_SIZE*y); };
    auto getBlock = [&](int wx, int y, int wz)->uint8_t{
        if(y < 0 || y >= CHUNK_HEIGHT) return 0;
        int ccx = (int)floorf((float)wx / CHUNK_SIZE);
        int ccz = (int)floorf((float)wz / CHUNK_SIZE);
        if(ccx < cminX || ccx > cmaxX || ccz < cminZ || ccz > cmaxZ) return 0; // treat outside as air
        int lx = (wx % CHUNK_SIZE + CHUNK_SIZE) % CHUNK_SIZE;
        int lz = (wz % CHUNK_SIZE + CHUNK_SIZE) % CHUNK_SIZE;
        const auto &v = cache[keyIdx(ccx, ccz)];
        return v[idx3(lx,y,lz)];
    };

    auto faceCenterAndNormal = [&](int axis, int k, int u, int v, bool back){
        XMFLOAT3 p; XMFLOAT3 n;
        if(axis==0){ // X faces, u->Z, v->Y
            float x = float(baseX + (back ? (k-1) : k));
            p = { x, float(v) + 0.5f, float(baseZ + u) + 0.5f };
            n = back ? XMFLOAT3{-1,0,0} : XMFLOAT3{1,0,0};
        } else if(axis==1){ // Y faces, u->X, v->Z
            float y = float(k);
            p = { float(baseX + u) + 0.5f, y, float(baseZ + v) + 0.5f };
            n = back ? XMFLOAT3{0,-1,0} : XMFLOAT3{0,1,0};
        } else { // Z faces, u->X, v->Y
            float z = float(baseZ + (back ? (k-1) : k));
            p = { float(baseX + u) + 0.5f, float(v) + 0.5f, z };
            n = back ? XMFLOAT3{0,0,-1} : XMFLOAT3{0,0,1};
        }
        return std::pair<XMFLOAT3,XMFLOAT3>(p,n);
    };

    // Utility to emit a quad given axis, slice k, rectangle [u0,u1) x [v0,v1) in (other two) axes
    std::function<void(int,int,int,int,int,int,uint8_t,bool,int)> emitQuad;
    emitQuad = [&](int axis, int k, int u0, int v0, int u1, int v1, uint8_t id, bool back, int depth){
        uint32_t baseColor = colorForId(id);
        XMFLOAT3 n;
        // axis: 0=X, 1=Y, 2=Z ; back=false => +axis, back=true => -axis
        if(axis==0) n = back ? XMFLOAT3{-1,0,0} : XMFLOAT3{1,0,0};
        else if(axis==1) n = back ? XMFLOAT3{0,-1,0} : XMFLOAT3{0,1,0};
        else n = back ? XMFLOAT3{0,0,-1} : XMFLOAT3{0,0,1};

    // Compute 4 corners
    XMFLOAT3 p[4];
        if(axis==0) {
            float x = float(baseX + k);
            float y0 = float(v0), y1 = float(v1);
            float z0 = float(baseZ + u0), z1 = float(baseZ + u1);
            if(!back){ // +X
                p[0] = {x, y0, z1}; p[1] = {x, y1, z1}; p[2] = {x, y1, z0}; p[3] = {x, y0, z0};
            } else {   // -X
                p[0] = {x, y0, z0}; p[1] = {x, y1, z0}; p[2] = {x, y1, z1}; p[3] = {x, y0, z1};
            }
        } else if(axis==1) {
            float y = float(k);
            float x0 = float(baseX + u0), x1 = float(baseX + u1);
            float z0 = float(baseZ + v0), z1 = float(baseZ + v1);
            if(!back){ // +Y
                p[0] = {x0, y, z1}; p[1] = {x0, y, z0}; p[2] = {x1, y, z0}; p[3] = {x1, y, z1};
            } else {   // -Y
                p[0] = {x0, y, z0}; p[1] = {x0, y, z1}; p[2] = {x1, y, z1}; p[3] = {x1, y, z0};
            }
        } else { // axis==2
            float z = float(baseZ + k);
            float x0 = float(baseX + u0), x1 = float(baseX + u1);
            float y0 = float(v0), y1 = float(v1);
            if(!back){ // +Z
                p[0] = {x0, y0, z}; p[1] = {x0, y1, z}; p[2] = {x1, y1, z}; p[3] = {x1, y0, z};
            } else {   // -Z
                p[0] = {x1, y0, z}; p[1] = {x1, y1, z}; p[2] = {x0, y1, z}; p[3] = {x0, y0, z};
            }
        }

    // Shadow factor sampling helpers
        auto biasNS = [&](const XMFLOAT3& q){
            // Bias along normal and slightly toward the sun to avoid hitting the source voxel edge
            XMFLOAT3 out{ q.x + n.x*settings.normalBias + sunDir.x*settings.sunBias,
                          q.y + n.y*settings.normalBias + sunDir.y*settings.sunBias,
                          q.z + n.z*settings.normalBias + sunDir.z*settings.sunBias };
            return out;
        };
        // Build two in-plane axes for sampling
        XMFLOAT3 t1, t2;
        if(axis==1){ t1 = XMFLOAT3{1,0,0}; t2 = XMFLOAT3{0,0,1}; }
        else if(axis==0){ t1 = XMFLOAT3{0,1,0}; t2 = XMFLOAT3{0,0,1}; }
        else { t1 = XMFLOAT3{1,0,0}; t2 = XMFLOAT3{0,1,0}; }
        auto addScaled = [](const XMFLOAT3& a, const XMFLOAT3& b, float s){ return XMFLOAT3{ a.x + b.x*s, a.y + b.y*s, a.z + b.z*s }; };
    auto avgShadow = [&](const XMFLOAT3& vtx){
            // Slightly inset toward face center
            XMFLOAT3 center{ (p[0].x + p[2].x)*0.5f, (p[0].y + p[2].y)*0.5f, (p[0].z + p[2].z)*0.5f };
            XMFLOAT3 q = XMFLOAT3{ vtx.x + (center.x - vtx.x)*settings.inset,
                                    vtx.y + (center.y - vtx.y)*settings.inset,
                                    vtx.z + (center.z - vtx.z)*settings.inset };
            const float off = settings.sampleRadius; // plane offsets in voxel units
            // Generate samples in a simple grid: if vertexSamples==4 use corners, if 8 include edges, if 1 just center
            float sum = 0.0f; int count = 0;
            auto sampleAt = [&](float su, float sv){
                XMFLOAT3 pos = addScaled(addScaled(q, t1, su), t2, sv);
                sum += VertexShadowDDA(getBlock, biasNS(pos), sunDir, settings.maxSteps); ++count;
            };
            if(settings.vertexSamples <= 1){
                sampleAt(0.0f, 0.0f);
            } else if(settings.vertexSamples <= 2){
                sampleAt( off,  off);
                sampleAt(-off, -off);
            } else if(settings.vertexSamples <= 4){
                sampleAt( off,  off);
                sampleAt(-off,  off);
                sampleAt( off, -off);
                sampleAt(-off, -off);
            } else {
                // 8 samples: 4 corners + 4 edge midpoints
                sampleAt( off,  off); sampleAt(-off,  off); sampleAt( off, -off); sampleAt(-off, -off);
                sampleAt( off, 0.0f); sampleAt(-off, 0.0f); sampleAt(0.0f,  off); sampleAt(0.0f, -off);
            }
            return (count>0)? (sum / (float)count) : 1.0f;
        };
    auto shadowAtPoint = [&](const XMFLOAT3& q){ return VertexShadowDDA(getBlock, biasNS(q), sunDir, settings.maxSteps); };
        auto randomFloat = [](uint32_t& state){ state ^= state << 13; state ^= state >> 17; state ^= state << 5; return (state & 0xFFFFFF) / 16777215.0f; };
        auto sampleSoftShadowAt = [&](const XMFLOAT3& q){
            if(settings.hardShadows || settings.sunSamples <= 1) return shadowAtPoint(q);
            // Build orthonormal basis around sunDir for cone sampling
            XMVECTOR d = XMVector3Normalize(XMLoadFloat3(&sunDir));
            XMVECTOR up = XMVectorSet(0,1,0,0);
            if(fabsf(XMVectorGetX(XMVector3Dot(d, up))) > 0.99f) up = XMVectorSet(1,0,0,0);
            XMVECTOR t1v = XMVector3Normalize(XMVector3Cross(up, d));
            XMVECTOR t2v = XMVector3Cross(d, t1v);
            XMFLOAT3 t1, t2, d3; XMStoreFloat3(&t1, t1v); XMStoreFloat3(&t2, t2v); XMStoreFloat3(&d3, d);
            uint32_t rng = (uint32_t)(k*73856093u) ^ (uint32_t)(u0*19349663u) ^ (uint32_t)(v0*83492791u);
            float sum = 0.0f; int cnt = 0;
            for(int i=0;i<settings.sunSamples;++i){
                float u = randomFloat(rng); float v = randomFloat(rng);
                // Uniform disk sample then map to small cone
                float r = sqrtf(u) * settings.sunAngularRadius;
                float a = 6.2831853f * v;
                float dx = r * cosf(a);
                float dy = r * sinf(a);
                XMFLOAT3 dir{ d3.x + t1.x*dx + t2.x*dy, d3.y + t1.y*dx + t2.y*dy, d3.z + t1.z*dx + t2.z*dy };
                XMVECTOR dn = XMVector3Normalize(XMLoadFloat3(&dir)); XMFLOAT3 dirN; XMStoreFloat3(&dirN, dn);
                // Temporarily override sun direction for this sample
                XMFLOAT3 qbiased = biasNS(q);
                sum += VertexShadowDDA(getBlock, qbiased, dirN, settings.maxSteps);
                ++cnt;
            }
            return (cnt>0)? sum / (float)cnt : shadowAtPoint(q);
        };
        // Fast path for initial load: no shadowing, no subdivision
        if(settings.fastBakeOnly){
            uint32_t aFull = 255;
            auto withA = [&](uint32_t rgb){ return (rgb & 0x00FFFFFFu) | (aFull<<24); };
            uint32_t c0 = withA(baseColor), c1 = c0, c2 = c0, c3 = c0;
            uint32_t start = (uint32_t)mesh.vertices.size();
            mesh.vertices.push_back({p[0], n, c0});
            mesh.vertices.push_back({p[1], n, c1});
            mesh.vertices.push_back({p[2], n, c2});
            mesh.vertices.push_back({p[3], n, c3});
            mesh.indices.push_back(start + 0);
            mesh.indices.push_back(start + 1);
            mesh.indices.push_back(start + 2);
            mesh.indices.push_back(start + 0);
            mesh.indices.push_back(start + 2);
            mesh.indices.push_back(start + 3);
            return;
        }

        // Corner samples (kept for variance check)
        float s0 = avgShadow(p[0]);
        float s1 = avgShadow(p[1]);
        float s2 = avgShadow(p[2]);
        float s3 = avgShadow(p[3]);
        // Center and edge midpoints for better split decisions
        XMFLOAT3 center{ (p[0].x + p[2].x)*0.5f, (p[0].y + p[2].y)*0.5f, (p[0].z + p[2].z)*0.5f };
        float sC = shadowAtPoint(center);
        XMFLOAT3 eTop   { (p[0].x + p[1].x)*0.5f, (p[0].y + p[1].y)*0.5f, (p[0].z + p[1].z)*0.5f };
        XMFLOAT3 eRight { (p[1].x + p[2].x)*0.5f, (p[1].y + p[2].y)*0.5f, (p[1].z + p[2].z)*0.5f };
        XMFLOAT3 eBottom{ (p[2].x + p[3].x)*0.5f, (p[2].y + p[3].y)*0.5f, (p[2].z + p[3].z)*0.5f };
        XMFLOAT3 eLeft  { (p[3].x + p[0].x)*0.5f, (p[3].y + p[0].y)*0.5f, (p[3].z + p[0].z)*0.5f };
        float sT = shadowAtPoint(eTop);
        float sR = shadowAtPoint(eRight);
        float sB = shadowAtPoint(eBottom);
        float sL = shadowAtPoint(eLeft);
        auto withAlpha = [&](uint32_t rgb, float s){
            // Clamp with a small minimum to preserve some detail but allow fairly dark shadows
            s = std::min(1.0f, std::max(settings.minShadow, s));
            uint8_t a = (uint8_t)std::round(s * 255.0f);
            // Our packColor stores bytes in memory as [R][G][B][A] (little-endian),
            // and the integer layout is (A<<24 | B<<16 | G<<8 | R).
            // Extract channels accordingly without swapping.
            uint8_t r = (rgb) & 0xFF;
            uint8_t g = (rgb >> 8) & 0xFF;
            uint8_t b = (rgb >> 16) & 0xFF;
            return packColor(r,g,b,a);
        };
    // Default: make the entire quad a single uniform shadow value (prevents gradients).
        // Hard shadows: conservative classification using min across several points on the quad
        float sUniform;
        if(settings.hardShadows){
            float sMin = std::min(std::min(std::min(std::min(std::min(s0,s1),std::min(s2,s3)), sC), std::min(sT,sB)), std::min(sL,sR));
            sUniform = (sMin > settings.mergeLitThreshold) ? 1.0f : settings.minShadow;
        } else {
            // Soft shadows: sample an area sun at the face center
            sUniform = sampleSoftShadowAt(center);
        }
    // Water should not receive baked shadow darkening; keep it bright for readability
    if(id == 8) sUniform = 1.0f;
    uint32_t c0 = withAlpha(baseColor, sUniform);
    uint32_t c1 = withAlpha(baseColor, sUniform);
    uint32_t c2 = withAlpha(baseColor, sUniform);
    uint32_t c3 = withAlpha(baseColor, sUniform);

    // Optional subdivision if shadow varies a lot across this rect
    if(settings.enableSubdivision && depth < settings.subdivMaxDepth){
        float smin = std::min(std::min(std::min(std::min(std::min(s0,s1),std::min(s2,s3)), sC), std::min(sT,sB)), std::min(sL,sR));
        float smax = std::max(std::max(std::max(std::max(std::max(s0,s1),std::max(s2,s3)), sC), std::max(sT,sB)), std::max(sL,sR));
            if((smax - smin) > settings.subdivThreshold){
                int du = u1 - u0;
                int dv = v1 - v0;
                if(du > 1 || dv > 1){
            // Choose split direction by which half differs more across edges
            float diffU = fabsf(sL - sR);
            float diffV = fabsf(sT - sB);
            if((diffU >= diffV && du > 1) || (du > 1 && dv <= 1)){
                        int um = (u0 + u1) / 2;
                        if(um <= u0) um = u0 + 1;
                        if(um >= u1) um = u1 - 1;
                        emitQuad(axis,k,u0,v0,um,v1,id,back,depth+1);
                        emitQuad(axis,k,um,v0,u1,v1,id,back,depth+1);
                        return;
            } else if(dv > 1){
                        int vm = (v0 + v1) / 2;
                        if(vm <= v0) vm = v0 + 1;
                        if(vm >= v1) vm = v1 - 1;
                        emitQuad(axis,k,u0,v0,u1,vm,id,back,depth+1);
                        emitQuad(axis,k,u0,vm,u1,v1,id,back,depth+1);
                        return;
                    }
                }
            }
        }

        uint32_t start = (uint32_t)mesh.vertices.size();
        mesh.vertices.push_back({p[0], n, c0});
        mesh.vertices.push_back({p[1], n, c1});
        mesh.vertices.push_back({p[2], n, c2});
        mesh.vertices.push_back({p[3], n, c3});
        // CCW triangles
        mesh.indices.push_back(start + 0);
        mesh.indices.push_back(start + 1);
        mesh.indices.push_back(start + 2);
        mesh.indices.push_back(start + 0);
        mesh.indices.push_back(start + 2);
        mesh.indices.push_back(start + 3);
    };

    // For each axis, build face mask per slice and greedy merge
    auto processAxis = [&](int axis){
        int du0, dv0, du1, dv1; // dimensions in mask
        int kMax; // number of slices (inclusive planes)
        // Axis-aligned sizes and accessors
    auto sample = [&](int x, int y, int z){ return getBlock(x,y,z); };

        if(axis==0){ du0 = sizeZ; dv0 = sizeY; kMax = sizeX; }
        else if(axis==1){ du0 = sizeX; dv0 = sizeZ; kMax = sizeY; }
        else { du0 = sizeX; dv0 = sizeY; kMax = sizeZ; }

        std::vector<MaskCell> mask(du0 * dv0);

        for(int k=0; k<=kMax; ++k){
            // Build mask between layer k-1 and k
            for(int v=0; v<dv0; ++v){
                for(int u=0; u<du0; ++u){
                    uint8_t a=0, b=0; // a: block on negative side, b: block on positive side
                    if(axis==0){
                        int x0 = baseX + (k-1);
                        int x1 = baseX + k;
                        int y = v;
                        int z = baseZ + u;
                        a = sample(x0, y, z);
                        b = sample(x1, y, z);
                    } else if(axis==1){
                        int y0 = k-1;
                        int y1 = k;
                        int x = baseX + u;
                        int z = baseZ + v;
                        a = sample(x, y0, z);
                        b = sample(x, y1, z);
                    } else { // axis==2
                        int z0 = baseZ + (k-1);
                        int z1 = baseZ + k;
                        int x = baseX + u;
                        int y = v;
                        a = sample(x, y, z0);
                        b = sample(x, y, z1);
                    }
                    MaskCell cell{0,false,1};
                    if(a!=0 && b==0){
                        // +axis face
                        auto pn = faceCenterAndNormal(axis, k, u, v, false);
                        XMFLOAT3 pc = pn.first; XMFLOAT3 n = pn.second;
                        XMFLOAT3 pb{ pc.x + n.x*settings.normalBias, pc.y + n.y*settings.normalBias, pc.z + n.z*settings.normalBias };
                        if(settings.fastBakeOnly){
                            cell = {a,false,1};
                        } else {
                            float s = VertexShadowDDA(getBlock, pb, sunDir, settings.maxSteps);
                            cell = {a,false,(uint8_t)((s>settings.mergeLitThreshold)?1:0)};
                        }
                    } else if(a==0 && b!=0){
                        // -axis face
                        auto pn = faceCenterAndNormal(axis, k, u, v, true);
                        XMFLOAT3 pc = pn.first; XMFLOAT3 n = pn.second;
                        XMFLOAT3 pb{ pc.x + n.x*settings.normalBias, pc.y + n.y*settings.normalBias, pc.z + n.z*settings.normalBias };
                        if(settings.fastBakeOnly){
                            cell = {b,true,1};
                        } else {
                            float s = VertexShadowDDA(getBlock, pb, sunDir, settings.maxSteps);
                            cell = {b,true,(uint8_t)((s>settings.mergeLitThreshold)?1:0)};
                        }
                    }
                    mask[u + du0*v] = cell;
                }
            }

            // Greedy merge on mask (u along width du0, v along height dv0)
            int u = 0, v = 0;
            std::vector<char> used(du0*dv0, 0);
            for(int v0=0; v0<dv0; ++v0){
                for(int u0=0; u0<du0; ++u0){
                    int idx0 = u0 + du0*v0;
                    if(used[idx0]) continue;
                    MaskCell m0 = mask[idx0];
                    if(m0.id==0) { used[idx0]=1; continue; }
                    // Compute width
                    int w=1;
                    while(u0 + w < du0){
                        MaskCell m1 = mask[(u0+w) + du0*v0];
                        if(m1.id==m0.id && m1.back==m0.back && m1.lit==m0.lit && !used[(u0+w)+du0*v0]) w++; else break;
                    }
                    // Compute height
                    int h=1;
                    while(v0 + h < dv0){
                        bool rowMatch = true;
                        for(int uu=0; uu<w; ++uu){
                            int i2 = (u0+uu) + du0*(v0+h);
                            MaskCell m2 = mask[i2];
                            if(!(m2.id==m0.id && m2.back==m0.back && m2.lit==m0.lit && !used[i2])) { rowMatch=false; break; }
                        }
                        if(!rowMatch) break;
                        h++;
                    }
                    // Mark used
                    for(int dv=0; dv<h; ++dv){
                        for(int du=0; du<w; ++du){
                            used[(u0+du) + du0*(v0+dv)] = 1;
                        }
                    }
                    // Emit, carrying the mask's binary lit flag to force uniform shading on this quad
                    if(axis==0) emitQuad(axis, k, u0, v0, u0+w, v0+h, m0.id, m0.back, 0);
                    else if(axis==1) emitQuad(axis, k, u0, v0, u0+w, v0+h, m0.id, m0.back, 0);
                    else emitQuad(axis, k, u0, v0, u0+w, v0+h, m0.id, m0.back, 0);
                }
            }
        }
    };

    processAxis(0);
    processAxis(1);
    processAxis(2);

    return mesh;
}
