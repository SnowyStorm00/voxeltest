#pragma once
#include <vector>
#include <array>
#include <cstdint>
#include <DirectXMath.h>

constexpr int CHUNK_SIZE = 16;
constexpr int CHUNK_HEIGHT = 64;

struct VoxelVertex {
    DirectX::XMFLOAT3 pos;
    DirectX::XMFLOAT3 normal;
    uint32_t color;
};

struct ChunkMesh {
    std::vector<VoxelVertex> vertices;
    std::vector<uint32_t> indices;
};

struct Chunk {
    std::array<uint8_t, CHUNK_SIZE*CHUNK_SIZE*CHUNK_HEIGHT> blocks{}; // 0 = air, >0 = solid color index
};

class VoxelWorld {
public:
    VoxelWorld(int radius);
    const Chunk& GetChunk(int cx, int cz) const;
    Chunk& GetChunk(int cx, int cz);

    bool InBounds(int cx, int cz) const;

    uint8_t GetBlock(int x, int y, int z) const;
    void SetBlock(int x, int y, int z, uint8_t id);

    void GenerateFlat();
    void GenerateProcedural(unsigned seed = 1337);

private:
    int m_radius; // chunks from origin
    std::vector<Chunk> m_chunks; // (2r+1)^2 grid
};
