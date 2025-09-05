#pragma once
#include <vector>
#include <array>
#include <cstdint>
#include <unordered_map>
#include <mutex>
#include <DirectXMath.h>

constexpr int CHUNK_SIZE = 16;
constexpr int CHUNK_HEIGHT = 128;

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
    enum class Biome : uint8_t { Plains=0, Hills=1, Mountains=2 };
    explicit VoxelWorld(unsigned seed = 1337);
    const Chunk& GetChunk(int cx, int cz) const;
    Chunk& GetChunk(int cx, int cz);

    // Infinite world: always in bounds
    bool InBounds(int, int) const { return true; }

    uint8_t GetBlock(int x, int y, int z) const;
    void SetBlock(int x, int y, int z, uint8_t id);

    // Query biome at world XZ position
    Biome GetBiomeAt(int wx, int wz) const;

    // Versioning for chunks: increments when blocks are modified in a chunk
    unsigned GetChunkVersion(int cx, int cz) const;

private:
    static long long Key(int cx, int cz);
    void EnsureChunkGenerated(int cx, int cz) const;
    void GenerateChunkData(Chunk& out, int cx, int cz) const;
    void ApplyModsIfAny(Chunk& out, long long key) const;

    unsigned m_seed = 1337;
    mutable std::unordered_map<long long, Chunk> m_chunks; // dynamic storage
    // Modified voxels overlay (persisted separately from procedural base)
    mutable std::unordered_map<long long, std::unordered_map<int, uint8_t>> m_mods;
    // Per-chunk version for invalidating cached meshes, etc.
    mutable std::unordered_map<long long, unsigned> m_versions;
    mutable std::mutex m_mutex;
};
