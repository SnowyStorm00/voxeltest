#pragma once
#include <vector>
#include <string>
#include <cstdint>
#include <DirectXMath.h>

struct ModelVertex {
    DirectX::XMFLOAT3 pos;
    DirectX::XMFLOAT3 normal;
};

struct ModelMesh {
    std::vector<ModelVertex> vertices;
    std::vector<uint32_t> indices;
};

// Loads a very simple Wavefront OBJ (positions, normals; triangulates polygons).
// - Supports lines like: v x y z, vn x y z, f v//vn or v/vt/vn or v// or v formats
// - Ignores materials and texcoords
// - Returns true on success.
bool LoadOBJ(const std::wstring& path, ModelMesh& out);
