#pragma once
#include <DirectXMath.h>

class Camera {
public:
    void SetLens(float fovY, float aspect, float zn, float zf);
    void LookAt(DirectX::FXMVECTOR pos, DirectX::FXMVECTOR target, DirectX::FXMVECTOR up);
    void Strafe(float d);
    void Walk(float d);
    void Pitch(float angle);
    void RotateY(float angle);

    void UpdateViewMatrix();

    DirectX::XMMATRIX View() const { return DirectX::XMLoadFloat4x4(&m_view); }
    DirectX::XMMATRIX Proj() const { return DirectX::XMLoadFloat4x4(&m_proj); }
    DirectX::XMMATRIX ViewProj() const { return DirectX::XMMatrixMultiply(View(), Proj()); }

    DirectX::XMFLOAT3 Position() const { return m_pos; }
    void SetPosition(float x, float y, float z) { m_pos = {x,y,z}; }
    DirectX::XMFLOAT3 Look() const { return m_look; }
    DirectX::XMFLOAT3 Right() const { return m_right; }

private:
    DirectX::XMFLOAT3 m_pos {0, 0, -5};
    DirectX::XMFLOAT3 m_right {1, 0, 0};
    DirectX::XMFLOAT3 m_up {0, 1, 0};
    DirectX::XMFLOAT3 m_look {0, 0, 1};

    DirectX::XMFLOAT4X4 m_view {};
    DirectX::XMFLOAT4X4 m_proj {};
};
