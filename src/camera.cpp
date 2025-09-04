#include "camera.h"
using namespace DirectX;

void Camera::SetLens(float fovY, float aspect, float zn, float zf) {
    XMStoreFloat4x4(&m_proj, XMMatrixPerspectiveFovLH(fovY, aspect, zn, zf));
}

void Camera::LookAt(FXMVECTOR pos, FXMVECTOR target, FXMVECTOR up) {
    XMStoreFloat3(&m_pos, pos);
    auto L = XMVector3Normalize(XMVectorSubtract(target, pos));
    auto R = XMVector3Normalize(XMVector3Cross(up, L));
    auto U = XMVector3Cross(L, R);
    XMStoreFloat3(&m_right, R);
    XMStoreFloat3(&m_up, U);
    XMStoreFloat3(&m_look, L);
    UpdateViewMatrix();
}

void Camera::Strafe(float d) {
    XMVECTOR s = XMVectorReplicate(d);
    XMVECTOR r = XMLoadFloat3(&m_right);
    XMVECTOR p = XMLoadFloat3(&m_pos);
    XMStoreFloat3(&m_pos, XMVectorAdd(p, XMVectorMultiply(s, r)));
}

void Camera::Walk(float d) {
    XMVECTOR s = XMVectorReplicate(d);
    XMVECTOR l = XMLoadFloat3(&m_look);
    XMVECTOR p = XMLoadFloat3(&m_pos);
    XMStoreFloat3(&m_pos, XMVectorAdd(p, XMVectorMultiply(s, l)));
}

void Camera::Pitch(float angle) {
    XMMATRIX R = XMMatrixRotationAxis(XMLoadFloat3(&m_right), angle);
    XMStoreFloat3(&m_up, XMVector3TransformNormal(XMLoadFloat3(&m_up), R));
    XMStoreFloat3(&m_look, XMVector3TransformNormal(XMLoadFloat3(&m_look), R));
}

void Camera::RotateY(float angle) {
    XMMATRIX R = XMMatrixRotationY(angle);
    XMStoreFloat3(&m_right, XMVector3TransformNormal(XMLoadFloat3(&m_right), R));
    XMStoreFloat3(&m_up, XMVector3TransformNormal(XMLoadFloat3(&m_up), R));
    XMStoreFloat3(&m_look, XMVector3TransformNormal(XMLoadFloat3(&m_look), R));
}

void Camera::UpdateViewMatrix() {
    XMVECTOR R = XMLoadFloat3(&m_right);
    XMVECTOR U = XMLoadFloat3(&m_up);
    XMVECTOR L = XMLoadFloat3(&m_look);
    XMVECTOR P = XMLoadFloat3(&m_pos);

    L = XMVector3Normalize(L);
    U = XMVector3Normalize(XMVector3Cross(L, R));
    R = XMVector3Cross(U, L);

    float x = -XMVectorGetX(XMVector3Dot(P, R));
    float y = -XMVectorGetX(XMVector3Dot(P, U));
    float z = -XMVectorGetX(XMVector3Dot(P, L));

    m_view = {
        XMVectorGetX(R), XMVectorGetX(U), XMVectorGetX(L), 0.0f,
        XMVectorGetY(R), XMVectorGetY(U), XMVectorGetY(L), 0.0f,
        XMVectorGetZ(R), XMVectorGetZ(U), XMVectorGetZ(L), 0.0f,
        x, y, z, 1.0f
    };
}
