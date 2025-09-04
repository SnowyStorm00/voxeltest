#include "d3d_app.h"
#include "camera.h"
#include "voxel_world.h"
#include "voxel_mesher.h"
#include <d3d11.h>
#include <d3dcompiler.h>
#include <DirectXMath.h>
#include <wrl/client.h>
#include <vector>
#include <memory>
#include <array>
#include <windowsx.h>
#include <unordered_map>
#include <unordered_set>
#include <algorithm>
#include <cmath>
#include <sstream>
#include <string>

using namespace DirectX;
using Microsoft::WRL::ComPtr;

struct VSConstants {
    XMMATRIX viewProj;
};

class VoxelApp : public D3DApp {
public:
    VoxelApp(HINSTANCE hInst) : D3DApp(hInst, 1280, 720, L"Voxel RPG") {}

public:
    bool Initialize() override {
        if(!D3DApp::Initialize()) return false;

    // Camera
    m_cam.SetLens(XM_PIDIV4, float(m_width)/float(m_height), 0.1f, 1000.0f);
    m_cam.LookAt(XMVectorSet(8,12,-20,1), XMVectorSet(8,8,8,1), XMVectorSet(0,1,0,0));
    m_cam.UpdateViewMatrix();

        // World
    m_world = std::make_unique<VoxelWorld>(m_viewRadius);
    m_world->GenerateProcedural(1337);

    // Seed chunk streaming queue based on current camera position
    UpdateCameraChunk();
    EnqueueVisibleChunks();

        // Create shaders and input layout
        if(!CreatePipeline()) return false;

        // Constant buffer
        D3D11_BUFFER_DESC cbd{};
        cbd.BindFlags = D3D11_BIND_CONSTANT_BUFFER;
        cbd.ByteWidth = sizeof(VSConstants);
        cbd.Usage = D3D11_USAGE_DEFAULT;
        GetDevice()->CreateBuffer(&cbd, nullptr, m_cbVP.GetAddressOf());

        // Rasterizer
    D3D11_RASTERIZER_DESC rs{}; rs.CullMode = D3D11_CULL_BACK; rs.FillMode = D3D11_FILL_SOLID; rs.DepthClipEnable = TRUE; rs.FrontCounterClockwise = TRUE;
        GetDevice()->CreateRasterizerState(&rs, m_rsSolid.GetAddressOf());

    D3D11_DEPTH_STENCIL_DESC dsd{}; dsd.DepthEnable = TRUE; dsd.DepthWriteMask = D3D11_DEPTH_WRITE_MASK_ALL; dsd.DepthFunc = D3D11_COMPARISON_LESS_EQUAL;
        GetDevice()->CreateDepthStencilState(&dsd, m_dss.GetAddressOf());

    // Overlay pipeline (2D quads)
    if(!CreateOverlayPipeline()) return false;

    // Initial title so users see something immediately
    SetWindowTextW(GetHwnd(), L"Voxel RPG - starting...");
    return true;
    }

    void Update(float dt) override {
        m_lastDt = dt;
        // WASD movement (+Space/Ctrl for vertical)
    auto isDown = [](int vk){ return (GetAsyncKeyState(vk) & 0x8000) != 0; };
    bool sprint = m_keys[VK_SHIFT] || isDown(VK_LSHIFT) || isDown(VK_RSHIFT);
        float speed = m_moveSpeed * (sprint ? 3.0f : 1.0f);
    if(m_keys['W'] || isDown('W')) m_cam.Walk( speed * dt);
    if(m_keys['S'] || isDown('S')) m_cam.Walk(-speed * dt);
    if(m_keys['A'] || isDown('A')) m_cam.Strafe(-speed * dt);
    if(m_keys['D'] || isDown('D')) m_cam.Strafe( speed * dt);

        // Vertical move in world space
        {
            auto p = m_cam.Position();
            if(m_keys[VK_SPACE] || isDown(VK_SPACE)) p.y += speed * dt;
            if(m_keys[VK_CONTROL] || isDown(VK_LCONTROL) || isDown(VK_RCONTROL)) p.y -= speed * dt;
            m_cam.SetPosition(p.x, p.y, p.z);
        }

        m_cam.UpdateViewMatrix();

        // Update streaming center when crossing chunk boundaries
        int prevCX = m_camChunkX, prevCZ = m_camChunkZ;
        UpdateCameraChunk();
        if(prevCX != m_camChunkX || prevCZ != m_camChunkZ) {
            EnqueueVisibleChunks();
            CleanupFarChunks();
        }

        // Build a few pending chunks per frame to avoid stalls
    int built = 0;
    const int MaxPerFrame = 64;
        while(built < MaxPerFrame && !m_pending.empty()) {
            auto [cx, cz] = m_pending.back(); m_pending.pop_back();
            if(!IsWithinView(cx, cz)) continue;
            long long key = ChunkKey(cx, cz);
            if(m_chunks.find(key) != m_chunks.end()) continue;
            GPUChunk g{};
            if(BuildChunkGPU(cx, cz, g)) {
                m_chunks.emplace(key, std::move(g));
            }
            ++built;
        }

    // FPS counter: instantaneous (current frame)
    m_fpsAccum += dt;
    m_fpsFrames += 1;
    float fpsNow = (m_lastDt > 1e-6f) ? (1.0f / m_lastDt) : 0.0f;
        {
            std::wstringstream ss;
            ss.setf(std::ios::fixed); ss.precision(1);
            ss << L"Voxel RPG - " << fpsNow << L" FPS";
            if(sprint) ss << L"  (x3)";
            SetWindowTextW(GetHwnd(), ss.str().c_str());
        }
    }

    void Render() override {
        float clear[4] = {0.53f, 0.81f, 0.92f, 1.0f};
        GetContext()->ClearRenderTargetView(GetRTV(), clear);
        GetContext()->ClearDepthStencilView(GetDSV(), D3D11_CLEAR_DEPTH|D3D11_CLEAR_STENCIL, 1.0f, 0);

        GetContext()->OMSetDepthStencilState(m_dss.Get(), 0);
        GetContext()->RSSetState(m_rsSolid.Get());
        GetContext()->IASetInputLayout(m_inputLayout.Get());
        GetContext()->IASetPrimitiveTopology(D3D11_PRIMITIVE_TOPOLOGY_TRIANGLELIST);
        GetContext()->VSSetShader(m_vs.Get(), nullptr, 0);
        GetContext()->PSSetShader(m_ps.Get(), nullptr, 0);

        VSConstants c{};
        c.viewProj = XMMatrixTranspose(m_cam.ViewProj());
        GetContext()->UpdateSubresource(m_cbVP.Get(), 0, nullptr, &c, 0, 0);
        GetContext()->VSSetConstantBuffers(0, 1, m_cbVP.GetAddressOf());

        for(const auto& kv : m_chunks) {
            const GPUChunk& g = kv.second;
            UINT stride = sizeof(VoxelVertex);
            UINT offset = 0;
            ID3D11Buffer* vb = g.vb.Get();
            GetContext()->IASetVertexBuffers(0, 1, &vb, &stride, &offset);
            GetContext()->IASetIndexBuffer(g.ib.Get(), DXGI_FORMAT_R32_UINT, 0);
            GetContext()->DrawIndexed(g.indexCount, 0, 0);
        }

    // Overlay: FPS at top-right
    DrawFPSOverlay();

        GetSwapChain()->Present(1, 0);
    }

    void OnResize(int w, int h) override {
        D3DApp::OnResize(w,h);
        m_cam.SetLens(XM_PIDIV4, float(w)/float(h), 0.1f, 1000.0f);
    }

private:
    struct GPUChunk { ComPtr<ID3D11Buffer> vb; ComPtr<ID3D11Buffer> ib; UINT indexCount=0; };

    bool CreatePipeline() {
        // Load/compile shaders
        UINT flags = D3DCOMPILE_ENABLE_STRICTNESS;
#ifdef _DEBUG
        flags |= D3DCOMPILE_DEBUG;
#endif
        ComPtr<ID3DBlob> vsb, psb, errors;
        HRESULT hr = D3DCompileFromFile(L"shaders/voxel.hlsl", nullptr, nullptr, "VSMain", "vs_5_0", flags, 0, &vsb, &errors);
        if(FAILED(hr)) return false;
        hr = D3DCompileFromFile(L"shaders/voxel.hlsl", nullptr, nullptr, "PSMain", "ps_5_0", flags, 0, &psb, &errors);
        if(FAILED(hr)) return false;
        GetDevice()->CreateVertexShader(vsb->GetBufferPointer(), vsb->GetBufferSize(), nullptr, m_vs.GetAddressOf());
        GetDevice()->CreatePixelShader(psb->GetBufferPointer(), psb->GetBufferSize(), nullptr, m_ps.GetAddressOf());

        D3D11_INPUT_ELEMENT_DESC il[] = {
            {"POSITION",0,DXGI_FORMAT_R32G32B32_FLOAT,0,0,D3D11_INPUT_PER_VERTEX_DATA,0},
            {"NORMAL",0,DXGI_FORMAT_R32G32B32_FLOAT,0,12,D3D11_INPUT_PER_VERTEX_DATA,0},
            {"COLOR",0,DXGI_FORMAT_R8G8B8A8_UNORM,0,24,D3D11_INPUT_PER_VERTEX_DATA,0},
        };
        GetDevice()->CreateInputLayout(il, _countof(il), vsb->GetBufferPointer(), vsb->GetBufferSize(), m_inputLayout.GetAddressOf());
        return true;
    }

    // ---------- Overlay (2D text) ----------
    struct OverlayVertex { DirectX::XMFLOAT2 pos; uint32_t color; };

    bool CreateOverlayPipeline() {
        UINT flags = D3DCOMPILE_ENABLE_STRICTNESS;
#ifdef _DEBUG
        flags |= D3DCOMPILE_DEBUG;
#endif
        ComPtr<ID3DBlob> vsb, psb, errors;
        if(FAILED(D3DCompileFromFile(L"shaders/overlay.hlsl", nullptr, nullptr, "VSMain", "vs_5_0", flags, 0, &vsb, &errors))) return false;
        if(FAILED(D3DCompileFromFile(L"shaders/overlay.hlsl", nullptr, nullptr, "PSMain", "ps_5_0", flags, 0, &psb, &errors))) return false;
        GetDevice()->CreateVertexShader(vsb->GetBufferPointer(), vsb->GetBufferSize(), nullptr, m_vsOverlay.GetAddressOf());
        GetDevice()->CreatePixelShader(psb->GetBufferPointer(), psb->GetBufferSize(), nullptr, m_psOverlay.GetAddressOf());

        D3D11_INPUT_ELEMENT_DESC il[] = {
            {"POSITION",0,DXGI_FORMAT_R32G32_FLOAT,0,0,D3D11_INPUT_PER_VERTEX_DATA,0},
            {"COLOR",0,DXGI_FORMAT_R8G8B8A8_UNORM,0,8,D3D11_INPUT_PER_VERTEX_DATA,0},
        };
        GetDevice()->CreateInputLayout(il, _countof(il), vsb->GetBufferPointer(), vsb->GetBufferSize(), m_layoutOverlay.GetAddressOf());

        // Dynamic vertex buffer (enough for small text overlay)
        D3D11_BUFFER_DESC vbd{}; vbd.Usage = D3D11_USAGE_DYNAMIC; vbd.ByteWidth = 65536 * sizeof(OverlayVertex); vbd.BindFlags = D3D11_BIND_VERTEX_BUFFER; vbd.CPUAccessFlags = D3D11_CPU_ACCESS_WRITE;
        if(FAILED(GetDevice()->CreateBuffer(&vbd, nullptr, m_vbOverlay.GetAddressOf()))) return false;

        // Constant buffer for screen size
        D3D11_BUFFER_DESC cbd{}; cbd.BindFlags = D3D11_BIND_CONSTANT_BUFFER; cbd.ByteWidth = 16; cbd.Usage = D3D11_USAGE_DEFAULT;
        GetDevice()->CreateBuffer(&cbd, nullptr, m_cbOverlay.GetAddressOf());

        // Depth disabled state
        D3D11_DEPTH_STENCIL_DESC dsd{}; dsd.DepthEnable = FALSE; dsd.DepthWriteMask = D3D11_DEPTH_WRITE_MASK_ZERO; dsd.DepthFunc = D3D11_COMPARISON_ALWAYS;
        GetDevice()->CreateDepthStencilState(&dsd, m_dssDisabled.GetAddressOf());

        // Rasterizer: no cull
        D3D11_RASTERIZER_DESC rs{}; rs.CullMode = D3D11_CULL_NONE; rs.FillMode = D3D11_FILL_SOLID; rs.DepthClipEnable = TRUE;
        GetDevice()->CreateRasterizerState(&rs, m_rsOverlay.GetAddressOf());

        return true;
    }

    static bool GetGlyph5x7(wchar_t ch, uint8_t rows[7]) {
        switch(ch) {
            case L'0': { uint8_t r[7]={0x1E,0x11,0x13,0x15,0x19,0x11,0x1E}; memcpy(rows,r,7); return true; }
            case L'1': { uint8_t r[7]={0x04,0x0C,0x14,0x04,0x04,0x04,0x1F}; memcpy(rows,r,7); return true; }
            case L'2': { uint8_t r[7]={0x1E,0x11,0x01,0x06,0x18,0x10,0x1F}; memcpy(rows,r,7); return true; }
            case L'3': { uint8_t r[7]={0x1E,0x11,0x01,0x0E,0x01,0x11,0x1E}; memcpy(rows,r,7); return true; }
            case L'4': { uint8_t r[7]={0x02,0x06,0x0A,0x12,0x1F,0x02,0x02}; memcpy(rows,r,7); return true; }
            case L'5': { uint8_t r[7]={0x1F,0x10,0x1E,0x01,0x01,0x11,0x1E}; memcpy(rows,r,7); return true; }
            case L'6': { uint8_t r[7]={0x0E,0x10,0x1E,0x11,0x11,0x11,0x0E}; memcpy(rows,r,7); return true; }
            case L'7': { uint8_t r[7]={0x1F,0x01,0x02,0x04,0x08,0x08,0x08}; memcpy(rows,r,7); return true; }
            case L'8': { uint8_t r[7]={0x0E,0x11,0x11,0x0E,0x11,0x11,0x0E}; memcpy(rows,r,7); return true; }
            case L'9': { uint8_t r[7]={0x0E,0x11,0x11,0x0F,0x01,0x02,0x0C}; memcpy(rows,r,7); return true; }
            case L'.': { uint8_t r[7]={0x00,0x00,0x00,0x00,0x00,0x06,0x06}; memcpy(rows,r,7); return true; }
            case L' ': { uint8_t r[7]={0,0,0,0,0,0,0}; memcpy(rows,r,7); return true; }
            case L'F': { uint8_t r[7]={0x1F,0x10,0x1E,0x10,0x10,0x10,0x10}; memcpy(rows,r,7); return true; }
            case L'P': { uint8_t r[7]={0x1E,0x11,0x11,0x1E,0x10,0x10,0x10}; memcpy(rows,r,7); return true; }
            case L'S': { uint8_t r[7]={0x0F,0x10,0x10,0x0E,0x01,0x01,0x1E}; memcpy(rows,r,7); return true; }
            case L'x': { uint8_t r[7]={0x00,0x11,0x0A,0x04,0x0A,0x11,0x00}; memcpy(rows,r,7); return true; }
            case L'(': { uint8_t r[7]={0x02,0x04,0x08,0x08,0x08,0x04,0x02}; memcpy(rows,r,7); return true; }
            case L')': { uint8_t r[7]={0x08,0x04,0x02,0x02,0x02,0x04,0x08}; memcpy(rows,r,7); return true; }
        }
        return false;
    }

    void DrawGlyphStringTopRight(const std::wstring& text, uint32_t color, int scale, int margin) {
        int charW = 6*scale; // 5 pixels + 1 spacing
        int totalW = (int)text.size() * charW;
        float startX = (float)(m_width - margin - totalW);
        float startY = (float)margin;
        float pxSize = (float)scale;
        m_overlayVerts.clear();
        float x = startX;
        for(wchar_t ch : text) {
            uint8_t rows[7];
            if(!GetGlyph5x7(ch, rows)) { x += charW; continue; }
            for(int r=0;r<7;++r){
                uint8_t bits = rows[r];
                for(int c=0;c<5;++c){
                    if(bits & (1 << (4-c))) {
                        float rx = x + c*pxSize;
                        float ry = startY + r*pxSize;
                        OverlayVertex v0{{rx,         ry        }, color};
                        OverlayVertex v1{{rx+pxSize,  ry        }, color};
                        OverlayVertex v2{{rx+pxSize,  ry+pxSize }, color};
                        OverlayVertex v3{{rx,         ry+pxSize }, color};
                        m_overlayVerts.push_back(v0);
                        m_overlayVerts.push_back(v1);
                        m_overlayVerts.push_back(v2);
                        m_overlayVerts.push_back(v0);
                        m_overlayVerts.push_back(v2);
                        m_overlayVerts.push_back(v3);
                    }
                }
            }
            x += (float)charW;
        }

        // Update screen size constant
        struct { float w; float h; float pad[2]; } cb{ (float)m_width, (float)m_height, {0,0} };
        GetContext()->UpdateSubresource(m_cbOverlay.Get(), 0, nullptr, &cb, 0, 0);

        if(!m_overlayVerts.empty()) {
            D3D11_MAPPED_SUBRESOURCE mapped{};
            if(SUCCEEDED(GetContext()->Map(m_vbOverlay.Get(), 0, D3D11_MAP_WRITE_DISCARD, 0, &mapped))) {
                size_t bytes = m_overlayVerts.size() * sizeof(OverlayVertex);
                memcpy(mapped.pData, m_overlayVerts.data(), bytes);
                GetContext()->Unmap(m_vbOverlay.Get(), 0);

                UINT stride = sizeof(OverlayVertex); UINT offset = 0; ID3D11Buffer* vb = m_vbOverlay.Get();
                GetContext()->OMSetDepthStencilState(m_dssDisabled.Get(), 0);
                GetContext()->RSSetState(m_rsOverlay.Get());
                GetContext()->IASetInputLayout(m_layoutOverlay.Get());
                GetContext()->IASetPrimitiveTopology(D3D11_PRIMITIVE_TOPOLOGY_TRIANGLELIST);
                GetContext()->VSSetShader(m_vsOverlay.Get(), nullptr, 0);
                GetContext()->PSSetShader(m_psOverlay.Get(), nullptr, 0);
                GetContext()->VSSetConstantBuffers(1, 1, m_cbOverlay.GetAddressOf());
                GetContext()->IASetVertexBuffers(0, 1, &vb, &stride, &offset);
                GetContext()->Draw((UINT)m_overlayVerts.size(), 0);

                // Restore 3D states
                GetContext()->OMSetDepthStencilState(m_dss.Get(), 0);
                GetContext()->RSSetState(m_rsSolid.Get());
                GetContext()->IASetInputLayout(m_inputLayout.Get());
                GetContext()->VSSetShader(m_vs.Get(), nullptr, 0);
                GetContext()->PSSetShader(m_ps.Get(), nullptr, 0);
            }
        }
    }

    void DrawFPSOverlay() {
        float fpsNow = (m_lastDt > 1e-6f) ? (1.0f / m_lastDt) : 0.0f;
        // Build text: "123.4 FPS" and add (x3) when sprinting
        auto isDown = [](int vk){ return (GetAsyncKeyState(vk) & 0x8000) != 0; };
        bool sprint = m_keys[VK_SHIFT] || isDown(VK_LSHIFT) || isDown(VK_RSHIFT);
        wchar_t buf[64]; swprintf(buf, 64, L"%.1f FPS%s", fpsNow, sprint ? L" (x3)" : L"");
        DrawGlyphStringTopRight(buf, 0xFFFFFFFFu, 2, 10);
    }

    bool BuildChunkGPU(int cx, int cz, GPUChunk& out) {
        if(!IsWithinView(cx, cz)) return false;
        ChunkMesh mesh = VoxelMesher::BuildChunkMesh(*m_world, cx, cz);
        if(mesh.indices.empty()) return false;

        D3D11_BUFFER_DESC vbd{}; vbd.ByteWidth = (UINT)(mesh.vertices.size()*sizeof(VoxelVertex)); vbd.BindFlags = D3D11_BIND_VERTEX_BUFFER; vbd.Usage = D3D11_USAGE_DEFAULT;
        D3D11_SUBRESOURCE_DATA vinit{ mesh.vertices.data(), 0, 0 };
        if(FAILED(GetDevice()->CreateBuffer(&vbd, &vinit, out.vb.GetAddressOf()))) return false;

        D3D11_BUFFER_DESC ibd{}; ibd.ByteWidth = (UINT)(mesh.indices.size()*sizeof(uint32_t)); ibd.BindFlags = D3D11_BIND_INDEX_BUFFER; ibd.Usage = D3D11_USAGE_DEFAULT;
        D3D11_SUBRESOURCE_DATA iinit{ mesh.indices.data(), 0, 0 };
        if(FAILED(GetDevice()->CreateBuffer(&ibd, &iinit, out.ib.GetAddressOf()))) return false;

        out.indexCount = (UINT)mesh.indices.size();
        return true;
    }

    static long long ChunkKey(int cx, int cz) {
        return ( (long long)cx << 32 ) ^ (unsigned long long)(unsigned int)cz;
    }

    void UpdateCameraChunk() {
        auto p = m_cam.Position();
        m_camChunkX = (int)floorf(p.x / CHUNK_SIZE);
        m_camChunkZ = (int)floorf(p.z / CHUNK_SIZE);
    }

    bool IsWithinView(int cx, int cz) const {
        int dx = cx - m_camChunkX;
        int dz = cz - m_camChunkZ;
        return (std::max(std::abs(dx), std::abs(dz)) <= m_viewRadius);
    }

    void EnqueueVisibleChunks() {
        // Generate list and sort by distance (closest first)
        std::vector<std::pair<int,int>> list;
        list.reserve((2*m_viewRadius+1)*(2*m_viewRadius+1));
        for(int dz=-m_viewRadius; dz<=m_viewRadius; ++dz){
            for(int dx=-m_viewRadius; dx<=m_viewRadius; ++dx){
                int cx = m_camChunkX + dx;
                int cz = m_camChunkZ + dz;
                long long key = ChunkKey(cx,cz);
                if(m_chunks.find(key) != m_chunks.end()) continue;
                list.emplace_back(cx,cz);
            }
        }
        auto sq = [](int a){ return a*a; };
        std::sort(list.begin(), list.end(), [&](auto&a, auto&b){
            int dax = a.first - m_camChunkX, daz = a.second - m_camChunkZ;
            int dbx = b.first - m_camChunkX, dbz = b.second - m_camChunkZ;
            return sq(dax)+sq(daz) > sq(dbx)+sq(dbz); // farthest first so we pop_back() closest
        });
        // Append to pending keeping uniqueness simple (could dedupe more)
        for(auto& c : list) m_pending.push_back(c);
    }

    void CleanupFarChunks() {
        std::vector<long long> toRemove;
        for(const auto& kv : m_chunks) {
            // unpack
            int cx = (int)(kv.first >> 32);
            int cz = (int)(kv.first & 0xffffffff);
            if(!IsWithinView(cx, cz)) {
                toRemove.push_back(kv.first);
            }
        }
        for(auto key : toRemove) m_chunks.erase(key);
    }

    LRESULT MsgProc(HWND hwnd, UINT msg, WPARAM wParam, LPARAM lParam) override {
        switch(msg) {
            case WM_KEYDOWN: {
                int vk = (int)wParam;
                if(vk >= 0 && vk < 256) m_keys[vk] = true;
                if(vk == VK_ESCAPE && m_mouseCaptured) {
                    ReleaseCapture();
                    ShowCursor(TRUE);
                    m_mouseCaptured = false;
                }
                break;
            }
            case WM_KEYUP: {
                int vk = (int)wParam;
                if(vk >= 0 && vk < 256) m_keys[vk] = false;
                break;
            }
            case WM_RBUTTONDOWN: {
                SetCapture(hwnd);
                ShowCursor(FALSE);
                m_mouseCaptured = true;
                m_lastMouse.x = GET_X_LPARAM(lParam);
                m_lastMouse.y = GET_Y_LPARAM(lParam);
                break;
            }
            case WM_RBUTTONUP: {
                if(m_mouseCaptured) {
                    ReleaseCapture();
                    ShowCursor(TRUE);
                    m_mouseCaptured = false;
                }
                break;
            }
            case WM_MOUSEMOVE: {
                if(m_mouseCaptured) {
                    int mx = GET_X_LPARAM(lParam);
                    int my = GET_Y_LPARAM(lParam);
                    int dx = mx - m_lastMouse.x;
                    int dy = my - m_lastMouse.y;
                    m_lastMouse.x = mx;
                    m_lastMouse.y = my;

                    float yaw = dx * m_mouseSens;
                    float pitch = -dy * m_mouseSens; // invert Y for natural feel

                    // clamp pitch to avoid flipping
                    float newPitch = m_pitchAngle + pitch;
                    const float limit = 1.2f; // ~69 deg
                    float applyPitch = pitch;
                    if(newPitch > limit) applyPitch = limit - m_pitchAngle;
                    else if(newPitch < -limit) applyPitch = -limit - m_pitchAngle;
                    m_pitchAngle += applyPitch;

                    m_cam.RotateY(yaw);
                    m_cam.Pitch(applyPitch);
                }
                break;
            }
        }
        return D3DApp::MsgProc(hwnd, msg, wParam, lParam);
    }

    Camera m_cam;
    std::unique_ptr<VoxelWorld> m_world;
    std::array<bool,256> m_keys{};
    bool m_mouseCaptured = false;
    POINT m_lastMouse{0,0};
    float m_moveSpeed = 10.0f; // units per second
    float m_mouseSens = 0.0025f; // radians per pixel
    float m_pitchAngle = 0.0f;
    // FPS
    float m_fpsAccum = 0.0f;
    int m_fpsFrames = 0;
    float m_lastDt = 0.0f;

    // Overlay pipeline
    ComPtr<ID3D11VertexShader> m_vsOverlay;
    ComPtr<ID3D11PixelShader> m_psOverlay;
    ComPtr<ID3D11InputLayout> m_layoutOverlay;
    ComPtr<ID3D11Buffer> m_vbOverlay;
    ComPtr<ID3D11Buffer> m_cbOverlay;
    ComPtr<ID3D11DepthStencilState> m_dssDisabled;
    ComPtr<ID3D11RasterizerState> m_rsOverlay;
    std::vector<OverlayVertex> m_overlayVerts;

    std::unordered_map<long long, GPUChunk> m_chunks;
    std::vector<std::pair<int,int>> m_pending;
    int m_viewRadius = 32;
    int m_camChunkX = 0, m_camChunkZ = 0;

    ComPtr<ID3D11VertexShader> m_vs;
    ComPtr<ID3D11PixelShader> m_ps;
    ComPtr<ID3D11InputLayout> m_inputLayout;
    ComPtr<ID3D11Buffer> m_cbVP;
    ComPtr<ID3D11RasterizerState> m_rsSolid;
    ComPtr<ID3D11DepthStencilState> m_dss;

    float m_angle = 0.0f;
};

int WINAPI wWinMain(HINSTANCE hInst, HINSTANCE, PWSTR, int) {
    VoxelApp app(hInst);
    if(!app.Initialize()) return -1;
    return app.Run();
}
