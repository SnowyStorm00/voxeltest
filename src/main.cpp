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
#include <thread>
#include <mutex>
#include <condition_variable>
#include <chrono>
#include <limits>
#include "obj_loader.h"

using namespace DirectX;
using Microsoft::WRL::ComPtr;

struct VSConstants {
    XMMATRIX viewProj;
    XMFLOAT2 jitter; // unused (TAA removed); kept for cbuffer layout stability with voxel.hlsl
    XMFLOAT2 pad;
    XMFLOAT3 cameraPos;
    float fogStart;
    float fogEnd;
    XMFLOAT3 fogColor;
    XMFLOAT3 sunDir;
    float    pad2;
};

class VoxelApp : public D3DApp {
public:
    VoxelApp(HINSTANCE hInst) : D3DApp(hInst, 1280, 720, L"Voxel RPG") {}
    ~VoxelApp() {
        // Stop workers and join
        {
            std::lock_guard<std::mutex> lk(m_jobMutex);
            m_stopWorkers = true;
        }
        m_jobCv.notify_all();
        for(auto &t : m_workers){ if(t.joinable()) t.join(); }
    }

public:
    bool Initialize() override {
        if(!D3DApp::Initialize()) return false;
    // High-precision timer setup for FPS capping
    LARGE_INTEGER f; QueryPerformanceFrequency(&f); m_qpcFreq = (long long)f.QuadPart;

    // Camera
    m_cam.SetLens(XM_PIDIV4, float(m_width)/float(m_height), 0.1f, 1000.0f);
    m_cam.LookAt(XMVectorSet(8,12,-20,1), XMVectorSet(8,8,8,1), XMVectorSet(0,1,0,0));
    m_cam.UpdateViewMatrix();

    // Defer world creation until user enters a seed
    m_needSeedInput = true;
    m_seedInput.clear();

    // Create shaders and input layout
    if(!CreatePipeline()) return false;
    // TAA removed

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

    // Sky pipeline
    if(!CreateSkyPipeline()) return false;

    // Overlay pipeline (2D quads)
    if(!CreateOverlayPipeline()) return false;

    // TAA removed: no post-process targets

    // Create MSAA targets if enabled
    if(m_msaaEnabled) CreateMsaaTargets();

    // Shadows via software ray test (no shadow map)

    // Load player model
    LoadPlayerModel();

    // Initial title so users see something immediately
    SetWindowTextW(GetHwnd(), L"Voxel RPG - enter seed");
    // Do not capture mouse while typing seed; capture after loading begins
    if(!m_showMenu && !m_needSeedInput) CaptureMouse();
    return true;
    }

    void Update(float dt) override {
    // Mark frame start time for accurate FPS capping (includes Update + Render)
    { LARGE_INTEGER t; QueryPerformanceCounter(&t); m_frameStartQpc = (long long)t.QuadPart; }
    // Clamp dt to reduce sensitivity to unstable frame times (e.g., vsync off)
    if(dt > 0.05f) dt = 0.05f; // cap ~20 FPS step max
    m_lastDt = dt;
        // Seed entry screen: block gameplay until seed chosen
        if(m_needSeedInput){
            float fpsNow = (m_lastDt > 1e-6f) ? (1.0f / m_lastDt) : 0.0f;
            std::wstringstream ss; ss.setf(std::ios::fixed); ss.precision(1);
            ss << L"Voxel RPG - enter seed  (" << fpsNow << L" FPS)";
            SetWindowTextW(GetHwnd(), ss.str().c_str());
            return;
        }
        // During initial loading, keep uploading meshes and skip gameplay update
        if(m_inInitialLoading){
            // Upload as many as possible to drive progress
            const int MaxUploadsPerFrame = 512;
            UploadReadyMeshesPerFrame(MaxUploadsPerFrame);
            // Finished?
            if(m_initialTotal > 0 && m_initialLoaded >= m_initialTotal){
                m_inInitialLoading = false;
                // Capture mouse to start gameplay
                if(!m_showMenu) CaptureMouse();
                // Now schedule HQ remeshes of all resident chunks in view
                std::vector<std::pair<int,int>> relist;
                for(int dz=-m_viewRadius; dz<=m_viewRadius; ++dz){
                    for(int dx=-m_viewRadius; dx<=m_viewRadius; ++dx){
                        relist.emplace_back(m_camChunkX+dx, m_camChunkZ+dz);
                    }
                }
                // Clear requested flags so they can be queued again
                m_requested.clear();
                QueueMeshingJobs(relist);
            }
            // Still loading: minimal work, keep title updated
            float fpsNow = (m_lastDt > 1e-6f) ? (1.0f / m_lastDt) : 0.0f;
            std::wstringstream ss; ss.setf(std::ios::fixed); ss.precision(1);
            int pct = (m_initialTotal>0) ? (int)std::round(100.0 * (double)m_initialLoaded / (double)m_initialTotal) : 0;
            ss << L"Voxel RPG - Loading " << pct << L"%  (" << fpsNow << L" FPS)";
            SetWindowTextW(GetHwnd(), ss.str().c_str());
            return;
        }

        // Freecam: bypass player physics; move camera directly when enabled
        if(m_freecam){
            auto isDown = [](int vk){ return (GetAsyncKeyState(vk) & 0x8000) != 0; };
            bool sprint = m_keys[VK_SHIFT] || isDown(VK_LSHIFT) || isDown(VK_RSHIFT);
            float speed = m_freecamSpeed * (sprint ? 3.0f : 1.0f);
            if(!m_showMenu){
                XMFLOAT3 look = m_cam.Look();
                XMFLOAT3 right = m_cam.Right();
                // Horizontal movement aligned to yaw; vertical via space/ctrl
                XMVECTOR f = XMVector3Normalize(XMVectorSet(look.x, 0.0f, look.z, 0.0f));
                XMVECTOR r = XMVector3Normalize(XMVectorSet(right.x, 0.0f, right.z, 0.0f));
                XMVECTOR dir = XMVectorZero();
                if(m_keys['W'] || isDown('W')) dir = XMVectorAdd(dir, f);
                if(m_keys['S'] || isDown('S')) dir = XMVectorSubtract(dir, f);
                if(m_keys['A'] || isDown('A')) dir = XMVectorSubtract(dir, r);
                if(m_keys['D'] || isDown('D')) dir = XMVectorAdd(dir, r);
                float len = XMVectorGetX(XMVector3Length(dir));
                if(len > 1e-4f) dir = XMVectorScale(dir, 1.0f / len);
                float vY = 0.0f;
                if(m_keys[VK_SPACE] || isDown(VK_SPACE)) vY += 1.0f;
                if(m_keys[VK_LCONTROL] || isDown(VK_LCONTROL) || m_keys[VK_RCONTROL] || isDown(VK_RCONTROL)) vY -= 1.0f;
                XMFLOAT3 cpos = m_cam.Position();
                XMVECTOR p = XMVectorSet(cpos.x, cpos.y, cpos.z, 1.0f);
                XMVECTOR delta = XMVectorAdd(XMVectorScale(dir, speed * dt), XMVectorSet(0, vY * speed * dt, 0, 0));
                p = XMVectorAdd(p, delta);
                XMFLOAT3 np; XMStoreFloat3(&np, p);
                m_cam.SetPosition(np.x, np.y, np.z);
                m_cam.UpdateViewMatrix();
            }
            // Hide player in freecam
            m_showPlayer = false;
        } else {
        // Check grounded against world blocks to make jump robust across frame times
        auto IsOnGroundNow = [&](){
            const float eps = 0.05f;
            float y = m_playerPos.y - m_playerHalf.y - eps;
            int yb = (int)std::floor(y);
            float minX = m_playerPos.x - m_playerHalf.x + 0.001f;
            float maxX = m_playerPos.x + m_playerHalf.x - 0.001f;
            float minZ = m_playerPos.z - m_playerHalf.z + 0.001f;
            float maxZ = m_playerPos.z + m_playerHalf.z - 0.001f;
            int x0 = (int)std::floor(minX), x1 = (int)std::floor(maxX);
            int z0 = (int)std::floor(minZ), z1 = (int)std::floor(maxZ);
            for(int z=z0; z<=z1; ++z){
                for(int x=x0; x<=x1; ++x){
                    if(IsSolid(x, yb, z)) return true;
                }
            }
            return false;
        };
        bool groundedNow = IsOnGroundNow();
        m_onGround = groundedNow;
        // Player movement with gravity + collision; camera anchored to player head
        auto isDown = [](int vk){ return (GetAsyncKeyState(vk) & 0x8000) != 0; };
        bool sprint = m_keys[VK_SHIFT] || isDown(VK_LSHIFT) || isDown(VK_RSHIFT);
        float speed = m_moveSpeed * (sprint ? 3.0f : 1.0f);
        if(!m_showMenu){
            // Determine desired horizontal velocity from camera facing (yaw only)
            XMFLOAT3 look = m_cam.Look();
            XMFLOAT3 right = m_cam.Right();
            XMVECTOR f = XMVector3Normalize(XMVectorSet(look.x, 0.0f, look.z, 0.0f));
            XMVECTOR r = XMVector3Normalize(XMVectorSet(right.x, 0.0f, right.z, 0.0f));
            XMVECTOR dir = XMVectorZero();
            if(m_keys['W'] || isDown('W')) dir = XMVectorAdd(dir, f);
            if(m_keys['S'] || isDown('S')) dir = XMVectorSubtract(dir, f);
            if(m_keys['A'] || isDown('A')) dir = XMVectorSubtract(dir, r);
            if(m_keys['D'] || isDown('D')) dir = XMVectorAdd(dir, r);
            float len = XMVectorGetX(XMVector3Length(dir));
            XMFLOAT3 hv{0,0,0};
            if(len > 1e-4f){
                XMVECTOR n = XMVectorScale(dir, speed / len);
                XMStoreFloat3(&hv, n);
            }
            // Apply horizontal desired velocity directly; keep vertical velocity from gravity
            m_playerVel.x = hv.x;
            m_playerVel.z = hv.z;
            // Jump (space) when on ground; edge triggered
            bool space = m_keys[VK_SPACE] || isDown(VK_SPACE);
            if(space && !m_jumpHeld && groundedNow){
                m_playerVel.y = 8.5f;
            }
            m_jumpHeld = space;
        }
        // Gravity
        m_playerVel.y += m_gravity * dt;
        if(m_playerVel.y < -80.0f) m_playerVel.y = -80.0f;
        // Integrate with collision
        CollideAndSlide(dt);
        // Third-person camera: position behind player horizontally, with slight over-the-shoulder offset
        XMFLOAT3 lookDir = m_cam.Look();
        XMFLOAT3 rightDir = m_cam.Right();
        // Horizontal-only follow to keep distance consistent
        XMVECTOR fwdH = XMVector3Normalize(XMVectorSet(lookDir.x, 0.0f, lookDir.z, 0.0f));
        XMVECTOR rightH = XMVector3Normalize(XMVectorSet(rightDir.x, 0.0f, rightDir.z, 0.0f));
        float camDist = m_thirdPersonDist;
        float shoulder = 0.6f; // over-the-shoulder
        float headRaise = 0.2f; // slightly above head
        XMVECTOR target = XMVectorSet(m_playerPos.x, m_playerPos.y + m_eyeHeight + headRaise, m_playerPos.z, 1.0f);
        // If target is inside a solid (tight corridors), nudge up to find a free spot
        {
            XMFLOAT3 tf; XMStoreFloat3(&tf, target);
            int bx = (int)std::floor(tf.x);
            int by = (int)std::floor(tf.y);
            int bz = (int)std::floor(tf.z);
            int safety = 6;
            while(IsSolid(bx, by, bz) && safety-- > 0){
                tf.y += 0.25f; by = (int)std::floor(tf.y);
            }
            target = XMVectorSet(tf.x, tf.y, tf.z, 1.0f);
        }
        XMVECTOR desired = XMVectorAdd(
            XMVectorSubtract(target, XMVectorScale(fwdH, camDist)),
            XMVectorScale(rightH, shoulder)
        );
        // Obstruction clamp: march from target outwards to desired, keep the farthest free point
        XMFLOAT3 desiredF; XMStoreFloat3(&desiredF, desired);
        XMFLOAT3 targetF; XMStoreFloat3(&targetF, target);
    XMFLOAT3 camPosF = desiredF; // default to far if unobstructed
        const float step = 0.2f;
        float total = std::max(0.001f, camDist);
    XMFLOAT3 lastFree = targetF;
        for(float s = 0.0f; s <= total; s += step){
            float u = s / total; // 0 at target, 1 at desired
            XMVECTOR p = XMVectorLerp(target, desired, u);
            XMFLOAT3 pf; XMStoreFloat3(&pf, p);
            int bx = (int)std::floor(pf.x);
            int by = (int)std::floor(pf.y);
            int bz = (int)std::floor(pf.z);
            if(IsSolid(bx, by, bz)){
                // hit obstruction; stop and use last free point
                // move a tiny bit toward desired to avoid exact overlap with target
                XMVECTOR last = XMLoadFloat3(&lastFree);
                XMVECTOR tinyForward = XMVectorScale(fwdH, 0.15f);
                XMVECTOR safe = XMVectorAdd(last, tinyForward);
                XMFLOAT3 safeF; XMStoreFloat3(&safeF, safe);
                camPosF = safeF;
                goto set_cam;
            }
            lastFree = pf;
        }
        camPosF = lastFree; // reached desired with no obstruction
set_cam:
        m_cam.SetPosition(camPosF.x, camPosF.y, camPosF.z);
        if(!m_showMenu) m_cam.UpdateViewMatrix();
    // Ensure player is visible in third-person mode
    m_showPlayer = true;
    }

        // Update streaming center when crossing chunk boundaries
        int prevCX = m_camChunkX, prevCZ = m_camChunkZ;
        UpdateCameraChunk();
        if(prevCX != m_camChunkX || prevCZ != m_camChunkZ) {
            EnqueueVisibleChunks();
            CleanupFarChunks();
        }

    // Upload a batch of ready meshes from worker threads
    const int MaxUploadsPerFrame = 128;
    UploadReadyMeshesPerFrame(MaxUploadsPerFrame);

    // FPS counter: instantaneous (current frame)
    m_fpsAccum += dt;
    m_fpsFrames += 1;
    float fpsNow = (m_lastDt > 1e-6f) ? (1.0f / m_lastDt) : 0.0f;
        {
            std::wstringstream ss;
            ss.setf(std::ios::fixed); ss.precision(1);
            ss << L"Voxel RPG - " << fpsNow << L" FPS";
            auto isDown2 = [](int vk){ return (GetAsyncKeyState(vk) & 0x8000) != 0; };
            bool sprintNow = m_keys[VK_SHIFT] || isDown2(VK_LSHIFT) || isDown2(VK_RSHIFT);
            if(sprintNow) ss << L"  (x3)";
            SetWindowTextW(GetHwnd(), ss.str().c_str());
        }
    }

    void Render() override {
    float clear[4] = {0.53f, 0.81f, 0.92f, 1.0f};
    ID3D11RenderTargetView* mainRTV = GetRTV();
    ID3D11DepthStencilView* mainDSV = GetDSV();
    if(!mainRTV || !mainDSV){
        return; // likely minimized or resizing
    }

    // Proactively clear PS SRVs to avoid SRV/RTV hazards at frame start
    { ID3D11ShaderResourceView* nullSrvs[8] = {nullptr,nullptr,nullptr,nullptr,nullptr,nullptr,nullptr,nullptr}; GetContext()->PSSetShaderResources(0, 8, nullSrvs); }
    float inner = float((m_viewRadius - 2) * CHUNK_SIZE);
    float outer = float((m_viewRadius + 1) * CHUNK_SIZE);
    if(inner < 0) inner = 0.0f;
    if(outer < inner + 1.0f) outer = inner + 1.0f;

    // Choose render target: MSAA if enabled, else backbuffer
    ID3D11RenderTargetView* rtv3D = (m_msaaEnabled && m_msaaRTV ? m_msaaRTV.Get() : mainRTV);
    ID3D11DepthStencilView* dsv3D = (m_msaaEnabled && m_msaaDSV ? m_msaaDSV.Get() : mainDSV);
    if(!rtv3D || !dsv3D) return; // invalid state, skip frame
    GetContext()->OMSetRenderTargets(1, &rtv3D, dsv3D);
    // Restore screen viewport
    {
        D3D11_VIEWPORT vp{0}; vp.TopLeftX = 0; vp.TopLeftY = 0; vp.Width = (float)m_width; vp.Height = (float)m_height; vp.MinDepth = 0.0f; vp.MaxDepth = 1.0f;
        GetContext()->RSSetViewports(1, &vp);
    }
    // Depth clear only; color will be covered by sky
    GetContext()->ClearDepthStencilView(dsv3D, D3D11_CLEAR_DEPTH|D3D11_CLEAR_STENCIL, 1.0f, 0);

        GetContext()->OMSetDepthStencilState(m_dss.Get(), 0);
        GetContext()->RSSetState(m_rsSolid.Get());
        GetContext()->IASetInputLayout(m_inputLayout.Get());
        GetContext()->IASetPrimitiveTopology(D3D11_PRIMITIVE_TOPOLOGY_TRIANGLELIST);
        GetContext()->VSSetShader(m_vs.Get(), nullptr, 0);
        GetContext()->PSSetShader(m_ps.Get(), nullptr, 0);

    // Jitter unused (TAA removed). Keep zero.
    VSConstants c{};
    c.viewProj = XMMatrixTranspose(m_cam.ViewProj());
    c.jitter = XMFLOAT2{0,0};
    c.cameraPos = m_cam.Position();
    // Fog: start near the inner edge and end at max view distance (in world units)
    c.fogStart = inner;
    c.fogEnd = outer;
    c.fogColor = XMFLOAT3(0.53f, 0.81f, 0.92f); // sky color
    c.sunDir = m_sunDir;
    GetContext()->UpdateSubresource(m_cbVP.Get(), 0, nullptr, &c, 0, 0);
    GetContext()->VSSetConstantBuffers(0, 1, m_cbVP.GetAddressOf());
    // PS also declares the same cbuffer in HLSL; bind it to avoid undefined state on some drivers
    GetContext()->PSSetConstantBuffers(0, 1, m_cbVP.GetAddressOf());
    // No shadow SRV bound

    // Draw sky background (uses its own cbuffer at slot b0). Ensure no shadow SRV bound.
    DrawSky();
    // Show seed input screen before world exists
    if(m_needSeedInput){
        DrawSeedInputScreen();
        HRESULT __presentHr = GetSwapChain()->Present(m_vsyncEnabled ? 1 : 0, 0);
        if(__presentHr == DXGI_ERROR_DEVICE_REMOVED || __presentHr == DXGI_ERROR_DEVICE_RESET){ PostMessage(GetHwnd(), WM_CLOSE, 0, 0); }
        return;
    }
    // If in loading phase, draw loading screen and skip world rendering
    if(m_inInitialLoading){
        DrawLoadingScreen();
        // Present with vsync according to setting
        HRESULT __presentHr = GetSwapChain()->Present(m_vsyncEnabled ? 1 : 0, 0);
        if(__presentHr == DXGI_ERROR_DEVICE_REMOVED || __presentHr == DXGI_ERROR_DEVICE_RESET){ PostMessage(GetHwnd(), WM_CLOSE, 0, 0); }
        return;
    }
    // Rebind voxel constants for subsequent world rendering
    GetContext()->VSSetConstantBuffers(0, 1, m_cbVP.GetAddressOf());
    GetContext()->PSSetConstantBuffers(0, 1, m_cbVP.GetAddressOf());
    // No shadow SRV bound

        // Build frustum planes from current VP for culling
        UpdateFrustumPlanes();

    for(const auto& kv : m_chunks) {
            // Frustum cull per chunk AABB
            int cx = (int)(kv.first >> 32);
            int cz = (int)(kv.first & 0xffffffff);
            if(!ChunkInFrustum(cx, cz)) continue;

            const GPUChunk& g = kv.second;
            UINT stride = sizeof(VoxelVertex);
            UINT offset = 0;
            ID3D11Buffer* vb = g.vb.Get();
            GetContext()->IASetVertexBuffers(0, 1, &vb, &stride, &offset);
            GetContext()->IASetIndexBuffer(g.ib.Get(), DXGI_FORMAT_R32_UINT, 0);
            GetContext()->DrawIndexed(g.indexCount, 0, 0);
        }

    // Draw player (3D) before MSAA resolve so it uses the same depth buffer
    if(!m_showMenu) DrawPlayer();

    // No shadow SRV bound

    // If we rendered with MSAA, resolve MSAA color to backbuffer before overlays
    if(m_msaaEnabled && m_msaaTex && GetBackBufferTexture()){
            GetContext()->ResolveSubresource(GetBackBufferTexture(), 0, m_msaaTex.Get(), 0, DXGI_FORMAT_R8G8B8A8_UNORM);
            // Bind backbuffer for overlays
            ID3D11RenderTargetView* rtvBB = mainRTV;
            GetContext()->OMSetRenderTargets(1, &rtvBB, mainDSV);
        }

    // Overlay: FPS at top-right
    DrawFPSOverlay();
        if(m_showMenu) DrawSettingsMenu();

        // Present with optional vsync
        HRESULT __presentHr = GetSwapChain()->Present(m_vsyncEnabled ? 1 : 0, 0);
        if(__presentHr == DXGI_ERROR_DEVICE_REMOVED || __presentHr == DXGI_ERROR_DEVICE_RESET){
            PostMessage(GetHwnd(), WM_CLOSE, 0, 0);
            return;
        }

        // Optional FPS cap using high-resolution timing (0 = uncapped)
        if(m_fpsCap > 0 && m_qpcFreq > 0){
            long long targetTicks = m_qpcFreq / (long long)m_fpsCap;
            LARGE_INTEGER nowLI; QueryPerformanceCounter(&nowLI);
            long long now = (long long)nowLI.QuadPart;
            long long endTicks = m_frameStartQpc + targetTicks;
            long long remain = endTicks - now;
            if(remain > 0){
                // Coarse sleep for most of the remaining time
                double remainMs = 1000.0 * (double)remain / (double)m_qpcFreq;
                if(remainMs > 2.0){
                    DWORD ms = (DWORD)(remainMs - 1.0);
                    if(ms > 0) Sleep(ms);
                }
                // Fine-grained yield until target reached
                for(;;){
                    QueryPerformanceCounter(&nowLI);
                    if((long long)nowLI.QuadPart >= endTicks) break;
                    Sleep(0);
                }
            }
        }
    }

    void OnResize(int w, int h) override {
        D3DApp::OnResize(w,h);
        m_cam.SetLens(XM_PIDIV4, float(w)/float(h), 0.1f, 1000.0f);
    // Recreate MSAA targets
    if(m_msaaEnabled) CreateMsaaTargets();
    // Update mouse clip region on resize
    if(m_mouseCaptured && !m_showMenu) UpdateMouseLockRect();
    }

private:
    struct GPUChunk { ComPtr<ID3D11Buffer> vb; ComPtr<ID3D11Buffer> ib; UINT indexCount=0; };
    struct CPUChunkMesh { int cx, cz; ChunkMesh mesh; };
    struct CachedMesh { ChunkMesh mesh; unsigned version=0; };

    bool CreatePipeline() {
        // Load/compile shaders
        UINT flags = D3DCOMPILE_ENABLE_STRICTNESS;
#ifdef _DEBUG
        flags |= D3DCOMPILE_DEBUG;
#endif
        ComPtr<ID3DBlob> vsb, psb, errors;
        HRESULT hr = D3DCompileFromFile(L"shaders/voxel.hlsl", nullptr, nullptr, "VSMain", "vs_5_0", flags, 0, &vsb, &errors);
        if(FAILED(hr)){
            if(errors) MessageBoxA(GetHwnd(), (const char*)errors->GetBufferPointer(), "VS Compile Error", MB_OK|MB_ICONERROR);
            return false;
        }
        hr = D3DCompileFromFile(L"shaders/voxel.hlsl", nullptr, nullptr, "PSMain", "ps_5_0", flags, 0, &psb, &errors);
        if(FAILED(hr)){
            if(errors) MessageBoxA(GetHwnd(), (const char*)errors->GetBufferPointer(), "PS Compile Error", MB_OK|MB_ICONERROR);
            return false;
        }
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

    // TAA removed: no fullscreen quad pipeline

    // ---------- Overlay (2D text) ----------
    struct OverlayVertex { DirectX::XMFLOAT2 pos; uint32_t color; };

    bool CreateOverlayPipeline() {
        UINT flags = D3DCOMPILE_ENABLE_STRICTNESS;
#ifdef _DEBUG
        flags |= D3DCOMPILE_DEBUG;
#endif
        ComPtr<ID3DBlob> vsb, psb, errors;
        if(FAILED(D3DCompileFromFile(L"shaders/overlay.hlsl", nullptr, nullptr, "VSMain", "vs_5_0", flags, 0, &vsb, &errors))){
            if(errors) MessageBoxA(GetHwnd(), (const char*)errors->GetBufferPointer(), "Overlay VS Compile Error", MB_OK|MB_ICONERROR);
            return false;
        }
        if(FAILED(D3DCompileFromFile(L"shaders/overlay.hlsl", nullptr, nullptr, "PSMain", "ps_5_0", flags, 0, &psb, &errors))){
            if(errors) MessageBoxA(GetHwnd(), (const char*)errors->GetBufferPointer(), "Overlay PS Compile Error", MB_OK|MB_ICONERROR);
            return false;
        }
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

    // ---------- Sky ----------
    struct SkyVertex { XMFLOAT2 pos; };
    struct SkyCB { XMMATRIX invViewProj; XMFLOAT3 sunDir; float sunSize; XMFLOAT3 skyHorizon; float pad0; XMFLOAT3 skyZenith; float pad1; };

    bool CreateSkyPipeline(){
        UINT flags = D3DCOMPILE_ENABLE_STRICTNESS;
#ifdef _DEBUG
        flags |= D3DCOMPILE_DEBUG;
#endif
        ComPtr<ID3DBlob> vsb, psb, errors;
        if(FAILED(D3DCompileFromFile(L"shaders/sky.hlsl", nullptr, nullptr, "VSMain", "vs_5_0", flags, 0, &vsb, &errors))){
            if(errors) MessageBoxA(GetHwnd(), (const char*)errors->GetBufferPointer(), "Sky VS Compile Error", MB_OK|MB_ICONERROR);
            return false;
        }
        if(FAILED(D3DCompileFromFile(L"shaders/sky.hlsl", nullptr, nullptr, "PSMain", "ps_5_0", flags, 0, &psb, &errors))){
            if(errors) MessageBoxA(GetHwnd(), (const char*)errors->GetBufferPointer(), "Sky PS Compile Error", MB_OK|MB_ICONERROR);
            return false;
        }
        GetDevice()->CreateVertexShader(vsb->GetBufferPointer(), vsb->GetBufferSize(), nullptr, m_vsSky.GetAddressOf());
        GetDevice()->CreatePixelShader(psb->GetBufferPointer(), psb->GetBufferSize(), nullptr, m_psSky.GetAddressOf());
        D3D11_INPUT_ELEMENT_DESC il[] = { {"POSITION",0,DXGI_FORMAT_R32G32_FLOAT,0,0,D3D11_INPUT_PER_VERTEX_DATA,0} };
        GetDevice()->CreateInputLayout(il, _countof(il), vsb->GetBufferPointer(), vsb->GetBufferSize(), m_layoutSky.GetAddressOf());
        // Fullscreen triangle
        SkyVertex tri[3] = { {{-1.0f,-1.0f}}, {{3.0f,-1.0f}}, {{-1.0f,3.0f}} };
        D3D11_BUFFER_DESC vbd{}; vbd.Usage=D3D11_USAGE_IMMUTABLE; vbd.BindFlags=D3D11_BIND_VERTEX_BUFFER; vbd.ByteWidth=sizeof(tri);
        D3D11_SUBRESOURCE_DATA vinit{ tri, 0, 0 };
        if(FAILED(GetDevice()->CreateBuffer(&vbd, &vinit, m_vbSky.GetAddressOf()))) return false;
        // Constant buffer
        D3D11_BUFFER_DESC cbd{}; cbd.BindFlags=D3D11_BIND_CONSTANT_BUFFER; cbd.ByteWidth=sizeof(SkyCB); cbd.Usage=D3D11_USAGE_DEFAULT;
        GetDevice()->CreateBuffer(&cbd, nullptr, m_cbSky.GetAddressOf());
        return true;
    }

    void DrawSky(){
        // Depth disabled, draw background
        GetContext()->OMSetDepthStencilState(m_dssDisabled.Get(), 0);
        GetContext()->RSSetState(m_rsSolid.Get());
        GetContext()->IASetInputLayout(m_layoutSky.Get());
        GetContext()->IASetPrimitiveTopology(D3D11_PRIMITIVE_TOPOLOGY_TRIANGLELIST);
        UINT stride = sizeof(SkyVertex), offset = 0; ID3D11Buffer* vb = m_vbSky.Get();
        GetContext()->IASetVertexBuffers(0,1,&vb,&stride,&offset);
        GetContext()->VSSetShader(m_vsSky.Get(), nullptr, 0);
        GetContext()->PSSetShader(m_psSky.Get(), nullptr, 0);
        // Build constants
        XMMATRIX invVP = XMMatrixInverse(nullptr, m_cam.ViewProj());
        XMMATRIX invVPT = XMMatrixTranspose(invVP);
        SkyCB cb{}; cb.invViewProj = invVPT; cb.sunDir = m_sunDir; cb.sunSize = m_sunSize;
        cb.skyHorizon = m_skyHorizon; cb.skyZenith = m_skyZenith;
        GetContext()->UpdateSubresource(m_cbSky.Get(), 0, nullptr, &cb, 0, 0);
        GetContext()->VSSetConstantBuffers(0,1,m_cbSky.GetAddressOf());
        GetContext()->PSSetConstantBuffers(0,1,m_cbSky.GetAddressOf());
        GetContext()->Draw(3,0);
        // Restore depth for 3D
        GetContext()->OMSetDepthStencilState(m_dss.Get(), 0);
        GetContext()->IASetInputLayout(m_inputLayout.Get());
        GetContext()->VSSetShader(m_vs.Get(), nullptr, 0);
        GetContext()->PSSetShader(m_ps.Get(), nullptr, 0);
    }

    static bool GetGlyph5x7(wchar_t ch, uint8_t rows[7]) {
        switch(ch) {
            // Uppercase letters (5x7)
            case L'A': { uint8_t r[7]={0x0E,0x11,0x11,0x1F,0x11,0x11,0x11}; memcpy(rows,r,7); return true; }
            case L'B': { uint8_t r[7]={0x1E,0x11,0x11,0x1E,0x11,0x11,0x1E}; memcpy(rows,r,7); return true; }
            case L'C': { uint8_t r[7]={0x0E,0x11,0x10,0x10,0x10,0x11,0x0E}; memcpy(rows,r,7); return true; }
            case L'D': { uint8_t r[7]={0x1E,0x11,0x11,0x11,0x11,0x11,0x1E}; memcpy(rows,r,7); return true; }
            case L'E': { uint8_t r[7]={0x1F,0x10,0x10,0x1E,0x10,0x10,0x1F}; memcpy(rows,r,7); return true; }
            case L'F': { uint8_t r[7]={0x1F,0x10,0x1E,0x10,0x10,0x10,0x10}; memcpy(rows,r,7); return true; }
            case L'G': { uint8_t r[7]={0x0E,0x11,0x10,0x17,0x11,0x11,0x0E}; memcpy(rows,r,7); return true; }
            case L'H': { uint8_t r[7]={0x11,0x11,0x11,0x1F,0x11,0x11,0x11}; memcpy(rows,r,7); return true; }
            case L'I': { uint8_t r[7]={0x1F,0x04,0x04,0x04,0x04,0x04,0x1F}; memcpy(rows,r,7); return true; }
            case L'J': { uint8_t r[7]={0x1F,0x01,0x01,0x01,0x11,0x11,0x0E}; memcpy(rows,r,7); return true; }
            case L'K': { uint8_t r[7]={0x11,0x12,0x14,0x18,0x14,0x12,0x11}; memcpy(rows,r,7); return true; }
            case L'L': { uint8_t r[7]={0x10,0x10,0x10,0x10,0x10,0x10,0x1F}; memcpy(rows,r,7); return true; }
            case L'M': { uint8_t r[7]={0x11,0x1B,0x15,0x11,0x11,0x11,0x11}; memcpy(rows,r,7); return true; }
            case L'N': { uint8_t r[7]={0x11,0x19,0x15,0x13,0x11,0x11,0x11}; memcpy(rows,r,7); return true; }
            case L'O': { uint8_t r[7]={0x0E,0x11,0x11,0x11,0x11,0x11,0x0E}; memcpy(rows,r,7); return true; }
            case L'Q': { uint8_t r[7]={0x0E,0x11,0x11,0x11,0x15,0x12,0x0D}; memcpy(rows,r,7); return true; }
            case L'R': { uint8_t r[7]={0x1E,0x11,0x11,0x1E,0x12,0x11,0x11}; memcpy(rows,r,7); return true; }
            case L'T': { uint8_t r[7]={0x1F,0x04,0x04,0x04,0x04,0x04,0x04}; memcpy(rows,r,7); return true; }
            case L'U': { uint8_t r[7]={0x11,0x11,0x11,0x11,0x11,0x11,0x0E}; memcpy(rows,r,7); return true; }
            case L'V': { uint8_t r[7]={0x11,0x11,0x11,0x11,0x0A,0x0A,0x04}; memcpy(rows,r,7); return true; }
            case L'W': { uint8_t r[7]={0x11,0x11,0x11,0x11,0x15,0x1B,0x11}; memcpy(rows,r,7); return true; }
            case L'X': { uint8_t r[7]={0x11,0x11,0x0A,0x04,0x0A,0x11,0x11}; memcpy(rows,r,7); return true; }
            case L'Y': { uint8_t r[7]={0x11,0x11,0x0A,0x04,0x04,0x04,0x04}; memcpy(rows,r,7); return true; }
            case L'Z': { uint8_t r[7]={0x1F,0x01,0x02,0x04,0x08,0x10,0x1F}; memcpy(rows,r,7); return true; }
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
            case L':': { uint8_t r[7]={0x00,0x06,0x06,0x00,0x06,0x06,0x00}; memcpy(rows,r,7); return true; }
            case L' ': { uint8_t r[7]={0,0,0,0,0,0,0}; memcpy(rows,r,7); return true; }
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
        // Biome label under FPS
        if(m_world){
            XMFLOAT3 p = m_cam.Position();
            int wx = (int)std::floor(p.x);
            int wz = (int)std::floor(p.z);
            auto biome = m_world->GetBiomeAt(wx, wz);
            const wchar_t* bname = L"PLAINS";
            if(biome == VoxelWorld::Biome::Hills) bname = L"HILLS";
            else if(biome == VoxelWorld::Biome::Mountains) bname = L"MOUNTAINS";
            DrawGlyphStringTopRight(bname, 0xFFB0FFC0u, 2, 36); // slightly lower margin
        }
    }

    // ---- Menu/UI helpers ----
    void DrawGlyphStringAt(float x, float y, const std::wstring& text, uint32_t color, int scale){
        int charW = 6*scale;
        float pxSize = (float)scale;
        m_overlayVerts.clear();
        float sx = x;
        for(wchar_t ch : text){
            uint8_t rows[7];
            if(!GetGlyph5x7(ch, rows)) { sx += (float)charW; continue; }
            for(int r=0;r<7;++r){
                uint8_t bits = rows[r];
                for(int c=0;c<5;++c){
                    if(bits & (1 << (4-c))){
                        float rx = sx + c*pxSize;
                        float ry = y + r*pxSize;
                        OverlayVertex v0{{rx,         ry        }, color};
                        OverlayVertex v1{{rx+pxSize,  ry        }, color};
                        OverlayVertex v2{{rx+pxSize,  ry+pxSize }, color};
                        OverlayVertex v3{{rx,         ry+pxSize }, color};
                        m_overlayVerts.push_back(v0); m_overlayVerts.push_back(v1); m_overlayVerts.push_back(v2);
                        m_overlayVerts.push_back(v0); m_overlayVerts.push_back(v2); m_overlayVerts.push_back(v3);
                    }
                }
            }
            sx += (float)charW;
        }
        struct { float w; float h; float pad[2]; } cb{ (float)m_width, (float)m_height, {0,0} };
        GetContext()->UpdateSubresource(m_cbOverlay.Get(), 0, nullptr, &cb, 0, 0);
        if(!m_overlayVerts.empty()){
            D3D11_MAPPED_SUBRESOURCE mapped{};
            if(SUCCEEDED(GetContext()->Map(m_vbOverlay.Get(), 0, D3D11_MAP_WRITE_DISCARD, 0, &mapped))){
                memcpy(mapped.pData, m_overlayVerts.data(), m_overlayVerts.size()*sizeof(OverlayVertex));
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
                // Restore
                GetContext()->OMSetDepthStencilState(m_dss.Get(), 0);
                GetContext()->RSSetState(m_rsSolid.Get());
                GetContext()->IASetInputLayout(m_inputLayout.Get());
                GetContext()->VSSetShader(m_vs.Get(), nullptr, 0);
                GetContext()->PSSetShader(m_ps.Get(), nullptr, 0);
            }
        }
    }

    void DrawSolidRect(float x, float y, float w, float h, uint32_t color){
        m_overlayVerts.clear();
        OverlayVertex v0{{x, y}, color};
        OverlayVertex v1{{x+w, y}, color};
        OverlayVertex v2{{x+w, y+h}, color};
        OverlayVertex v3{{x, y+h}, color};
        m_overlayVerts.push_back(v0); m_overlayVerts.push_back(v1); m_overlayVerts.push_back(v2);
        m_overlayVerts.push_back(v0); m_overlayVerts.push_back(v2); m_overlayVerts.push_back(v3);
        struct { float w; float h; float pad[2]; } cb{ (float)m_width, (float)m_height, {0,0} };
        GetContext()->UpdateSubresource(m_cbOverlay.Get(), 0, nullptr, &cb, 0, 0);
        D3D11_MAPPED_SUBRESOURCE mapped{};
        if(SUCCEEDED(GetContext()->Map(m_vbOverlay.Get(), 0, D3D11_MAP_WRITE_DISCARD, 0, &mapped))){
            memcpy(mapped.pData, m_overlayVerts.data(), m_overlayVerts.size()*sizeof(OverlayVertex));
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
            // Restore
            GetContext()->OMSetDepthStencilState(m_dss.Get(), 0);
            GetContext()->RSSetState(m_rsSolid.Get());
            GetContext()->IASetInputLayout(m_inputLayout.Get());
            GetContext()->VSSetShader(m_vs.Get(), nullptr, 0);
            GetContext()->PSSetShader(m_ps.Get(), nullptr, 0);
        }
    }

    int NextCap(int current){ static int caps[] = {0,30,60,90,120,144,165,240,360}; int n=(int)(sizeof(caps)/sizeof(caps[0])); for(int i=0;i<n;++i){ if(caps[i]==current) return caps[(i+1)%n]; } return 60; }
    int PrevCap(int current){ static int caps[] = {0,30,60,90,120,144,165,240,360}; int n=(int)(sizeof(caps)/sizeof(caps[0])); for(int i=0;i<n;++i){ if(caps[i]==current) return caps[(i-1+n)%n]; } return 60; }

    void DrawSettingsMenu(){
        // Panel
    float w = 560.0f, h = 504.0f;
        float x = (m_width - w) * 0.5f;
        float y = (m_height - h) * 0.4f;
        DrawSolidRect(x, y, w, h, 0xAA000000u);

        int scaleTitle = 3;
        int scaleItem = 2;
        float ty = y + 20.0f;
        float tx = x + 20.0f;
    DrawGlyphStringAt(tx, ty, L"SETTINGS", 0xFFFFFFFFu, scaleTitle);
        ty += 40.0f;

    std::wstring vs = L"VSYNC: "; vs += (m_vsyncEnabled ? L"ON" : L"OFF");
        uint32_t color0 = (m_menuIndex==0 ? 0xFFFFFF00u : 0xFFFFFFFFu);
        DrawGlyphStringAt(tx, ty, vs, color0, scaleItem); ty += 28.0f;

    std::wstring msaa = L"MSAA: "; if(!m_msaaEnabled) msaa += L"OFF"; else { wchar_t sbuf[16]; swprintf(sbuf,16,L"%dx", m_msaaSamples); msaa += sbuf; }
        uint32_t colorMS = (m_menuIndex==1 ? 0xFFFFFF00u : 0xFFFFFFFFu);
        DrawGlyphStringAt(tx, ty, msaa, colorMS, scaleItem); ty += 28.0f;

    wchar_t capBuf[64]; swprintf(capBuf, 64, L"FPS CAP: %d%s", m_fpsCap, (m_fpsCap==0? L" (UNCAPPED)" : L""));
        uint32_t color1 = (m_menuIndex==2 ? 0xFFFFFF00u : 0xFFFFFFFFu);
        DrawGlyphStringAt(tx, ty, capBuf, color1, scaleItem); ty += 28.0f;

        // Render distance (in chunks)
        wchar_t rdBuf[64]; swprintf(rdBuf, 64, L"RENDER DIST: %d", m_viewRadius);
        uint32_t color2 = (m_menuIndex==3 ? 0xFFFFFF00u : 0xFFFFFFFFu);
        DrawGlyphStringAt(tx, ty, rdBuf, color2, scaleItem); ty += 28.0f;

    // Shadow preset and sliders
    static const wchar_t* presetNames[] = { L"FASTEST", L"FAST", L"BALANCED", L"ACCURATE" };
    int pi = ResolvePresetIndex();
    const wchar_t* pname = (pi>=0? presetNames[pi] : L"CUSTOM");
    uint32_t colorP = (m_menuIndex==4 ? 0xFFFFFF00u : 0xFFFFFFFFu);
    DrawGlyphStringAt(tx, ty, std::wstring(L"SHADOW PRESET: ") + pname, colorP, scaleItem); ty += 28.0f;

    // Samples slider
    uint32_t colorS = (m_menuIndex==5 ? 0xFFFFFF00u : 0xFFFFFFFFu);
    wchar_t sBuf[64]; swprintf(sBuf,64,L"SAMPLES: %d", std::max(1,std::min(8,m_mesherSettings.vertexSamples)));
    DrawGlyphStringAt(tx, ty, sBuf, colorS, scaleItem); ty += 28.0f;

    // Radius slider
    uint32_t colorR = (m_menuIndex==6 ? 0xFFFFFF00u : 0xFFFFFFFFu);
    wchar_t rBuf[64]; swprintf(rBuf,64,L"RADIUS: %.2f", m_mesherSettings.sampleRadius);
    DrawGlyphStringAt(tx, ty, rBuf, colorR, scaleItem); ty += 28.0f;

    // Distance slider (steps)
    uint32_t colorD = (m_menuIndex==7 ? 0xFFFFFF00u : 0xFFFFFFFFu);
    wchar_t dBuf[64]; swprintf(dBuf,64,L"DIST: %d", m_mesherSettings.maxSteps);
    DrawGlyphStringAt(tx, ty, dBuf, colorD, scaleItem); ty += 28.0f;

    // Subdivision toggle
    uint32_t colorT = (m_menuIndex==8 ? 0xFFFFFF00u : 0xFFFFFFFFu);
    DrawGlyphStringAt(tx, ty, std::wstring(L"SUBDIV: ") + (m_mesherSettings.enableSubdivision?L"ON":L"OFF"), colorT, scaleItem); ty += 28.0f;

    // Subdivision threshold
    uint32_t colorTh = (m_menuIndex==9 ? 0xFFFFFF00u : 0xFFFFFFFFu);
    wchar_t thBuf[64]; swprintf(thBuf,64,L"SUBDIV THRESH: %.2f", m_mesherSettings.subdivThreshold);
    DrawGlyphStringAt(tx, ty, thBuf, colorTh, scaleItem); ty += 28.0f;

    // Hard shadows toggle
    uint32_t colorHS = (m_menuIndex==10 ? 0xFFFFFF00u : 0xFFFFFFFFu);
    DrawGlyphStringAt(tx, ty, std::wstring(L"HARD SHADOWS: ") + (m_mesherSettings.hardShadows?L"ON":L"OFF"), colorHS, scaleItem); ty += 28.0f;

    // Sun samples for soft shadows
    uint32_t colorSS = (m_menuIndex==11 ? 0xFFFFFF00u : 0xFFFFFFFFu);
    wchar_t ssBuf[64]; swprintf(ssBuf,64,L"SUN SAMPLES: %d", std::max(1,m_mesherSettings.sunSamples));
    DrawGlyphStringAt(tx, ty, ssBuf, colorSS, scaleItem); ty += 28.0f;

    // Instruction line omitted
    }

    bool BuildChunkGPU(int cx, int cz, GPUChunk& out) {
        if(!IsWithinView(cx, cz)) return false;
        MesherSettings s = m_mesherSettings;
        if(m_inInitialLoading){
            // Fast path: full-bright, no subdivision, minimal work
            s.fastBakeOnly = true;
            s.enableSubdivision = false;
            s.vertexSamples = 1; s.maxSteps = 16; s.sunSamples = 1;
        } else {
            s.fastBakeOnly = false;
        }
        ChunkMesh mesh = VoxelMesher::BuildChunkMesh(*m_world, cx, cz, m_sunDir, s);
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
        // Generate list and sort by distance (closest first). Frustum-cull to avoid meshing unseen chunks.
    UpdateFrustumPlanes();
    std::vector<std::pair<int,int>> list;
        list.reserve((2*m_viewRadius+1)*(2*m_viewRadius+1));
        for(int dz=-m_viewRadius; dz<=m_viewRadius; ++dz){
            for(int dx=-m_viewRadius; dx<=m_viewRadius; ++dx){
                int cx = m_camChunkX + dx;
                int cz = m_camChunkZ + dz;
                if(!ChunkInFrustum(cx, cz)) continue;
                long long key = ChunkKey(cx,cz);
        if(m_chunks.find(key) != m_chunks.end()) continue;
        if(m_requested.count(key)) continue;
        list.emplace_back(cx,cz);
            }
        }
        auto sq = [](int a){ return a*a; };
        std::sort(list.begin(), list.end(), [&](auto&a, auto&b){
            int dax = a.first - m_camChunkX, daz = a.second - m_camChunkZ;
            int dbx = b.first - m_camChunkX, dbz = b.second - m_camChunkZ;
            return sq(dax)+sq(daz) > sq(dbx)+sq(dbz); // farthest first so we pop_back() closest
        });
    // Track requested and queue jobs to worker threads
    for(auto& c : list) m_requested.insert(ChunkKey(c.first,c.second));
    QueueMeshingJobs(list);
    }

    void CleanupFarChunks() {
        // Hysteresis: unload only when beyond viewRadius + margin
        const int margin = 2;
        auto beyond = [&](int cx, int cz){
            int dx = cx - m_camChunkX; int dz = cz - m_camChunkZ;
            return (std::max(std::abs(dx), std::abs(dz)) > (m_viewRadius + margin));
        };
        std::vector<long long> toRemove;
        for(const auto& kv : m_chunks) {
            int cx = (int)(kv.first >> 32);
            int cz = (int)(kv.first & 0xffffffff);
            if(beyond(cx, cz)) toRemove.push_back(kv.first);
        }
        for(auto key : toRemove){
            m_chunks.erase(key);
            m_requested.erase(key);
            m_lru.erase(key);
        }
        // LRU cap
        const size_t MaxResident = 2048;
        if(m_chunks.size() > MaxResident){
            std::vector<std::pair<long long, uint64_t>> candidates;
            candidates.reserve(m_lru.size());
            for(auto &p : m_lru){
                long long key = p.first; uint64_t stamp = p.second;
                int cx = (int)(key >> 32); int cz = (int)(key & 0xffffffff);
                if(!IsWithinView(cx, cz)) candidates.emplace_back(key, stamp);
            }
            std::sort(candidates.begin(), candidates.end(), [](auto&a, auto&b){ return a.second < b.second; });
            size_t need = m_chunks.size() - MaxResident;
            for(size_t i=0; i<need && i<candidates.size(); ++i){
                long long key = candidates[i].first;
                m_chunks.erase(key);
                m_requested.erase(key);
                m_lru.erase(key);
            }
        }
    }

    // ---------- Multithreaded meshing ----------
    void EnsureWorkers() {
        if(!m_workers.empty()) return;
        unsigned n = std::max(1u, std::thread::hardware_concurrency());
        m_stopWorkers = false;
        m_workers.reserve(n);
        for(unsigned i=0;i<n;++i){
            m_workers.emplace_back([this]{ this->WorkerLoop(); });
        }
    }

    void WorkerLoop() {
        for(;;){
            std::pair<int,int> job;
            {
                std::unique_lock<std::mutex> lk(m_jobMutex);
                m_jobCv.wait(lk, [&]{ return m_stopWorkers || !m_jobs.empty(); });
                if(m_stopWorkers && m_jobs.empty()) return;
                job = m_jobs.back(); m_jobs.pop_back();
            }
            // Build mesh on CPU
            int cx = job.first, cz = job.second;
            if(!IsWithinView(cx, cz)) continue;
            // Mesh cache: reuse if version unchanged
        unsigned v = m_world->GetChunkVersion(cx, cz);
            ChunkMesh mesh;
            {
                std::lock_guard<std::mutex> lk(m_cacheMutex);
                auto key = ChunkKey(cx, cz);
                auto it = m_meshCache.find(key);
                if(it != m_meshCache.end() && it->second.version == v){
                    mesh = it->second.mesh; // copy small vectors; fine for now
                } else {
            mesh = VoxelMesher::BuildChunkMesh(*m_world, cx, cz, m_sunDir, m_mesherSettings);
                    if(!mesh.indices.empty()){
                        m_meshCache[key] = CachedMesh{ mesh, v };
                        m_cacheLru[key] = ++m_frameCounter;
                        // Cap the cache
                        const size_t MaxCache = 4096;
                        if(m_meshCache.size() > MaxCache){
                            std::vector<std::pair<long long,uint64_t>> items;
                            items.reserve(m_cacheLru.size());
                            for(auto &p : m_cacheLru) items.emplace_back(p.first, p.second);
                            std::sort(items.begin(), items.end(), [](auto&a, auto&b){ return a.second < b.second; });
                            size_t need = m_meshCache.size() - MaxCache;
                            for(size_t i=0;i<need && i<items.size(); ++i){
                                m_meshCache.erase(items[i].first);
                                m_cacheLru.erase(items[i].first);
                            }
                        }
                    }
                }
            }
            if(mesh.indices.empty()) continue;
            // Push result
            {
                std::lock_guard<std::mutex> lk(m_readyMutex);
                m_ready.push_back({cx, cz, std::move(mesh)});
            }
        }
    }

    void QueueMeshingJobs(const std::vector<std::pair<int,int>>& list) {
        EnsureWorkers();
        {
            std::lock_guard<std::mutex> lk(m_jobMutex);
            for(auto& c : list) m_jobs.push_back(c);
        }
        m_jobCv.notify_all();
    }

    void UploadReadyMeshesPerFrame(int maxUploads) {
        int uploaded=0;
        std::vector<CPUChunkMesh> local;
        {
            std::lock_guard<std::mutex> lk(m_readyMutex);
            // Move a batch to avoid holding the lock
            while(uploaded < maxUploads && !m_ready.empty()){
                local.push_back(std::move(m_ready.back()));
                m_ready.pop_back();
                ++uploaded;
            }
        }
    for(auto& cm : local){
            long long key = ChunkKey(cm.cx, cm.cz);
            if(m_chunks.find(key) != m_chunks.end()) continue;
            GPUChunk g{};
            // Upload GPU buffers
            D3D11_BUFFER_DESC vbd{}; vbd.ByteWidth = (UINT)(cm.mesh.vertices.size()*sizeof(VoxelVertex)); vbd.BindFlags = D3D11_BIND_VERTEX_BUFFER; vbd.Usage = D3D11_USAGE_DEFAULT;
            D3D11_SUBRESOURCE_DATA vinit{ cm.mesh.vertices.data(), 0, 0 };
            if(FAILED(GetDevice()->CreateBuffer(&vbd, &vinit, g.vb.GetAddressOf()))) continue;
            D3D11_BUFFER_DESC ibd{}; ibd.ByteWidth = (UINT)(cm.mesh.indices.size()*sizeof(uint32_t)); ibd.BindFlags = D3D11_BIND_INDEX_BUFFER; ibd.Usage = D3D11_USAGE_DEFAULT;
            D3D11_SUBRESOURCE_DATA iinit{ cm.mesh.indices.data(), 0, 0 };
            if(FAILED(GetDevice()->CreateBuffer(&ibd, &iinit, g.ib.GetAddressOf()))) continue;
            g.indexCount = (UINT)cm.mesh.indices.size();
            m_chunks.emplace(key, std::move(g));
            // Mark as satisfied
            m_requested.erase(key);
            m_lru[key] = ++m_frameCounter;
            // Track initial loading progress
            if(m_inInitialLoading){
                if(m_initialLoadSet.erase(key) > 0){
                    ++m_initialLoaded;
                }
            }
        }
    }

    // ---------- Frustum culling ----------
    struct Plane { XMFLOAT4 p; };
    Plane m_frustum[6]; // left, right, bottom, top, near, far

    static void NormalizePlane(Plane& pl){
        XMVECTOR v = XMLoadFloat4(&pl.p);
        XMVECTOR n = XMPlaneNormalize(v);
        XMStoreFloat4(&pl.p, n);
    }

    void UpdateFrustumPlanes(){
        XMMATRIX vp = m_cam.ViewProj();
        // Extract planes from VP (row-major per XMMATRIX is fine via direct access)
        XMFLOAT4X4 m; XMStoreFloat4x4(&m, vp);
        // Left:  m14 + m11, Right: m14 - m11
        m_frustum[0].p = { m._14 + m._11, m._24 + m._21, m._34 + m._31, m._44 + m._41 };
        m_frustum[1].p = { m._14 - m._11, m._24 - m._21, m._34 - m._31, m._44 - m._41 };
        // Bottom: m14 + m12, Top: m14 - m12
        m_frustum[2].p = { m._14 + m._12, m._24 + m._22, m._34 + m._32, m._44 + m._42 };
        m_frustum[3].p = { m._14 - m._12, m._24 - m._22, m._34 - m._32, m._44 - m._42 };
        // Near:   m13, Far: m14 - m13
        m_frustum[4].p = { m._13, m._23, m._33, m._43 };
        m_frustum[5].p = { m._14 - m._13, m._24 - m._23, m._34 - m._33, m._44 - m._43 };
        for(int i=0;i<6;++i) NormalizePlane(m_frustum[i]);
    }

    bool AabbInFrustum(const XMFLOAT3& mn, const XMFLOAT3& mx) const {
        for(int i=0;i<6;++i){
            const XMFLOAT4& p = m_frustum[i].p;
            // Compute positive vertex for plane
            float x = (p.x >= 0) ? mx.x : mn.x;
            float y = (p.y >= 0) ? mx.y : mn.y;
            float z = (p.z >= 0) ? mx.z : mn.z;
            if(p.x*x + p.y*y + p.z*z + p.w < 0) return false;
        }
        return true;
    }

    bool ChunkInFrustum(int cx, int cz) const {
        float x0 = float(cx * CHUNK_SIZE);
        float z0 = float(cz * CHUNK_SIZE);
        XMFLOAT3 mn{ x0, 0.0f, z0 };
        XMFLOAT3 mx{ x0 + CHUNK_SIZE, float(CHUNK_HEIGHT), z0 + CHUNK_SIZE };
        return AabbInFrustum(mn, mx);
    }

    LRESULT MsgProc(HWND hwnd, UINT msg, WPARAM wParam, LPARAM lParam) override {
        switch(msg) {
        case WM_CHAR: {
                if(m_needSeedInput){
                    wchar_t ch = (wchar_t)wParam;
                    if(ch >= L'0' && ch <= L'9'){
                        if(m_seedInput.size() < 16) m_seedInput.push_back(ch);
                    } else if(ch == L'-'){
                        if(m_seedInput.empty()) m_seedInput.push_back(ch);
                    }
                    return 0;
                }
                break;
            }
        case WM_KEYDOWN: {
                int vk = (int)wParam;
                if(vk >= 0 && vk < 256) m_keys[vk] = true;
                if(m_needSeedInput){
                    if(vk == VK_RETURN){
                        unsigned seed = 1337u;
                        try{
                            if(!m_seedInput.empty()){
                                unsigned long v = std::stoul(std::wstring(m_seedInput.begin(), m_seedInput.end()));
                                seed = (unsigned)v;
                            }
                        } catch(...){ seed = 1337u; }
                        m_world = std::make_unique<VoxelWorld>(seed);
                        InitializePlayer();
                        UpdateCameraChunk();
                        EnqueueVisibleChunks();
                        PrepareInitialLoad();
                        m_needSeedInput = false;
                        return 0;
                    } else if(vk == VK_BACK){
                        if(!m_seedInput.empty()) m_seedInput.pop_back();
                        return 0;
                    } else if(vk == VK_ESCAPE){
                        m_showMenu = !m_showMenu;
                        if(m_showMenu) ReleaseMouse(); else CaptureMouse();
                        return 0;
                    }
                } else if(vk == VK_ESCAPE){
                    // Toggle settings menu
                    m_showMenu = !m_showMenu;
            if(m_showMenu) ReleaseMouse(); else CaptureMouse();
                } else if(vk == VK_OEM_3) {
                    // Toggle freecam with backtick (`)
                    m_freecam = !m_freecam;
                    m_showPlayer = !m_freecam;
                } else if(m_showMenu) {
                    // Menu navigation when open
                    if(vk == VK_UP){ m_menuIndex = (m_menuIndex + (m_menuItems-1)) % m_menuItems; }
                    else if(vk == VK_DOWN){ m_menuIndex = (m_menuIndex + 1) % m_menuItems; }
                    else if(vk == VK_LEFT || vk == VK_RIGHT || vk == VK_RETURN){
                        if(m_menuIndex == 0){ // vsync toggle
                            m_vsyncEnabled = !m_vsyncEnabled;
                        } else if(m_menuIndex == 1){ // MSAA toggle/cycle
                            if(vk == VK_RETURN){ m_msaaEnabled = !m_msaaEnabled; CreateMsaaTargets(); }
                            else if(vk == VK_LEFT){ if(m_msaaSamples>2) m_msaaSamples = (m_msaaSamples==8?4:2); m_msaaEnabled = true; CreateMsaaTargets(); }
                            else if(vk == VK_RIGHT){ if(m_msaaSamples<8) m_msaaSamples = (m_msaaSamples==2?4:8); m_msaaEnabled = true; CreateMsaaTargets(); }
                        } else if(m_menuIndex == 2){ // fps cap
                            if(vk == VK_RETURN){ m_fpsCap = 0; }
                            else if(vk == VK_LEFT){ m_fpsCap = PrevCap(m_fpsCap); }
                            else if(vk == VK_RIGHT){ m_fpsCap = NextCap(m_fpsCap); }
                        } else if(m_menuIndex == 3){ // render distance
                            const int maxR = 64; // cap per requirements
                            if(vk == VK_LEFT && m_viewRadius > 1) { m_viewRadius--; }
                            else if(vk == VK_RIGHT && m_viewRadius < maxR) { m_viewRadius++; }
                            // Reset queues to reflect new radius
                            m_requested.clear();
                            CleanupFarChunks();
                            EnqueueVisibleChunks();
                        } else if(m_menuIndex == 4){ // shadow preset
                            int pi = ResolvePresetIndex();
                            if(vk == VK_LEFT){ ApplyShadowPreset(std::max(0, (pi>=0?pi:2) - 1)); }
                            else if(vk == VK_RIGHT){ ApplyShadowPreset(std::min(3, (pi>=0?pi:2) + 1)); }
                            else if(vk == VK_RETURN){ ApplyShadowPreset(2); }
                        } else if(m_menuIndex == 5){ // samples
                            int s = m_mesherSettings.vertexSamples;
                            if(vk == VK_LEFT){ s = (s==8?4:(s==4?2:1)); }
                            else if(vk == VK_RIGHT){ s = (s==1?2:(s==2?4:8)); }
                            m_mesherSettings.vertexSamples = s; OnShadowSettingsChanged();
                        } else if(m_menuIndex == 6){ // radius
                            float step = 0.01f;
                            if(vk == VK_LEFT) m_mesherSettings.sampleRadius = std::max(0.05f, m_mesherSettings.sampleRadius - step);
                            else if(vk == VK_RIGHT) m_mesherSettings.sampleRadius = std::min(0.50f, m_mesherSettings.sampleRadius + step);
                            OnShadowSettingsChanged();
                        } else if(m_menuIndex == 7){ // distance
                            int stepN = 16;
                            if(vk == VK_LEFT) m_mesherSettings.maxSteps = std::max(16, m_mesherSettings.maxSteps - stepN);
                            else if(vk == VK_RIGHT) m_mesherSettings.maxSteps = std::min(256, m_mesherSettings.maxSteps + stepN);
                            OnShadowSettingsChanged();
                        } else if(m_menuIndex == 8){ // subdiv toggle
                            if(vk == VK_RETURN || vk == VK_LEFT || vk == VK_RIGHT){ m_mesherSettings.enableSubdivision = !m_mesherSettings.enableSubdivision; OnShadowSettingsChanged(); }
                        } else if(m_menuIndex == 9){ // subdiv threshold
                            float step = 0.05f;
                            if(vk == VK_LEFT) m_mesherSettings.subdivThreshold = std::max(0.05f, m_mesherSettings.subdivThreshold - step);
                            else if(vk == VK_RIGHT) m_mesherSettings.subdivThreshold = std::min(0.80f, m_mesherSettings.subdivThreshold + step);
                            OnShadowSettingsChanged();
                        } else if(m_menuIndex == 10){ // hard shadows toggle
                            if(vk == VK_RETURN || vk == VK_LEFT || vk == VK_RIGHT){ m_mesherSettings.hardShadows = !m_mesherSettings.hardShadows; OnShadowSettingsChanged(); }
                        } else if(m_menuIndex == 11){ // sun samples
                            int s = std::max(1, m_mesherSettings.sunSamples);
                            if(vk == VK_LEFT){ s = (s<=1?1:(s<=2?1:(s<=4?2:(s<=8?4:(s<=16?8:16))))); }
                            else if(vk == VK_RIGHT){ s = (s<1?1:(s<2?2:(s<4?4:(s<8?8:16)))); }
                            m_mesherSettings.sunSamples = s; OnShadowSettingsChanged();
                        }
                    }
                }
                break;
            }
            case WM_KEYUP: {
                int vk = (int)wParam;
                if(vk >= 0 && vk < 256) m_keys[vk] = false;
                break;
            }
            // Right-click capture no longer required; mouse is locked by default
            case WM_MOUSEMOVE: {
                if(m_mouseCaptured && !m_showMenu) {
                    int mx = GET_X_LPARAM(lParam);
                    int my = GET_Y_LPARAM(lParam);
                    if(m_warpingMouse){
                        m_warpingMouse = false;
                        m_lastMouse.x = mx; m_lastMouse.y = my;
                        break;
                    }
                    int dx = mx - m_lastMouse.x;
                    int dy = my - m_lastMouse.y;
                    m_lastMouse.x = mx;
                    m_lastMouse.y = my;

                    float yaw = dx * m_mouseSens;
                    float pitch = dy * m_mouseSens; // inverted per user: moving mouse up looks down and vice versa

                    // clamp pitch to avoid flipping
                    float newPitch = m_pitchAngle + pitch;
                    const float limit = 1.2f; // ~69 deg
                    float applyPitch = pitch;
                    if(newPitch > limit) applyPitch = limit - m_pitchAngle;
                    else if(newPitch < -limit) applyPitch = -limit - m_pitchAngle;
                    m_pitchAngle += applyPitch;

                    m_cam.RotateY(yaw);
                    m_cam.Pitch(applyPitch);

                    // Recenter cursor to avoid hitting window edges
                    RECT rc{}; GetClientRect(hwnd, &rc);
                    POINT tl{rc.left, rc.top}; ClientToScreen(hwnd, &tl);
                    int cx = tl.x + (rc.right - rc.left)/2;
                    int cy = tl.y + (rc.bottom - rc.top)/2;
                    m_warpingMouse = true;
                    SetCursorPos(cx, cy);
                    m_lastMouse.x = (rc.right - rc.left)/2;
                    m_lastMouse.y = (rc.bottom - rc.top)/2;
                }
                break;
            }
            case WM_SETFOCUS: {
                if(!m_showMenu) CaptureMouse();
                break;
            }
            case WM_KILLFOCUS: {
                ReleaseMouse();
                break;
            }
        }
        return D3DApp::MsgProc(hwnd, msg, wParam, lParam);
    }

    Camera m_cam;
    std::unique_ptr<VoxelWorld> m_world;
    std::array<bool,256> m_keys{};
    bool m_mouseCaptured = false;
    bool m_warpingMouse = false;
    POINT m_lastMouse{0,0};
    float m_moveSpeed = 10.0f; // units per second
    float m_mouseSens = 0.0025f; // radians per pixel
    float m_pitchAngle = 0.0f;
    // FPS
    float m_fpsAccum = 0.0f;
    int m_fpsFrames = 0;
    float m_lastDt = 0.0f;
    // High-precision frame timing for FPS cap
    long long m_qpcFreq = 0;
    long long m_frameStartQpc = 0;

    // Overlay pipeline
    ComPtr<ID3D11VertexShader> m_vsOverlay;
    ComPtr<ID3D11PixelShader> m_psOverlay;
    ComPtr<ID3D11InputLayout> m_layoutOverlay;
    ComPtr<ID3D11Buffer> m_vbOverlay;
    ComPtr<ID3D11Buffer> m_cbOverlay;
    ComPtr<ID3D11DepthStencilState> m_dssDisabled;
    ComPtr<ID3D11RasterizerState> m_rsOverlay;
    std::vector<OverlayVertex> m_overlayVerts;

    // Sky pipeline
    ComPtr<ID3D11VertexShader> m_vsSky;
    ComPtr<ID3D11PixelShader>  m_psSky;
    ComPtr<ID3D11InputLayout>  m_layoutSky;
    ComPtr<ID3D11Buffer>       m_vbSky;
    ComPtr<ID3D11Buffer>       m_cbSky;
    XMFLOAT3 m_sunDir{0.5f, 1.0f, 0.2f};
    float    m_sunSize = 0.009f; // ~0.5 degrees
    XMFLOAT3 m_skyHorizon{0.75f, 0.85f, 0.95f};
    XMFLOAT3 m_skyZenith{0.25f, 0.45f, 0.85f};

    // Player model
    ModelMesh m_playerMesh;
    ComPtr<ID3D11Buffer> m_playerVB;
    ComPtr<ID3D11Buffer> m_playerIB;
    UINT m_playerIndexCount = 0;
    ComPtr<ID3D11VertexShader> m_vsModel;
    ComPtr<ID3D11PixelShader>  m_psModel;
    ComPtr<ID3D11InputLayout>  m_layoutModel;
    ComPtr<ID3D11Buffer>       m_cbModel;
    bool m_showPlayer = true; // show in third-person
    bool m_freecam = false;   // free-fly camera toggle
    float m_freecamSpeed = 15.0f;
    // Model fit
    float m_modelScale = 1.0f;
    float m_modelMinY = 0.0f;
    float m_modelMaxY = 1.0f;
    // Player physics
    XMFLOAT3 m_playerPos{8.0f, 60.0f, 8.0f};
    XMFLOAT3 m_playerVel{0,0,0};
    XMFLOAT3 m_playerHalf{0.3f, 0.9f, 0.3f}; // half extents (width 0.6, height 1.8)
    float    m_eyeHeight = 0.75f; // camera offset from player origin
    float    m_gravity   = -25.0f;
    bool     m_onGround  = false;
    bool     m_jumpHeld  = false;
    float    m_thirdPersonDist = 6.0f; // camera distance behind player

    std::unordered_map<long long, GPUChunk> m_chunks;
    std::unordered_set<long long> m_requested;
    std::unordered_map<long long, uint64_t> m_lru;
    // CPU mesh cache for unchanged chunks
    std::unordered_map<long long, CachedMesh> m_meshCache;
    std::unordered_map<long long, uint64_t> m_cacheLru;
    std::mutex m_cacheMutex;
    // Meshing workers
    std::vector<std::thread> m_workers;
    std::vector<std::pair<int,int>> m_jobs;
    std::mutex m_jobMutex;
    std::condition_variable m_jobCv;
    bool m_stopWorkers = false;
    // Ready meshes from workers
    std::vector<CPUChunkMesh> m_ready;
    std::mutex m_readyMutex;
    int m_viewRadius = 24; // default per requirements
    int m_camChunkX = 0, m_camChunkZ = 0;
    uint64_t m_frameCounter = 0;

    ComPtr<ID3D11VertexShader> m_vs;
    ComPtr<ID3D11PixelShader> m_ps;
    ComPtr<ID3D11InputLayout> m_inputLayout;
    ComPtr<ID3D11Buffer> m_cbVP;
    ComPtr<ID3D11RasterizerState> m_rsSolid;
    ComPtr<ID3D11DepthStencilState> m_dss;

    float m_angle = 0.0f;

    // Settings
    bool m_vsyncEnabled = true;
    int m_fpsCap = 0; // 0 = uncapped
    bool m_showMenu = false;
    int m_menuIndex = 0; // 0: vsync, 1: MSAA, 2: FPS Cap, 3: Render Distance, 4: Shadow Preset, 5: Samples, 6: Radius, 7: Distance, 8: Subdiv Toggle, 9: Subdiv Threshold, 10: Hard Shadows, 11: Sun Samples
    int m_menuItems = 12;

    // Startup seed entry
    bool m_needSeedInput = false;
    std::wstring m_seedInput;

    // TAA removed

    // Meshing / shadow settings
    MesherSettings m_mesherSettings{};

    // --------- Initial Loading ---------
    bool m_inInitialLoading = false;
    int  m_initialTotal = 0;
    int  m_initialLoaded = 0;
    std::unordered_set<long long> m_initialLoadSet; // keys to wait for

    void PrepareInitialLoad(){
        // Build list of all chunks within current view radius around camera, regardless of frustum
        EnsureWorkers();
        UpdateCameraChunk();
        std::vector<std::pair<int,int>> list;
        list.reserve((2*m_viewRadius+1)*(2*m_viewRadius+1));
        for(int dz=-m_viewRadius; dz<=m_viewRadius; ++dz){
            for(int dx=-m_viewRadius; dx<=m_viewRadius; ++dx){
                int cx = m_camChunkX + dx;
                int cz = m_camChunkZ + dz;
                long long key = ChunkKey(cx,cz);
                if(m_chunks.find(key) != m_chunks.end()) continue;
                if(m_requested.count(key)) continue;
                list.emplace_back(cx,cz);
                m_initialLoadSet.insert(key);
            }
        }
        m_initialTotal = (int)m_initialLoadSet.size();
        m_initialLoaded = 0;
    if(m_initialTotal > 0){
            QueueMeshingJobs(list);
            m_inInitialLoading = true;
            // Ensure mouse is released during loading so the cursor is visible if needed
            ReleaseMouse();
        }
    }

    void DrawLoadingScreen(){
        // Draw a centered panel with a progress bar using overlay pipeline helpers
        float w = 520.0f, h = 140.0f;
        float x = (m_width - w) * 0.5f;
        float y = (m_height - h) * 0.55f;
        DrawSolidRect(x, y, w, h, 0xBB000000u);
        // Title
        DrawGlyphStringAt(x + 20.0f, y + 16.0f, L"LOADING WORLD", 0xFFFFFFFFu, 3);
        // Bar background
        float bx = x + 20.0f; float bw = w - 40.0f; float by = y + 72.0f; float bh = 24.0f;
        DrawSolidRect(bx, by, bw, bh, 0xFF333333u);
        // Progress fill
        float p = (m_initialTotal>0) ? (float)m_initialLoaded / (float)m_initialTotal : 0.0f;
        p = std::max(0.0f, std::min(1.0f, p));
        float fw = bw * p;
        DrawSolidRect(bx, by, fw, bh, 0xFF66CC66u);
        // Percent text
        wchar_t pctBuf[32]; int pct = (int)std::round(p * 100.0f); swprintf(pctBuf, 32, L"%d%%", pct);
        // Center percent over the bar (approximate centering)
        float tx = x + (w - (float)6*2 /*char width*/ * (float)wcslen(pctBuf)) * 0.5f; // coarse center
        DrawGlyphStringAt(tx, by + 32.0f, pctBuf, 0xFFFFFFFFu, 2);
    }

    void DrawSeedInputScreen(){
        float w = 560.0f, h = 180.0f;
        float x = (m_width - w) * 0.5f;
        float y = (m_height - h) * 0.45f;
        DrawSolidRect(x, y, w, h, 0xCC000000u);
        DrawGlyphStringAt(x + 20.0f, y + 16.0f, L"ENTER SEED", 0xFFFFFFFFu, 3);
        float bx = x + 20.0f, bw = w - 40.0f, by = y + 72.0f, bh = 28.0f;
        DrawSolidRect(bx, by, bw, bh, 0xFF222222u);
        std::wstring label = m_seedInput.empty() ? L"(numbers)" : m_seedInput;
        DrawGlyphStringAt(bx + 8.0f, by + 6.0f, label, 0xFFE0F0FFu, 2);
        DrawGlyphStringAt(x + 20.0f, y + h - 44.0f, L"ENTER: START   BACKSPACE: DELETE   ESC: MENU", 0xFFB0FFC0u, 2);
    }


    void OnShadowSettingsChanged(){
        // Invalidate caches and queued results; rebuild visible chunks
        {
            std::lock_guard<std::mutex> lk(m_cacheMutex);
            m_meshCache.clear();
            m_cacheLru.clear();
        }
        m_chunks.clear();
        m_requested.clear();
        {
            std::lock_guard<std::mutex> lk(m_readyMutex);
            m_ready.clear();
        }
        EnqueueVisibleChunks();
    }

    void ApplyShadowPreset(int idx){
        // 0 FASTEST, 1 FAST, 2 BALANCED, 3 ACCURATE
        if(idx <= 0){
            m_mesherSettings.maxSteps = 48;
            m_mesherSettings.vertexSamples = 1;
            m_mesherSettings.sampleRadius = 0.10f;
            m_mesherSettings.normalBias = 0.02f;
            m_mesherSettings.sunBias = 0.02f;
            m_mesherSettings.inset = 0.05f;
            m_mesherSettings.minShadow = 0.10f;
        } else if(idx == 1){
            m_mesherSettings.maxSteps = 64;
            m_mesherSettings.vertexSamples = 2;
            m_mesherSettings.sampleRadius = 0.12f;
            m_mesherSettings.normalBias = 0.02f;
            m_mesherSettings.sunBias = 0.02f;
            m_mesherSettings.inset = 0.06f;
            m_mesherSettings.minShadow = 0.08f;
        } else if(idx == 2){
            m_mesherSettings.maxSteps = 96;
            m_mesherSettings.vertexSamples = 4;
            m_mesherSettings.sampleRadius = 0.15f;
            m_mesherSettings.normalBias = 0.02f;
            m_mesherSettings.sunBias = 0.02f;
            m_mesherSettings.inset = 0.08f;
            m_mesherSettings.minShadow = 0.05f;
        } else { // 3 ACCURATE
            m_mesherSettings.maxSteps = 144;
            m_mesherSettings.vertexSamples = 8;
            m_mesherSettings.sampleRadius = 0.18f;
            m_mesherSettings.normalBias = 0.02f;
            m_mesherSettings.sunBias = 0.02f;
            m_mesherSettings.inset = 0.10f;
            m_mesherSettings.minShadow = 0.04f;
            m_mesherSettings.enableSubdivision = true;
            m_mesherSettings.subdivThreshold = 0.18f;
        }
        OnShadowSettingsChanged();
    }

    int ResolvePresetIndex() const {
        // Roughly match against presets; if none match exactly, return -1 (CUSTOM)
        struct P{int steps,samples; float radius,inset,minS;};
        P p[4] = {
            {48,1,0.10f,0.05f,0.10f},
            {64,2,0.12f,0.06f,0.08f},
            {96,4,0.15f,0.08f,0.05f},
            {144,8,0.18f,0.10f,0.04f}
        };
        for(int i=0;i<4;++i){
            if(m_mesherSettings.maxSteps==p[i].steps && m_mesherSettings.vertexSamples==p[i].samples &&
               fabsf(m_mesherSettings.sampleRadius - p[i].radius) < 1e-3f && fabsf(m_mesherSettings.inset - p[i].inset) < 1e-3f &&
               fabsf(m_mesherSettings.minShadow - p[i].minS) < 1e-3f) return i;
        }
        return -1;
    }

    // MSAA
    bool m_msaaEnabled = false;
    int  m_msaaSamples = 4;
    ComPtr<ID3D11Texture2D> m_msaaTex;
    ComPtr<ID3D11RenderTargetView> m_msaaRTV;
    ComPtr<ID3D11Texture2D> m_msaaDepthTex;
    ComPtr<ID3D11DepthStencilView> m_msaaDSV;

    void CreateMsaaTargets(){
        // Release old
        m_msaaTex.Reset(); m_msaaRTV.Reset(); m_msaaDepthTex.Reset(); m_msaaDSV.Reset();
        if(!m_msaaEnabled) return;
        // Check quality levels and adjust samples if needed
        UINT quality = 0; int samples = m_msaaSamples;
        while(samples > 1){
            if(SUCCEEDED(GetDevice()->CheckMultisampleQualityLevels(DXGI_FORMAT_R8G8B8A8_UNORM, samples, &quality)) && quality > 0) break;
            samples /= 2;
        }
        if(samples <= 1){ m_msaaEnabled = false; return; }
        m_msaaSamples = samples; UINT q = quality - 1;
        // Color
        D3D11_TEXTURE2D_DESC td{}; td.Width=m_width; td.Height=m_height; td.MipLevels=1; td.ArraySize=1; td.Format=DXGI_FORMAT_R8G8B8A8_UNORM; td.SampleDesc.Count=m_msaaSamples; td.SampleDesc.Quality=q; td.Usage=D3D11_USAGE_DEFAULT; td.BindFlags=D3D11_BIND_RENDER_TARGET;
        if(FAILED(GetDevice()->CreateTexture2D(&td, nullptr, m_msaaTex.GetAddressOf()))){ m_msaaEnabled=false; return; }
        GetDevice()->CreateRenderTargetView(m_msaaTex.Get(), nullptr, m_msaaRTV.GetAddressOf());
        // Depth
        D3D11_TEXTURE2D_DESC dd{}; dd.Width=m_width; dd.Height=m_height; dd.MipLevels=1; dd.ArraySize=1; dd.Format=DXGI_FORMAT_D24_UNORM_S8_UINT; dd.SampleDesc.Count=m_msaaSamples; dd.SampleDesc.Quality=q; dd.Usage=D3D11_USAGE_DEFAULT; dd.BindFlags=D3D11_BIND_DEPTH_STENCIL;
        if(FAILED(GetDevice()->CreateTexture2D(&dd, nullptr, m_msaaDepthTex.GetAddressOf()))){ m_msaaEnabled=false; m_msaaTex.Reset(); m_msaaRTV.Reset(); return; }
        GetDevice()->CreateDepthStencilView(m_msaaDepthTex.Get(), nullptr, m_msaaDSV.GetAddressOf());
    }

    // TAA removed: no post-process helpers

    // --------- Player Model ----------
    bool CreateModelPipeline(){
        UINT flags = D3DCOMPILE_ENABLE_STRICTNESS;
#ifdef _DEBUG
        flags |= D3DCOMPILE_DEBUG;
#endif
        ComPtr<ID3DBlob> vsb, psb, errors;
        if(FAILED(D3DCompileFromFile(L"shaders/model.hlsl", nullptr, nullptr, "VSMain", "vs_5_0", flags, 0, &vsb, &errors))){
            if(errors) MessageBoxA(GetHwnd(), (const char*)errors->GetBufferPointer(), "Model VS Compile Error", MB_OK|MB_ICONERROR);
            return false;
        }
        if(FAILED(D3DCompileFromFile(L"shaders/model.hlsl", nullptr, nullptr, "PSMain", "ps_5_0", flags, 0, &psb, &errors))){
            if(errors) MessageBoxA(GetHwnd(), (const char*)errors->GetBufferPointer(), "Model PS Compile Error", MB_OK|MB_ICONERROR);
            return false;
        }
        GetDevice()->CreateVertexShader(vsb->GetBufferPointer(), vsb->GetBufferSize(), nullptr, m_vsModel.GetAddressOf());
        GetDevice()->CreatePixelShader(psb->GetBufferPointer(), psb->GetBufferSize(), nullptr, m_psModel.GetAddressOf());
        D3D11_INPUT_ELEMENT_DESC il[] = {
            {"POSITION",0,DXGI_FORMAT_R32G32B32_FLOAT,0,0, D3D11_INPUT_PER_VERTEX_DATA,0},
            {"NORMAL",  0,DXGI_FORMAT_R32G32B32_FLOAT,0,12,D3D11_INPUT_PER_VERTEX_DATA,0},
        };
        GetDevice()->CreateInputLayout(il, _countof(il), vsb->GetBufferPointer(), vsb->GetBufferSize(), m_layoutModel.GetAddressOf());
        // cbuffer (matrix + color block + sun block + optional light matrix)
    D3D11_BUFFER_DESC cbd{}; cbd.BindFlags = D3D11_BIND_CONSTANT_BUFFER; cbd.ByteWidth = 160; cbd.Usage = D3D11_USAGE_DEFAULT;
        GetDevice()->CreateBuffer(&cbd, nullptr, m_cbModel.GetAddressOf());
        return true;
    }

    void LoadPlayerModel(){
        if(!CreateModelPipeline()) return;
        // Ensure assets are available; expect assets/models/player.obj
        std::wstring path = L"assets/models/player.obj";
        ModelMesh mm;
        if(!LoadOBJ(path, mm)){
            // Fallback: skip
            return;
        }
        m_playerMesh = std::move(mm);
        if(m_playerMesh.indices.empty()) return;
        // Create buffers
        D3D11_BUFFER_DESC vbd{}; vbd.ByteWidth = (UINT)(m_playerMesh.vertices.size()*sizeof(ModelVertex)); vbd.BindFlags=D3D11_BIND_VERTEX_BUFFER; vbd.Usage=D3D11_USAGE_DEFAULT;
        D3D11_SUBRESOURCE_DATA vinit{ m_playerMesh.vertices.data(), 0, 0 };
        if(FAILED(GetDevice()->CreateBuffer(&vbd, &vinit, m_playerVB.GetAddressOf()))) return;
        D3D11_BUFFER_DESC ibd{}; ibd.ByteWidth = (UINT)(m_playerMesh.indices.size()*sizeof(uint32_t)); ibd.BindFlags=D3D11_BIND_INDEX_BUFFER; ibd.Usage=D3D11_USAGE_DEFAULT;
        D3D11_SUBRESOURCE_DATA iinit{ m_playerMesh.indices.data(), 0, 0 };
        if(FAILED(GetDevice()->CreateBuffer(&ibd, &iinit, m_playerIB.GetAddressOf()))) return;
        m_playerIndexCount = (UINT)m_playerMesh.indices.size();

        // Compute bounds and fit to desired height (about 1.8m)
        float minY = std::numeric_limits<float>::infinity();
        float maxY = -std::numeric_limits<float>::infinity();
        for(const auto& v : m_playerMesh.vertices){
            minY = std::min(minY, v.pos.y);
            maxY = std::max(maxY, v.pos.y);
        }
        float modelH = std::max(0.0001f, maxY - minY);
        float desiredH = 1.8f; // match collision height
        m_modelScale = desiredH / modelH;
        m_modelMinY = minY;
        m_modelMaxY = maxY;
    }

    void DrawPlayer(){
        if(!m_showPlayer || m_playerIndexCount==0) return;
    // Determine model yaw from camera facing
        XMFLOAT3 look = m_cam.Look();
        float yaw = std::atan2f(look.x, look.z);
    // Compute scale to fit collision height and align base to player feet
    float scale = m_modelScale;
    float footY = m_playerPos.y - m_playerHalf.y; // ground contact point
    float baseOffsetY = -m_modelMinY * scale;     // raise model so minY touches y=0 in local
    float worldY = footY + baseOffsetY;
    // World matrix with scale and yaw rotation
    XMMATRIX S = XMMatrixScaling(scale, scale, scale);
        XMMATRIX R = XMMatrixRotationY(yaw);
    XMMATRIX T = XMMatrixTranslation(m_playerPos.x, worldY, m_playerPos.z);
        XMMATRIX W = S * R * T;
        XMMATRIX WVP = XMMatrixTranspose(W * m_cam.ViewProj());
    // Model cbuffer layout matches shaders/model.hlsl
    struct CBModel { XMMATRIX wvp; XMFLOAT3 color; float pad; XMFLOAT3 sunDir; float pad2; XMMATRIX lightVP; } cb;
    cb.wvp = WVP; cb.color = XMFLOAT3(0.9f,0.85f,0.8f); cb.pad = 0.0f; cb.sunDir = m_sunDir; cb.pad2 = 0.0f; cb.lightVP = XMMatrixIdentity();
        GetContext()->UpdateSubresource(m_cbModel.Get(), 0, nullptr, &cb, 0, 0);

        UINT stride = sizeof(ModelVertex), offset = 0; ID3D11Buffer* vb = m_playerVB.Get();
        GetContext()->OMSetDepthStencilState(m_dss.Get(), 0);
        GetContext()->RSSetState(m_rsSolid.Get());
        GetContext()->IASetInputLayout(m_layoutModel.Get());
        GetContext()->IASetPrimitiveTopology(D3D11_PRIMITIVE_TOPOLOGY_TRIANGLELIST);
        GetContext()->VSSetShader(m_vsModel.Get(), nullptr, 0);
        GetContext()->PSSetShader(m_psModel.Get(), nullptr, 0);
    GetContext()->VSSetConstantBuffers(0, 1, m_cbModel.GetAddressOf());
    GetContext()->PSSetConstantBuffers(0, 1, m_cbModel.GetAddressOf());
        GetContext()->IASetVertexBuffers(0, 1, &vb, &stride, &offset);
        GetContext()->IASetIndexBuffer(m_playerIB.Get(), DXGI_FORMAT_R32_UINT, 0);
    GetContext()->DrawIndexed(m_playerIndexCount, 0, 0);
        // restore voxel shaders for any subsequent draws
        GetContext()->IASetInputLayout(m_inputLayout.Get());
        GetContext()->VSSetShader(m_vs.Get(), nullptr, 0);
        GetContext()->PSSetShader(m_ps.Get(), nullptr, 0);
    }


    // --------- Mouse capture helpers ---------
    void UpdateMouseLockRect(){
        RECT rc{}; GetClientRect(GetHwnd(), &rc);
        POINT tl{rc.left, rc.top}; POINT br{rc.right, rc.bottom};
        ClientToScreen(GetHwnd(), &tl); ClientToScreen(GetHwnd(), &br);
        RECT clip{ tl.x, tl.y, br.x, br.y };
        ClipCursor(&clip);
    }
    void CaptureMouse(){
        if(m_mouseCaptured) return;
        ShowCursor(FALSE);
        m_mouseCaptured = true;
        UpdateMouseLockRect();
        // center cursor
        RECT rc{}; GetClientRect(GetHwnd(), &rc);
        POINT tl{rc.left, rc.top}; ClientToScreen(GetHwnd(), &tl);
        int cx = tl.x + (rc.right - rc.left)/2;
        int cy = tl.y + (rc.bottom - rc.top)/2;
        m_warpingMouse = true;
        SetCursorPos(cx, cy);
        m_lastMouse.x = (rc.right - rc.left)/2;
        m_lastMouse.y = (rc.bottom - rc.top)/2;
    }
    void ReleaseMouse(){
        if(!m_mouseCaptured) return;
        ClipCursor(nullptr);
        ShowCursor(TRUE);
        m_mouseCaptured = false;
    }

    // ---------- Player init and collision ----------
    void InitializePlayer(){
        // Find ground at initial XZ and place player slightly above
        int wx = (int)std::floor(m_playerPos.x);
        int wz = (int)std::floor(m_playerPos.z);
        int groundY = FindGroundY(wx, wz);
        if(groundY >= 0){ m_playerPos.y = (float)groundY + 2.0f; }
    // Camera position will be set by Update() in third-person
    }

    int FindGroundY(int wx, int wz){
        // search downward for first solid block from top of chunk height
        for(int y = CHUNK_HEIGHT-1; y >= 0; --y){
            if(m_world->GetBlock(wx, y, wz) != 0) return y;
        }
        return -1;
    }

    bool IsSolid(int x, int y, int z){
        if(y < 0) return true; // treat below world as solid
        if(y >= CHUNK_HEIGHT) return false;
        return m_world->GetBlock(x,y,z) != 0;
    }

    void CollideAndSlide(float dt){
        m_onGround = false;
        XMFLOAT3 pos = m_playerPos;
        XMFLOAT3 vel = m_playerVel;
        const float eps = 0.001f;
        auto sweepAxis = [&](int axis){
            float* pComp = (axis==0? &pos.x : (axis==1? &pos.y : &pos.z));
            float  vComp = (axis==0? vel.x : (axis==1? vel.y : vel.z));
            if(vComp == 0.0f) return;
            float newP = *pComp + vComp * dt;
            // Compute AABB after move along this axis
            float minX = (axis==0? newP - m_playerHalf.x : pos.x - m_playerHalf.x);
            float maxX = (axis==0? newP + m_playerHalf.x : pos.x + m_playerHalf.x);
            float minY = (axis==1? newP - m_playerHalf.y : pos.y - m_playerHalf.y);
            float maxY = (axis==1? newP + m_playerHalf.y : pos.y + m_playerHalf.y);
            float minZ = (axis==2? newP - m_playerHalf.z : pos.z - m_playerHalf.z);
            float maxZ = (axis==2? newP + m_playerHalf.z : pos.z + m_playerHalf.z);
            int x0 = (int)std::floor(minX); int x1 = (int)std::floor(maxX);
            int y0 = (int)std::floor(minY); int y1 = (int)std::floor(maxY);
            int z0 = (int)std::floor(minZ); int z1 = (int)std::floor(maxZ);
            bool collided = false;
            for(int z=z0; z<=z1; ++z){
                for(int y=y0; y<=y1; ++y){
                    for(int x=x0; x<=x1; ++x){
                        if(IsSolid(x,y,z)){
                            collided = true;
                            if(axis==0){
                                if(vComp > 0) newP = (float)x - m_playerHalf.x - eps;
                                else          newP = (float)(x + 1) + m_playerHalf.x + eps;
                                vel.x = 0.0f;
                            } else if(axis==1){
                                if(vComp > 0) newP = (float)y - m_playerHalf.y - eps;
                                else          newP = (float)(y + 1) + m_playerHalf.y + eps, m_onGround = true;
                                vel.y = 0.0f;
                            } else {
                                if(vComp > 0) newP = (float)z - m_playerHalf.z - eps;
                                else          newP = (float)(z + 1) + m_playerHalf.z + eps;
                                vel.z = 0.0f;
                            }
                            // Recompute AABB range because we snapped; update bounds to avoid repeated overlaps
                            minX = (axis==0? newP - m_playerHalf.x : pos.x - m_playerHalf.x);
                            maxX = (axis==0? newP + m_playerHalf.x : pos.x + m_playerHalf.x);
                            minY = (axis==1? newP - m_playerHalf.y : pos.y - m_playerHalf.y);
                            maxY = (axis==1? newP + m_playerHalf.y : pos.y + m_playerHalf.y);
                            minZ = (axis==2? newP - m_playerHalf.z : pos.z - m_playerHalf.z);
                            maxZ = (axis==2? newP + m_playerHalf.z : pos.z + m_playerHalf.z);
                            x0 = (int)std::floor(minX); x1 = (int)std::floor(maxX);
                            y0 = (int)std::floor(minY); y1 = (int)std::floor(maxY);
                            z0 = (int)std::floor(minZ); z1 = (int)std::floor(maxZ);
                        }
                    }
                }
            }
            *pComp = newP;
        };
        // Sweep Y then X then Z
        sweepAxis(1);
        sweepAxis(0);
        sweepAxis(2);
        m_playerPos = pos;
        m_playerVel = vel;
    }
};

int WINAPI wWinMain(HINSTANCE hInst, HINSTANCE, PWSTR, int) {
    VoxelApp app(hInst);
    if(!app.Initialize()) return -1;
    int code = app.Run();
    return code;
}
