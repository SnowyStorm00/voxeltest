#include "d3d_app.h"
#include <stdexcept>
#include <cassert>

using Microsoft::WRL::ComPtr;

D3DApp::D3DApp(HINSTANCE hInstance, int width, int height, const std::wstring& title)
    : m_hInstance(hInstance), m_width(width), m_height(height), m_title(title) {}

D3DApp::~D3DApp() {}

bool D3DApp::Initialize() {
    if(!InitWindow()) return false;
    if(!InitD3D()) return false;
    return true;
}

int D3DApp::Run() {
    MSG msg = {0};
    LARGE_INTEGER freq, prev, now;
    QueryPerformanceFrequency(&freq);
    QueryPerformanceCounter(&prev);

    while(msg.message != WM_QUIT) {
        if(PeekMessage(&msg, nullptr, 0, 0, PM_REMOVE)) {
            TranslateMessage(&msg);
            DispatchMessage(&msg);
        } else {
            QueryPerformanceCounter(&now);
            float dt = float(now.QuadPart - prev.QuadPart) / float(freq.QuadPart);
            prev = now;
            Update(dt);
            Render();
        }
    }
    return (int)msg.wParam;
}

bool D3DApp::InitWindow() {
    WNDCLASSEXW wc = { sizeof(WNDCLASSEXW) };
    wc.style = CS_HREDRAW | CS_VREDRAW;
    wc.lpfnWndProc = StaticWndProc;
    wc.hInstance = m_hInstance;
    wc.hCursor = LoadCursor(nullptr, IDC_ARROW);
    wc.lpszClassName = L"VoxelWindowClass";
    RegisterClassExW(&wc);

    RECT rc = {0, 0, m_width, m_height};
    AdjustWindowRect(&rc, WS_OVERLAPPEDWINDOW, FALSE);

    m_hWnd = CreateWindowExW(0, wc.lpszClassName, m_title.c_str(),
        WS_OVERLAPPEDWINDOW, CW_USEDEFAULT, CW_USEDEFAULT,
        rc.right - rc.left, rc.bottom - rc.top, nullptr, nullptr, m_hInstance, this);

    ShowWindow(m_hWnd, SW_SHOW);
    return m_hWnd != nullptr;
}

bool D3DApp::InitD3D() {
    UINT createDeviceFlags = 0;
#ifdef _DEBUG
    createDeviceFlags |= D3D11_CREATE_DEVICE_DEBUG;
#endif

    DXGI_SWAP_CHAIN_DESC sd{};
    sd.BufferCount = 1;
    sd.BufferDesc.Width = m_width;
    sd.BufferDesc.Height = m_height;
    sd.BufferDesc.Format = DXGI_FORMAT_R8G8B8A8_UNORM;
    sd.BufferUsage = DXGI_USAGE_RENDER_TARGET_OUTPUT;
    sd.OutputWindow = m_hWnd;
    sd.SampleDesc.Count = 1;
    sd.Windowed = TRUE;

    D3D_FEATURE_LEVEL featureLevels[] = { D3D_FEATURE_LEVEL_11_0 };
    D3D_FEATURE_LEVEL createdFL;

    HRESULT hr = D3D11CreateDeviceAndSwapChain(
        nullptr, D3D_DRIVER_TYPE_HARDWARE, nullptr, createDeviceFlags,
        featureLevels, 1, D3D11_SDK_VERSION, &sd, m_swapChain.GetAddressOf(),
        m_device.GetAddressOf(), &createdFL, m_context.GetAddressOf());
    if(FAILED(hr)) return false;

    ComPtr<ID3D11Texture2D> backBuffer;
    hr = m_swapChain->GetBuffer(0, IID_PPV_ARGS(&backBuffer));
    if(FAILED(hr)) return false;
    hr = m_device->CreateRenderTargetView(backBuffer.Get(), nullptr, m_rtv.GetAddressOf());
    if(FAILED(hr)) return false;

    D3D11_TEXTURE2D_DESC depthDesc{};
    depthDesc.Width = m_width;
    depthDesc.Height = m_height;
    depthDesc.MipLevels = 1;
    depthDesc.ArraySize = 1;
    depthDesc.Format = DXGI_FORMAT_D24_UNORM_S8_UINT;
    depthDesc.SampleDesc.Count = 1;
    depthDesc.Usage = D3D11_USAGE_DEFAULT;
    depthDesc.BindFlags = D3D11_BIND_DEPTH_STENCIL;
    hr = m_device->CreateTexture2D(&depthDesc, nullptr, m_depth.GetAddressOf());
    if(FAILED(hr)) return false;
    hr = m_device->CreateDepthStencilView(m_depth.Get(), nullptr, m_dsv.GetAddressOf());
    if(FAILED(hr)) return false;

    m_context->OMSetRenderTargets(1, m_rtv.GetAddressOf(), m_dsv.Get());

    D3D11_VIEWPORT vp{};
    vp.Width = (FLOAT)m_width;
    vp.Height = (FLOAT)m_height;
    vp.MinDepth = 0.0f;
    vp.MaxDepth = 1.0f;
    m_context->RSSetViewports(1, &vp);

    return true;
}

void D3DApp::OnResize(int width, int height) {
    m_width = width; m_height = height;
    if(!m_swapChain) return;

    m_context->OMSetRenderTargets(0, nullptr, nullptr);
    m_rtv.Reset();
    m_dsv.Reset();
    m_depth.Reset();

    m_swapChain->ResizeBuffers(0, width, height, DXGI_FORMAT_UNKNOWN, 0);

    Microsoft::WRL::ComPtr<ID3D11Texture2D> backBuffer;
    m_swapChain->GetBuffer(0, IID_PPV_ARGS(&backBuffer));
    m_device->CreateRenderTargetView(backBuffer.Get(), nullptr, m_rtv.GetAddressOf());

    D3D11_TEXTURE2D_DESC depthDesc{};
    depthDesc.Width = width;
    depthDesc.Height = height;
    depthDesc.MipLevels = 1;
    depthDesc.ArraySize = 1;
    depthDesc.Format = DXGI_FORMAT_D24_UNORM_S8_UINT;
    depthDesc.SampleDesc.Count = 1;
    depthDesc.Usage = D3D11_USAGE_DEFAULT;
    depthDesc.BindFlags = D3D11_BIND_DEPTH_STENCIL;
    m_device->CreateTexture2D(&depthDesc, nullptr, m_depth.GetAddressOf());
    m_device->CreateDepthStencilView(m_depth.Get(), nullptr, m_dsv.GetAddressOf());

    m_context->OMSetRenderTargets(1, m_rtv.GetAddressOf(), m_dsv.Get());

    D3D11_VIEWPORT vp{};
    vp.Width = (FLOAT)width;
    vp.Height = (FLOAT)height;
    vp.MinDepth = 0.0f;
    vp.MaxDepth = 1.0f;
    m_context->RSSetViewports(1, &vp);
}

LRESULT D3DApp::MsgProc(HWND hwnd, UINT msg, WPARAM wParam, LPARAM lParam) {
    switch(msg) {
        case WM_SIZE:
            OnResize(LOWORD(lParam), HIWORD(lParam));
            return 0;
        default:
            return DefWindowProc(hwnd, msg, wParam, lParam);
    }
}

LRESULT CALLBACK D3DApp::StaticWndProc(HWND hwnd, UINT msg, WPARAM wParam, LPARAM lParam) {
    if(msg == WM_NCCREATE) {
        auto cs = reinterpret_cast<CREATESTRUCT*>(lParam);
        auto app = reinterpret_cast<D3DApp*>(cs->lpCreateParams);
        SetWindowLongPtr(hwnd, GWLP_USERDATA, (LONG_PTR)app);
        app->m_hWnd = hwnd;
    }
    auto app = reinterpret_cast<D3DApp*>(GetWindowLongPtr(hwnd, GWLP_USERDATA));
    if(app) return app->MsgProc(hwnd, msg, wParam, lParam);
    return DefWindowProc(hwnd, msg, wParam, lParam);
}
