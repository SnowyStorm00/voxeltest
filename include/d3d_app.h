#pragma once
#ifndef NOMINMAX
#define NOMINMAX
#endif
#include <windows.h>
#include <d3d11.h>
#include <wrl/client.h>
#include <string>

class D3DApp {
public:
    D3DApp(HINSTANCE hInstance, int width, int height, const std::wstring& title);
    ~D3DApp();

    virtual bool Initialize();
    int Run();

    HWND GetHwnd() const { return m_hWnd; }
    ID3D11Device* GetDevice() const { return m_device.Get(); }
    ID3D11DeviceContext* GetContext() const { return m_context.Get(); }
    IDXGISwapChain* GetSwapChain() const { return m_swapChain.Get(); }
    ID3D11RenderTargetView* GetRTV() const { return m_rtv.Get(); }
    ID3D11DepthStencilView* GetDSV() const { return m_dsv.Get(); }

protected:
    virtual void Update(float dt) {}
    virtual void Render() {}
    virtual void OnResize(int width, int height);
    virtual LRESULT MsgProc(HWND hwnd, UINT msg, WPARAM wParam, LPARAM lParam);

private:
    bool InitWindow();
    bool InitD3D();

    static LRESULT CALLBACK StaticWndProc(HWND hwnd, UINT msg, WPARAM wParam, LPARAM lParam);

protected:
    HINSTANCE m_hInstance;
    int m_width;
    int m_height;
    std::wstring m_title;

    HWND m_hWnd = nullptr;
    Microsoft::WRL::ComPtr<ID3D11Device> m_device;
    Microsoft::WRL::ComPtr<ID3D11DeviceContext> m_context;
    Microsoft::WRL::ComPtr<IDXGISwapChain> m_swapChain;
    Microsoft::WRL::ComPtr<ID3D11RenderTargetView> m_rtv;
    Microsoft::WRL::ComPtr<ID3D11Texture2D> m_depth;
    Microsoft::WRL::ComPtr<ID3D11DepthStencilView> m_dsv;
};
