cbuffer VSConstants : register(b0) {
    float4x4 gWorldViewProj;
    float3   gColor;
    float    pad0;
    float3   gSunDir;
    float    pad1;
    float4x4 gLightViewProj;
}

Texture2D<float> gShadowMap : register(t0);
SamplerComparisonState gShadowSamp : register(s0);

struct VSIn { float3 pos : POSITION; float3 normal : NORMAL; };
struct VSOut { float4 pos : SV_POSITION; float3 n : NORMAL; };

VSOut VSMain(VSIn i){
    VSOut o; o.pos = mul(float4(i.pos,1.0), gWorldViewProj); o.n = i.normal; return o;
}

float4 PSMain(VSOut i) : SV_TARGET {
    float3 N = normalize(i.n);
    float3 L = normalize(gSunDir);
    float NdotL = max(dot(N, L), 0.0);
    // Assume model in world space near camera; approximate pos from viewProj not available here; use N.L only and ambient
    float hemi = 0.5 + 0.5 * saturate(N.y);
    float3 skyAmb = float3(0.35,0.45,0.6) * 0.35;
    float3 grdAmb = float3(0.05,0.04,0.03);
    float3 ambient = lerp(grdAmb, skyAmb, hemi);
    float3 c = gColor * (ambient + 0.85*NdotL);
    return float4(c, 1.0);
}
