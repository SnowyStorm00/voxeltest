Texture2D gCurr : register(t0);
Texture2D gPrev : register(t1);
Texture2D gVel  : register(t2);
SamplerState gSamp : register(s0);

cbuffer TAAConst : register(b0)
{
    float2 gInvResolution; // 1/width, 1/height
    float  gAlpha;         // blend factor for history (e.g., 0.1..0.2)
    float  gPad;
    float2 gDeltaUV;       // previousJitterUV - currentJitterUV
    float2 gPad2;
}

struct VSIn { float2 pos:POSITION; float2 uv:TEXCOORD0; };
struct PSIn { float4 pos:SV_POSITION; float2 uv:TEXCOORD0; };

PSIn VSMain(VSIn vin){ PSIn o; o.pos=float4(vin.pos,0,1); o.uv=vin.uv; return o; }

float Luma(float3 c){ return dot(c, float3(0.299, 0.587, 0.114)); }

float3 NeighborhoodMin(Texture2D tex, float2 uv){
    float3 m = tex.SampleLevel(gSamp, uv, 0).rgb;
    [unroll] for(int y=-1;y<=1;++y){ [unroll] for(int x=-1;x<=1;++x){
        float2 o = float2(x,y)*gInvResolution;
        m = min(m, tex.SampleLevel(gSamp, uv+o, 0).rgb);
    }}
    return m;
}
float3 NeighborhoodMax(Texture2D tex, float2 uv){
    float3 M = tex.SampleLevel(gSamp, uv, 0).rgb;
    [unroll] for(int y=-1;y<=1;++y){ [unroll] for(int x=-1;x<=1;++x){
        float2 o = float2(x,y)*gInvResolution;
        M = max(M, tex.SampleLevel(gSamp, uv+o, 0).rgb);
    }}
    return M;
}

// Catmull-Rom helpers
float w0(float a){ return ((-a + 2.0*a*a - a*a*a) * 0.5); }
float w1(float a){ return ((2.0 - 5.0*a + 4.0*a*a - a*a*a) * 0.5); }
float w2(float a){ return ((-a + 4.0*a*a - 3.0*a*a*a) * 0.5); }
float w3(float a){ return ((-a*a + a*a*a) * 0.5); }

float3 SampleCatmullRom(Texture2D tex, float2 uv){
    float2 texSize = 1.0 / gInvResolution;
    float2 pos = uv * texSize - 0.5;
    float2 pf = frac(pos);
    int2 pi = (int2)floor(pos);
    // base pixel center
    float2 baseUV = (float2(pi) + 0.5) / texSize;
    // Sample 4x4 neighborhood
    float3 c00 = tex.SampleLevel(gSamp, baseUV + gInvResolution*float2(-1,-1), 0).rgb;
    float3 c10 = tex.SampleLevel(gSamp, baseUV + gInvResolution*float2( 0,-1), 0).rgb;
    float3 c20 = tex.SampleLevel(gSamp, baseUV + gInvResolution*float2( 1,-1), 0).rgb;
    float3 c30 = tex.SampleLevel(gSamp, baseUV + gInvResolution*float2( 2,-1), 0).rgb;
    float3 c01 = tex.SampleLevel(gSamp, baseUV + gInvResolution*float2(-1, 0), 0).rgb;
    float3 c11 = tex.SampleLevel(gSamp, baseUV + gInvResolution*float2( 0, 0), 0).rgb;
    float3 c21 = tex.SampleLevel(gSamp, baseUV + gInvResolution*float2( 1, 0), 0).rgb;
    float3 c31 = tex.SampleLevel(gSamp, baseUV + gInvResolution*float2( 2, 0), 0).rgb;
    float3 c02 = tex.SampleLevel(gSamp, baseUV + gInvResolution*float2(-1, 1), 0).rgb;
    float3 c12 = tex.SampleLevel(gSamp, baseUV + gInvResolution*float2( 0, 1), 0).rgb;
    float3 c22 = tex.SampleLevel(gSamp, baseUV + gInvResolution*float2( 1, 1), 0).rgb;
    float3 c32 = tex.SampleLevel(gSamp, baseUV + gInvResolution*float2( 2, 1), 0).rgb;
    float3 c03 = tex.SampleLevel(gSamp, baseUV + gInvResolution*float2(-1, 2), 0).rgb;
    float3 c13 = tex.SampleLevel(gSamp, baseUV + gInvResolution*float2( 0, 2), 0).rgb;
    float3 c23 = tex.SampleLevel(gSamp, baseUV + gInvResolution*float2( 1, 2), 0).rgb;
    float3 c33 = tex.SampleLevel(gSamp, baseUV + gInvResolution*float2( 2, 2), 0).rgb;
    float3 row0 = c00*w0(pf.x) + c10*w1(pf.x) + c20*w2(pf.x) + c30*w3(pf.x);
    float3 row1 = c01*w0(pf.x) + c11*w1(pf.x) + c21*w2(pf.x) + c31*w3(pf.x);
    float3 row2 = c02*w0(pf.x) + c12*w1(pf.x) + c22*w2(pf.x) + c32*w3(pf.x);
    float3 row3 = c03*w0(pf.x) + c13*w1(pf.x) + c23*w2(pf.x) + c33*w3(pf.x);
    float3 col  = row0*w0(pf.y) + row1*w1(pf.y) + row2*w2(pf.y) + row3*w3(pf.y);
    return saturate(col);
}

float4 PSMain(PSIn i) : SV_Target
{
    float2 uv = i.uv;
    float3 curr = gCurr.Sample(gSamp, uv).rgb;
    // Reproject previous using jitter delta
    float2 uvPrev = uv + gDeltaUV;
    float3 prev = SampleCatmullRom(gPrev, uvPrev);

    float3 nmin = NeighborhoodMin(gCurr, uv);
    float3 nmax = NeighborhoodMax(gCurr, uv);
    prev = clamp(prev, nmin, nmax);
    // Luma clamp for extra stability
    float lmin = Luma(nmin);
    float lmax = Luma(nmax);
    float lp = Luma(prev);
    float lpc = clamp(lp, lmin, lmax);
    float3 prevDir = (lp > 1e-5) ? (prev / lp) : float3(1,1,1);
    prev = prevDir * lpc;

    // Adaptive per-pixel alpha based on local contrast
    float contrast = lmax - lmin;
    float alphaLocal = lerp(gAlpha, gAlpha*0.4, saturate(contrast*4.0));
    float3 outc = lerp(curr, prev, alphaLocal);
    return float4(outc, 1);
}
