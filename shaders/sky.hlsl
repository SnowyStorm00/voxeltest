cbuffer SkyCB : register(b0)
{
    float4x4 gInvViewProj; // for ray direction
    float3   gSunDir;      // normalized
    float    gSunSize;     // radians, apparent radius
    float3   gSkyHorizon;  // horizon color
    float    gPad0;
    float3   gSkyZenith;   // zenith color
    float    gPad1;
}

struct VSIn { float2 pos:POSITION; };
struct PSIn { float4 pos:SV_POSITION; float3 ray:TEXCOORD0; };

PSIn VSMain(VSIn vin)
{
    // Fullscreen tri in clip space
    float2 p = vin.pos; // expects {-1,-1}, {3,-1}, {-1,3}
    PSIn o; o.pos = float4(p, 0, 1);
    // reconstruct ray dir by unprojecting NDC
    float2 ndc = p;
    float4 h0 = mul(float4(ndc, 0, 1), gInvViewProj); // near plane
    float4 h1 = mul(float4(ndc, 1, 1), gInvViewProj); // far plane
    float3 w0 = h0.xyz / max(h0.w, 1e-5);
    float3 w1 = h1.xyz / max(h1.w, 1e-5);
    o.ray = normalize(w1 - w0);
    return o;
}

float smoothstepf(float a, float b, float x){
    float t = saturate((x-a)/(b-a)); return t*t*(3.0-2.0*t);
}

float4 PSMain(PSIn i) : SV_Target
{
    // Simple gradient sky + sun disk
    float3 up = float3(0,1,0);
    float t = saturate(0.5 + 0.5*dot(i.ray, up));
    float3 sky = lerp(gSkyHorizon, gSkyZenith, t);

    // Sun disk intensity
    float cosAng = dot(normalize(i.ray), normalize(gSunDir));
    float ang = acos(saturate(cosAng));
    float core = 1.0 - smoothstepf(gSunSize*0.6, gSunSize, ang);
    float glow = 1.0 - smoothstepf(gSunSize*1.5, gSunSize*4.0, ang);
    float3 sunCol = float3(1.0, 0.95, 0.85);
    float3 col = sky + sunCol * (core*6.0 + glow*0.3);
    return float4(col, 1.0);
}
