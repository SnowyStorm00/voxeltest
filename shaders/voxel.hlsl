cbuffer VSConstants : register(b0)
{
    float4x4 gViewProj;
    // gJitter unused (TAA removed); keep for constant buffer layout compatibility
    float2 gJitter;
    float2 gPad;
    float3 gCameraPos;
    float  gFogStart;
    float  gFogEnd;
    float3 gFogColor;
    float3 gSunDir; // normalized sun direction
    float  gPad2;
};

struct VSIn { float3 pos:POSITION; float3 normal:NORMAL; float4 color:COLOR; };
struct PSIn { float4 pos:SV_POSITION; float3 normal:NORMAL; float4 color:COLOR; float3 worldPos:TEXCOORD0; };

PSIn VSMain(VSIn vin)
{
    PSIn o;
    o.pos = mul(float4(vin.pos,1), gViewProj);
    // TAA removed: no jitter offset
    o.normal = vin.normal;
    o.color = vin.color;
    o.worldPos = vin.pos;
    return o;
}

float4 PSMain(PSIn pin) : SV_Target
{
    float3 N = normalize(pin.normal);
    float3 L = normalize(gSunDir);
    float NdotL = saturate(dot(N, L));
    // Software shadow baked in alpha channel (0..1)
    float shadow = pin.color.a;
    // Hemi ambient (simple bounce): up vs down
    float hemi = 0.5 + 0.5 * saturate(N.y);
    float3 skyAmb = gFogColor * 0.35;
    float3 grdAmb = float3(0.05, 0.04, 0.03);
    float3 ambient = lerp(grdAmb, skyAmb, hemi);
    float3 base = pin.color.rgb;
    float3 direct = (0.85 * NdotL) * shadow;
    float3 lit = base * (ambient + direct);
    // Distance-based fog
    float dist = distance(pin.worldPos, gCameraPos);
    float denom = max(gFogEnd - gFogStart, 1e-3);
    float t = saturate((dist - gFogStart) / denom);
    float3 col = lerp(lit, gFogColor, t);
    return float4(col, 1);
}
