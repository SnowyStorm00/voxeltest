cbuffer VSConstants : register(b0)
{
    float4x4 gViewProj;
};

struct VSIn { float3 pos:POSITION; float3 normal:NORMAL; float4 color:COLOR; };
struct PSIn { float4 pos:SV_POSITION; float3 normal:NORMAL; float4 color:COLOR; };

PSIn VSMain(VSIn vin)
{
    PSIn o;
    o.pos = mul(float4(vin.pos,1), gViewProj);
    o.normal = vin.normal;
    o.color = vin.color;
    return o;
}

float4 PSMain(PSIn pin) : SV_Target
{
    float3 L = normalize(float3(0.5, 1, 0.2));
    float NdotL = saturate(dot(normalize(pin.normal), L));
    float3 base = pin.color.rgb;
    float3 col = base * (0.2 + 0.8*NdotL);
    return float4(col, 1);
}
