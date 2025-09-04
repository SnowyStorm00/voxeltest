cbuffer VSConstants : register(b1)
{
    float gScreenW;
    float gScreenH;
    float2 pad;
}

struct VSIn { float2 pos:POSITION; float4 color:COLOR; };
struct PSIn { float4 pos:SV_POSITION; float4 color:COLOR; };

PSIn VSMain(VSIn vin)
{
    PSIn o;
    // Convert pixel to NDC
    float2 ndc = float2((vin.pos.x / gScreenW) * 2.0 - 1.0, (vin.pos.y / gScreenH) * -2.0 + 1.0);
    o.pos = float4(ndc, 0.0, 1.0);
    o.color = vin.color;
    return o;
}

float4 PSMain(PSIn pin) : SV_Target
{
    return pin.color;
}
