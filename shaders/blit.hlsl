Texture2D gTex : register(t0);
SamplerState gSamp : register(s0);

cbuffer VSConst : register(b0)
{
    float2 gScaleBias; // not used
    float2 gPad;
}

struct VSIn { float2 pos:POSITION; float2 uv:TEXCOORD0; };
struct PSIn { float4 pos:SV_POSITION; float2 uv:TEXCOORD0; };

PSIn VSMain(VSIn vin){ PSIn o; o.pos=float4(vin.pos,0,1); o.uv=vin.uv; return o; }
float4 PSMain(PSIn i):SV_Target { return gTex.Sample(gSamp, i.uv); }
