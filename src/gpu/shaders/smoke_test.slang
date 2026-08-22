// GPU compute smoke test (Stage 1, checkpoint 1).
//
// This shader does not touch the renderer's scene data yet - its only job is
// to prove that the full toolchain (Slang -> SPIR-V -> MoltenVK -> Metal on
// this machine, or -> native Vulkan on Windows/Linux) works end-to-end by
// procedurally painting a test pattern into a buffer that the host then
// reads back and saves as a PNG.

struct PushConstants {
    uint width;
    uint height;
};

[[vk::push_constant]]
ConstantBuffer<PushConstants> pc;

[[vk::binding(0, 0)]]
RWStructuredBuffer<float4> output;

[shader("compute")]
[numthreads(16, 16, 1)]
void main(uint3 tid: SV_DispatchThreadID)
{
    if (tid.x >= pc.width || tid.y >= pc.height)
        return;

    float2 uv = float2(tid.x, tid.y) / float2(float(pc.width), float(pc.height));

    // A simple radial "lens flare" style pattern: colorful UV gradient with
    // a soft ring, chosen only because it makes a toolchain regression
    // (e.g. a broken dispatch, wrong buffer stride, or NaN from a bad
    // push-constant read) immediately, visually obvious in the output PNG.
    float2 centered = uv - float2(0.5, 0.5);
    float dist       = length(centered);
    float ring       = smoothstep(0.32, 0.30, dist) - smoothstep(0.30, 0.28, dist);

    float3 base_color = float3(uv.x, uv.y, 1.0 - 0.5 * (uv.x + uv.y));
    float3 ring_color  = float3(1.0, 0.85, 0.2);
    float3 color       = lerp(base_color, ring_color, saturate(ring * 4.0));

    uint idx     = tid.y * pc.width + tid.x;
    output[idx] = float4(color, 1.0);
}
