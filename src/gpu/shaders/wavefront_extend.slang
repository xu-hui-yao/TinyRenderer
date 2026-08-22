// Stage 2, wavefront kernel 2/5: extend (closest-hit BVH traversal).
//
// Dispatched indirectly with groupCount = ceil(counters[0] / GROUP_SIZE)
// (see indirect_args[0..2], set by wavefront_raygen for bounce 0 and by
// wavefront_advance for every subsequent bounce). Each thread pops one
// pixel from the current bounce's active queue, traces its ray against the
// flattened BVH, and writes the (slot-indexed) hit result for
// wavefront_shade to consume. This is deliberately the ONLY thing this
// kernel does - separating "what did the ray hit" from "what do we do
// about it" is the core of the wavefront architecture: it keeps this
// kernel's register footprint small and divergence-free (every thread runs
// the identical BVH-traversal loop, regardless of what material/light it
// eventually hits), which is exactly what a megakernel's fused loop can't
// guarantee once BSDF-specific branches are interleaved with traversal.
#include "common/wavefront_common.slang"

[shader("compute")]
[numthreads(256, 1, 1)]
void main(uint3 tid: SV_DispatchThreadID)
{
    uint slot = tid.x;
    if (slot >= counters[0]) return;

    uint pixel = read_active(slot);

    float4 ro4 = path_ray_o[pixel];
    float4 rd4 = path_ray_d[pixel];

    HitResult hit = trace_ray(ro4.xyz, rd4.xyz, ro4.w, rd4.w);

    hit_tri[slot] = hit.hit ? hit.tri : INVALID_INDEX;
    hit_tuv[slot] = float4(hit.t, hit.u, hit.v, 0.0);
}
