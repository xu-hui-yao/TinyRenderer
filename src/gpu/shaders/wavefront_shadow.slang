// Stage 2, wavefront kernel 4/5: shadow-ray (NEE visibility) tracing.
//
// Dispatched indirectly with groupCount = ceil(counters[2] / GROUP_SIZE),
// where counters[2] is the shadow-ray count wavefront_shade just built for
// THIS bounce (snapshotted into indirect_args[3..5] by
// wavefront_prepare_shadow_indirect, which runs between shade and this
// kernel WITHOUT resetting counters[2] itself - see that file's header for
// why the ordering here matters). Each thread traces exactly one queued
// shadow ray and, if unoccluded, adds its pre-computed NEE contribution
// into accum_buffer. Kept in its own dispatch (rather than traced inline
// inside wavefront_shade) so shadow-ray tracing - a second, cheaper
// BVH-traversal workload with a different divergence profile (any-hit
// early-exit vs closest-hit) - doesn't share a kernel invocation with the
// much heavier material evaluation code in shade.
#include "common/wavefront_common.slang"

[shader("compute")]
[numthreads(256, 1, 1)]
void main(uint3 tid: SV_DispatchThreadID)
{
    uint slot = tid.x;
    if (slot >= counters[2]) return;

    float4 o4 = shadow_o[slot];
    float4 d4 = shadow_d[slot];
    uint pixel = shadow_pixel[slot];
    float4 contrib = shadow_contrib[slot];

    bool occluded = trace_shadow_ray(o4.xyz, d4.xyz, 0.0, o4.w);
    if (!occluded) {
        // Safe as a plain += : wavefront_shade only ever queues at most one
        // shadow ray per active pixel per bounce, so `pixel` here is unique
        // across all threads of this dispatch, and (per wavefront_common.
        // slang's file header) can never coincide with a pixel that
        // wavefront_shade itself already wrote to accum_buffer this bounce.
        float4 acc = accum_buffer[pixel];
        acc.xyz += contrib.xyz;
        accum_buffer[pixel] = acc;
    }
}
