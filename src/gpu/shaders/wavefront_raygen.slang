// Stage 2, wavefront kernel 1/5: ray generation.
//
// Dispatched once per sample-per-pixel (spp) with a FIXED size of
// ceil(width*height / GROUP_SIZE) groups (no indirect dispatch needed here -
// every pixel is always active at the start of a new spp). For each pixel:
// generates the primary camera ray (identical math to path_trace.slang's
// megakernel), initializes its persistent path state, seeds it into the
// bounce-0 active queue (always active_indices_a - see file header of
// wavefront_common.slang), and increments its running sample count. Thread 0
// additionally (re)initializes the counters/indirect-args buffers for the
// upcoming bounce-0 extend/shade dispatch.
#include "common/wavefront_common.slang"

[shader("compute")]
[numthreads(256, 1, 1)]
void main(uint3 tid: SV_DispatchThreadID)
{
    uint pixel = tid.x;
    if (pixel >= pc.total_pixels) return;

    uint px = pixel % camera.width;
    uint py = pixel / camera.width;

    Rng rng = rng_init(pixel, pc.spp_index, 0u);

    float2 jitter = rng_next2(rng);
    float fx = float(px) + jitter.x;
    float fy = float(py) + jitter.y;
    float inv_w = 1.0 / float(camera.width);
    float inv_h = 1.0 / float(camera.height);

    float3 near_p = transform_point(camera.s2c0, camera.s2c1, camera.s2c2, camera.s2c3, float3(fx * inv_w, fy * inv_h, 0.0));
    float3 d = normalize(near_p);
    float inv_z = 1.0 / d.z;

    float3 ray_o = transform_point(camera.c2w0, camera.c2w1, camera.c2w2, camera.c2w3, float3(0.0, 0.0, 0.0));
    float3 ray_d = transform_vector(camera.c2w0, camera.c2w1, camera.c2w2, camera.c2w3, d);
    float ray_tmin = camera.near_clip * inv_z;
    float ray_tmax = camera.far_clip * inv_z;

    path_ray_o[pixel] = float4(ray_o, ray_tmin);
    path_ray_d[pixel] = float4(ray_d, ray_tmax);
    path_throughput[pixel] = float4(1.0, 1.0, 1.0, 1.0); // .xyz=throughput=1, .w=eta=1
    path_rng[pixel] = rng.state;
    path_prev_bsdf[pixel] = float2(1.0, 1.0); // pdf=1, delta=true (matches megakernel's initial state)
    path_prev_p[pixel] = float4(ray_o, 0.0);

    active_indices_a[pixel] = pixel;

    // Running sample count lives in accum_buffer[pixel].w; incremented
    // exactly once per pixel per spp, decoupled from whether this sample
    // ends up contributing any radiance (matches the megakernel's `+1.0`).
    // Safe as a plain (non-atomic) read-modify-write: raygen has exactly
    // one thread per pixel, and runs-and-completes (via a pipeline barrier)
    // strictly before any shade/shadow dispatch of THIS spp touches
    // accum_buffer, and strictly after all of the PREVIOUS spp's.
    float4 acc = accum_buffer[pixel];
    acc.w += 1.0;
    accum_buffer[pixel] = acc;

    if (pixel == 0u) {
        counters[0] = pc.total_pixels; // bounce-0 extend/shade will process all pixels
        counters[1] = 0u;
        counters[2] = 0u;
        uint groups = (pc.total_pixels + GROUP_SIZE - 1u) / GROUP_SIZE;
        indirect_args[0] = groups; indirect_args[1] = 1u; indirect_args[2] = 1u;
        indirect_args[3] = 0u; indirect_args[4] = 1u; indirect_args[5] = 1u;
    }
}
