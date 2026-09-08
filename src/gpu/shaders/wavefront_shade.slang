// Stage 2, wavefront kernel 3/5: shade.
//
// Dispatched with the SAME indirect args (and therefore the same slot range
// tid.x < counters[0]) as the wavefront_extend dispatch that just ran, so
// hit_tri[slot]/hit_tuv[slot] line up with read_active(slot)'s path state.
// This is a line-by-line port of path_trace.slang's per-bounce loop body
// (which is itself a port of src/integrators/path.cpp's li() loop; see that
// file's header for the full CPU-parity notes), restructured to:
//   - accumulate direct emission (with MIS against the PREVIOUS bounce's
//     BSDF sample, carried in path_prev_bsdf/path_prev_p) directly here,
//   - queue a shadow ray (rather than tracing it inline) for NEE, so
//     shadow-ray tracing can be its own dedicated dispatch/kernel, and
//   - if the path survives (didn't terminate, passed Russian roulette):
//     write its updated persistent state back to path_ray_o/d/throughput/
//     rng/prev_bsdf/prev_p (indexed by PIXEL) and atomically append it to
//     the NEXT bounce's active queue (indexed by a fresh SLOT).
#include "common/wavefront_common.slang"

[shader("compute")]
[numthreads(256, 1, 1)]
void main(uint3 tid: SV_DispatchThreadID)
{
    uint slot = tid.x;
    if (slot >= counters[0]) return;

    uint pixel = read_active(slot);

    uint tri = hit_tri[slot];
    float4 tuv = hit_tuv[slot];
    bool hit_geometry = tri != INVALID_INDEX;
    // Mirrors Scene::ray_intersect (scene.cpp): a ray that hits no geometry
    // still counts as an "intersection" when the scene has an environment
    // map - the CPU fabricates a hit at t = M_MAX_FLOAT with
    // mesh_id = M_INVALID_INDEX so the direct-emission block below can pick
    // the environment emitter up. Without this, the envmap branch of that
    // block is unreachable and the sky never contributes.
    bool is_intersect = hit_geometry || pc.environment_light_id != INVALID_INDEX;

    float4 ro4 = path_ray_o[pixel];
    float4 rd4 = path_ray_d[pixel];
    float3 ray_d = rd4.xyz;
    float4 th4 = path_throughput[pixel];
    float3 throughput = th4.xyz;
    float eta = th4.w;
    Rng rng; rng.state = path_rng[pixel];
    float2 prev_bsdf = path_prev_bsdf[pixel];
    float prev_bsdf_pdf = prev_bsdf.x;
    bool prev_bsdf_delta = prev_bsdf.y != 0.0;
    float3 prev_p = path_prev_p[pixel].xyz;

    uint hit_mesh_material = INVALID_INDEX;
    uint hit_mesh_light = INVALID_INDEX;
    float3 hit_p = float3(0, 0, 0), hit_n_geo = float3(0, 0, 1), hit_shading_n = float3(0, 0, 1);
    float3 hit_dp_du = float3(1, 0, 0), hit_dp_dv = float3(0, 1, 0);
    float2 hit_uv = float2(0, 0);
    float3 hit_wi_world = -ray_d;

    if (hit_geometry) {
        hit_mesh_material = triangle_material_id[tri];
        hit_mesh_light = triangle_light_id[tri];

        uint i0 = indices[tri * 3 + 0], i1 = indices[tri * 3 + 1], i2 = indices[tri * 3 + 2];
        float3 p0 = positions[i0].xyz, p1 = positions[i1].xyz, p2 = positions[i2].xyz;
        float3 bary = float3(1.0 - tuv.y - tuv.z, tuv.y, tuv.z);
        hit_p = p0 * bary.x + p1 * bary.y + p2 * bary.z;
        hit_uv = uvs[i0] * bary.x + uvs[i1] * bary.y + uvs[i2] * bary.z;

        hit_n_geo = normalize(cross(p1 - p0, p2 - p0));
        float3 n0 = normals[i0].xyz, n1 = normals[i1].xyz, n2 = normals[i2].xyz;
        float nlen2 = dot(n0, n0) + dot(n1, n1) + dot(n2, n2);
        hit_shading_n = nlen2 > 0.0 ? normalize(n0 * bary.x + n1 * bary.y + n2 * bary.z) : hit_n_geo;

        float2 uv0 = uvs[i0], uv1 = uvs[i1], uv2 = uvs[i2];
        float2 duv1 = uv1 - uv0, duv2 = uv2 - uv0;
        float3 edge1 = p1 - p0, edge2 = p2 - p0;
        float detuv = duv1.x * duv2.y - duv1.y * duv2.x;
        if (abs(detuv) < 1e-8) {
            hit_dp_du = edge1; hit_dp_dv = edge2;
        } else {
            float invd = 1.0 / detuv;
            hit_dp_du = (edge1 * duv2.y - edge2 * duv1.y) * invd;
            hit_dp_dv = (-edge1 * duv2.x + edge2 * duv1.x) * invd;
        }
    }

    bool hit_mesh_valid = is_intersect && hit_mesh_material != INVALID_INDEX;

    // ---------------------- Direct emission ----------------------
    // Safe as a plain += (no atomics): within this dispatch, each `pixel`
    // appears in exactly one thread (read_active(slot) is a permutation of
    // the active set), so no two threads ever touch the same
    // accum_buffer[pixel]. Across kernels within this SAME bounce, a pixel
    // that accumulates emission here (and therefore returns below, without
    // being queued for NEE) can never ALSO be touched by wavefront_shadow
    // this bounce - see wavefront_common.slang's file header for the full
    // argument.
    if (is_intersect && (!hit_mesh_valid || hit_mesh_light != INVALID_INDEX)) {
        float3 emission = float3(0, 0, 0);
        float em_pdf = 0.0;
        if (hit_mesh_light != INVALID_INDEX) {
            GPULightGPU L = lights[hit_mesh_light];
            bool front = dot(hit_wi_world, hit_shading_n) > 0.0;
            if (front) emission = tex_eval(L.radiance_tex, hit_uv);

            if (!prev_bsdf_delta) {
                float3 ds_d = normalize(hit_p - prev_p);
                float ds_dist = length(hit_p - prev_p);
                em_pdf = pdf_emitter_direction(hit_mesh_light, ds_d, hit_shading_n, ds_dist, pc.num_lights);
            }
        } else if (!hit_mesh_valid && pc.environment_light_id != INVALID_INDEX) {
            emission = eval_envmap(pc.environment_light_id, ray_d);
            if (!prev_bsdf_delta) {
                em_pdf = pdf_emitter_direction(pc.environment_light_id, ray_d, -ray_d, 1.0, pc.num_lights);
            }
        }
        float w = mis_weight(prev_bsdf_pdf, em_pdf);
        float4 acc = accum_buffer[pixel];
        acc.xyz += throughput * emission * w;
        accum_buffer[pixel] = acc;
    }

    bool active_next = (pc.depth + 1u < pc.max_depth) && is_intersect && hit_mesh_valid;
    if (!active_next) return; // path terminates: simply absent from the next active queue

    SurfaceHit sh;
    sh.uv = hit_uv;
    sh.shading_n = hit_shading_n;
    sh.shading_s = normalize(hit_dp_du - hit_shading_n * dot(hit_dp_du, hit_shading_n));
    if (length(sh.shading_s) < 1e-6) {
        float3 up = abs(hit_shading_n.z) < 0.999 ? float3(0, 0, 1) : float3(1, 0, 0);
        sh.shading_s = normalize(cross(up, hit_shading_n));
    }
    sh.shading_t = cross(hit_shading_n, sh.shading_s);
    sh.dp_du = hit_dp_du;
    sh.dp_dv = hit_dp_dv;
    sh.wi = frame_to_local(sh, hit_wi_world);

    // ---------------------- Emitter sampling (NEE) ----------------------
    bool active_em = bsdf_is_smooth(hit_mesh_material);
    DirSample ds;
    float3 em_weight = float3(0, 0, 0);
    float3 wo_world = float3(0, 0, 0);
    float3 wo_local = float3(0, 0, 0);

    if (active_em) {
        ds = sample_emitter_direction(hit_p, rng, pc.num_lights, em_weight);
        active_em = active_em && ds.pdf != 0.0;
        wo_world = ds.d;
        wo_local = frame_to_local(sh, wo_world);
    }

    float sample1 = rng_next(rng);
    float2 sample2 = rng_next2(rng);

    float3 bsdf_val = bsdf_eval(hit_mesh_material, sh, wo_local);
    float bsdf_pdf_v = bsdf_pdf(hit_mesh_material, sh, wo_local);
    float3 bsdf_weight;
    BSDFSample bs = bsdf_sample(hit_mesh_material, sh, sample1, sample2, bsdf_weight);

    // Queue a shadow ray instead of tracing it inline - wavefront_shadow
    // will trace it and (if unoccluded) add its contribution to
    // accum_buffer[pixel] in a separate, dedicated dispatch.
    if (active_em) {
        float mis_em = ds.delta ? 1.0 : mis_weight(ds.pdf, bsdf_pdf_v);
        float3 shadow_origin = offset_p(hit_p, hit_shading_n, wo_world);
        float3 to_target = ds.p - shadow_origin;
        float shadow_dist = length(to_target);
        float3 shadow_dir = to_target / max(shadow_dist, 1e-8);
        float3 contrib = throughput * bsdf_val * em_weight * mis_em;

        uint shadow_slot;
        InterlockedAdd(counters[2], 1u, shadow_slot);
        shadow_o[shadow_slot] = float4(shadow_origin, shadow_dist * (1.0 - 1e-4));
        shadow_d[shadow_slot] = float4(shadow_dir, 0.0);
        shadow_pixel[shadow_slot] = pixel;
        shadow_contrib[shadow_slot] = float4(contrib, 0.0);
    }

    // ---------------------- BSDF sampling / path update ----------------------
    float3 next_wo_world = frame_to_world(sh, bs.wo);
    float3 next_o = offset_p(hit_p, hit_shading_n, next_wo_world);

    throughput *= bsdf_weight;
    eta *= bs.eta;

    float throughput_max = max(max(throughput.x, throughput.y), throughput.z);
    float rr_prob = min(throughput_max * eta * eta, 0.95);
    bool rr_active = (pc.depth + 1u) >= pc.rr_depth;
    bool rr_continue = rng_next(rng) < rr_prob;
    if (rr_active) throughput /= rr_prob;
    bool survives = (!rr_active || rr_continue) && throughput_max != 0.0;

    if (survives) {
        path_ray_o[pixel] = float4(next_o, 0.0);
        path_ray_d[pixel] = float4(next_wo_world, 3.402823466e+30);
        path_throughput[pixel] = float4(throughput, eta);
        path_rng[pixel] = rng.state;
        path_prev_bsdf[pixel] = float2(bs.pdf, bs.delta ? 1.0 : 0.0);
        path_prev_p[pixel] = float4(hit_p, 0.0);

        uint next_slot;
        InterlockedAdd(counters[1], 1u, next_slot);
        write_next_active(next_slot, pixel);
    }
}
