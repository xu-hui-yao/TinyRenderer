// Stage 2: shared code for the wavefront path tracer.
//
// This is #include'd (not compiled standalone - it has no [shader] entry
// point, so it must live outside shaders/*.slang's non-recursive glob; see
// src/gpu/CMakeLists.txt) by each of the five wavefront_*.slang stage
// kernels.
//
// The actual shading math (struct layouts, RNG, warp/fresnel/microfacet
// helpers, texture evaluation, every BSDF, emitter sampling, BVH traversal)
// lives in scene_common.slang, which this file includes and which is shared
// verbatim with the megakernel (path_trace.slang) and real-time
// (rt_path_trace.slang) backends. It used to be duplicated here as a
// byte-identical ~1300-line copy of path_trace.slang's; that duplication is
// gone, so a BSDF fix now applies to every backend at once and the two can no
// longer silently diverge. See scene_common.slang's header for the CPU-parity
// notes (RNG choice, envmap sampling simplification, etc.).
//
// What remains in THIS file is only what is specific to the wavefront
// architecture: the additional bindings and push-constants needed to carry
// per-path state and work queues BETWEEN separate kernel dispatches instead
// of within one big loop.
//
// ============================================================================
// Wavefront architecture overview
// ============================================================================
//
// Per-pixel path state (ray, throughput, rng, etc.) is persistent across
// bounces in device buffers indexed by PIXEL index (not by "queue slot"),
// so each kernel can look up/update a path's state without needing it
// passed through a queue payload. A pixel's path is represented in the
// CURRENT bounce's *work queue* only while it is still alive; queues store
// just the pixel index (a uint), making them cheap to build/compact.
//
// Two `active_indices_*` buffers are used as ping-pong queues: on bounce D,
// one holds the list of pixels still alive entering this bounce (the
// "read" queue) and the other is atomically appended to by wavefront_shade
// to build the list of pixels alive AFTER this bounce (the "write" queue).
// `pc.read_is_a` (flipped by the host every bounce) selects which physical
// buffer plays which role, avoiding a copy between bounces.
//
// Per-bounce sequence (see src/gpu/wavefront_main.cpp for the exact
// dispatch/barrier orchestration):
//   1. wavefront_extend   - BVH-traverse each active path's current ray,
//                           write the hit result (or miss) to a slot-indexed
//                           hit buffer.
//   2. wavefront_shade    - resolve material, accumulate emission (with
//                           MIS), evaluate+sample the BSDF, and for paths
//                           that survive: atomically append to the NEXT
//                           active queue AND (if applicable) atomically
//                           append a shadow ray to the NEE queue.
//   3. wavefront_prepare_shadow_indirect - tiny (1-thread) kernel: snapshot
//                           the shadow queue's size into an indirect
//                           dispatch-args buffer for step 4. Must run
//                           strictly between shade and shadow (see its own
//                           file for why the ordering here is load-bearing).
//   4. wavefront_shadow   - trace each queued shadow ray; if unoccluded,
//                           add its pre-computed NEE contribution into
//                           accum_buffer.
//   5. wavefront_advance  - tiny (1-thread) kernel: snapshot the next-active
//                           queue's size into the indirect dispatch-args
//                           buffer for bounce D+1's extend/shade, then reset
//                           the atomic counters for D+1's shade to build into.
//
// accum_buffer[pixel] is updated via plain (non-atomic) `+=` in both
// wavefront_shade (emission) and wavefront_shadow (NEE), which is safe
// WITHOUT atomics because within a single bounce, each pixel appears in
// at most one of {a terminating emission-hit in shade, a queued shadow ray}
// - never both - and every queue (by construction) contains each currently-
// active pixel at most once. See the conversation history / PR description
// for the full data-race argument.

#include "scene_common.slang"

// Threads per 1D workgroup for every wavefront kernel except raygen (which
// dispatches over pixels directly) and the two 1-thread housekeeping
// kernels. Must match the group size the host uses to compute indirect
// dispatch arguments (see wavefront_main.cpp's GROUP_SIZE).
static const uint GROUP_SIZE = 256u;

// ============================================================================
// Bindings
// ============================================================================
// Bindings 1-15 (the shared flattened scene) come from scene_common.slang, in
// the same order the megakernel uses, so GPUScene upload code is shared
// between all backends (see src/gpu/gpu_scene_upload.cpp).

[[vk::binding(0, 0)]] RWStructuredBuffer<float4> accum_buffer; // running sum of radiance (.xyz) + sample count (.w)

// ---- Per-pixel persistent path state (indexed by PIXEL index) ----
[[vk::binding(16, 0)]] RWStructuredBuffer<float4> path_ray_o;       // .xyz = origin,    .w = t_min
[[vk::binding(17, 0)]] RWStructuredBuffer<float4> path_ray_d;       // .xyz = direction, .w = t_max
[[vk::binding(18, 0)]] RWStructuredBuffer<float4> path_throughput;  // .xyz = throughput,.w = eta
[[vk::binding(19, 0)]] RWStructuredBuffer<uint> path_rng;           // Rng.state
[[vk::binding(20, 0)]] RWStructuredBuffer<float2> path_prev_bsdf;   // .x = pdf, .y = delta (0/1)
[[vk::binding(21, 0)]] RWStructuredBuffer<float4> path_prev_p;      // .xyz = previous vertex position

// ---- Hit results (indexed by QUEUE SLOT, aligned with the active-index
// list read by that bounce's extend/shade dispatch) ----
[[vk::binding(22, 0)]] RWStructuredBuffer<uint> hit_tri;    // INVALID_INDEX => miss
[[vk::binding(23, 0)]] RWStructuredBuffer<float4> hit_tuv;  // .x=t .y=u .z=v

// ---- Active-path queues (ping-pong; see file header) ----
[[vk::binding(24, 0)]] RWStructuredBuffer<uint> active_indices_a;
[[vk::binding(25, 0)]] RWStructuredBuffer<uint> active_indices_b;

// counters[0] = size of the CURRENT bounce's active queue (read by extend/shade)
// counters[1] = size of the NEXT bounce's active queue (atomically built by shade)
// counters[2] = size of this bounce's shadow-ray queue (atomically built by shade)
[[vk::binding(26, 0)]] RWStructuredBuffer<uint> counters;

// ---- Shadow-ray (NEE) queue, indexed by QUEUE SLOT ----
[[vk::binding(27, 0)]] RWStructuredBuffer<float4> shadow_o;      // .xyz = origin,    .w = t_max
[[vk::binding(28, 0)]] RWStructuredBuffer<float4> shadow_d;      // .xyz = direction
[[vk::binding(29, 0)]] RWStructuredBuffer<uint> shadow_pixel;    // owning pixel index
[[vk::binding(30, 0)]] RWStructuredBuffer<float4> shadow_contrib; // .xyz = radiance to add if unoccluded

// Indirect dispatch args: [0..2] = VkDispatchIndirectCommand for the
// active-queue-sized dispatches (extend/shade), [3..5] = for the
// shadow-queue-sized dispatch. Written by wavefront_advance /
// wavefront_prepare_shadow_indirect respectively.
[[vk::binding(31, 0)]] RWStructuredBuffer<uint> indirect_args;

struct PushConstants {
    uint spp_index;             // seeds the RNG in raygen
    uint max_depth;
    uint rr_depth;
    uint num_lights;
    uint environment_light_id;  // INVALID_INDEX if none
    uint depth;                 // current bounce index (0-based) for this extend/shade dispatch
    uint read_is_a;             // 1 => read active_indices_a/write _b this bounce, 0 => vice versa
    uint total_pixels;          // width*height, used by raygen's bounds check
};
[[vk::push_constant]] ConstantBuffer<PushConstants> pc;

// Ping-pong queue accessors - centralizing the `read_is_a` branch here means
// every stage kernel just calls these two functions instead of repeating the
// branch (and risking the two ever falling out of sync with each other).
uint read_active(uint slot) { return pc.read_is_a != 0u ? active_indices_a[slot] : active_indices_b[slot]; }
void write_next_active(uint slot, uint pixel)
{
    if (pc.read_is_a != 0u) active_indices_b[slot] = pixel; else active_indices_a[slot] = pixel;
}
