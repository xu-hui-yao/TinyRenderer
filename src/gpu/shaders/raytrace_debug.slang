// Stage 1, checkpoint 2: primary-ray BVH traversal + geometric-normal
// visualization.
//
// Verifies that GPUScene's flattened geometry + BVH (see core/gpu_scene.h,
// Scene::build_gpu_scene(), BVHAccel::export_gpu_bvh()) can be uploaded to
// and correctly traversed on the GPU. This is NOT a path tracer yet - there
// is no material/light shading here, just "does the traced geometry look
// like the input scene" (visualized as the hit triangle's geometric
// normal). Full shading (materials, lights, multiple bounces) is the next
// checkpoint.
//
// Camera ray generation mirrors PerspectiveCamera::sample_ray() exactly
// (see src/sensors/perspective.cpp): same sample_to_camera / camera_to_world
// matrices, applied with the same row-major convention (see GPUCamera's
// comment in core/gpu_scene.h for why the matrices are decomposed into row
// vectors instead of using a native float4x4 - this sidesteps any ambiguity
// about this shading language's default matrix layout/multiplication
// order). BVH traversal (explicit stack, near-child-first ordering by split
// axis) mirrors BVHAccel::ray_intersect() (see src/accelerate/bvh.cpp).

struct GPUCameraGPU {
    float4 c2w0, c2w1, c2w2, c2w3; // camera_to_world, row-major (row i in c2wI)
    float4 s2c0, s2c1, s2c2, s2c3; // sample_to_camera, row-major
    uint width, height;
    float near_clip, far_clip;
};

// Mirrors core/gpu_scene.h's GPUBVHNode exactly (32 bytes, tightly packed:
// see that header's static_assert). meta's low 30 bits are the leaf
// primitive count (0 => interior node); its top 2 bits are the split axis.
struct GPUBVHNodeGPU {
    float3 bmin;
    uint offset;
    float3 bmax;
    uint meta;
};

static const uint BVH_META_COUNT_MASK = 0x3FFFFFFFu;
static const uint BVH_META_AXIS_SHIFT  = 30u;

[[vk::binding(0, 0)]]
RWStructuredBuffer<float4> output_image;

[[vk::binding(1, 0)]]
ConstantBuffer<GPUCameraGPU> camera;

[[vk::binding(2, 0)]]
StructuredBuffer<float4> positions; // xyz used, w unused (padding to dodge any vec3-array packing ambiguity)

[[vk::binding(3, 0)]]
StructuredBuffer<uint> indices; // 3 per triangle, global vertex indices

[[vk::binding(4, 0)]]
StructuredBuffer<GPUBVHNodeGPU> bvh_nodes;

[[vk::binding(5, 0)]]
StructuredBuffer<uint> bvh_primitives; // global triangle ids, indexed by leaf node offset/count

// Row-major 4x4 * homogeneous point (w=1), returns xyz after the perspective divide.
float3 transform_point(float4 r0, float4 r1, float4 r2, float4 r3, float3 p)
{
    float4 hp = float4(p, 1.0);
    float w = dot(r3, hp);
    return float3(dot(r0, hp), dot(r1, hp), dot(r2, hp)) / w;
}

// Row-major 4x4 * direction vector (w=0); no perspective divide.
float3 transform_vector(float4 r0, float4 r1, float4 r2, float4 r3, float3 v)
{
    float4 hv = float4(v, 0.0);
    return float3(dot(r0, hv), dot(r1, hv), dot(r2, hv));
}

// Slab test against an axis-aligned box. rd_inv = 1/ray_direction (may
// contain +-inf along axes where rd is 0; relies on IEEE-754 min/max
// handling this correctly, same trick used by most GPU ray tracers).
bool ray_aabb(float3 ro, float3 rd_inv, float3 bmin, float3 bmax, float t_min, float t_max,
             out float near_t, out float far_t)
{
    float3 t0 = (bmin - ro) * rd_inv;
    float3 t1 = (bmax - ro) * rd_inv;
    float3 tsmall = min(t0, t1);
    float3 tbig   = max(t0, t1);
    near_t = max(max(tsmall.x, tsmall.y), max(tsmall.z, t_min));
    far_t  = min(min(tbig.x, tbig.y), min(tbig.z, t_max));
    return near_t <= far_t;
}

// Standard Moeller-Trumbore ray-triangle intersection.
bool ray_triangle(float3 ro, float3 rd, float3 p0, float3 p1, float3 p2, out float t, out float u, out float v)
{
    float3 e1 = p1 - p0;
    float3 e2 = p2 - p0;
    float3 pvec = cross(rd, e2);
    float det = dot(e1, pvec);
    t = 0.0; u = 0.0; v = 0.0;
    if (abs(det) < 1e-9)
        return false;
    float inv_det = 1.0 / det;
    float3 tvec = ro - p0;
    u = dot(tvec, pvec) * inv_det;
    if (u < 0.0 || u > 1.0)
        return false;
    float3 qvec = cross(tvec, e1);
    v = dot(rd, qvec) * inv_det;
    if (v < 0.0 || u + v > 1.0)
        return false;
    t = dot(e2, qvec) * inv_det;
    return t > 1e-6;
}

[shader("compute")]
[numthreads(16, 16, 1)]
void main(uint3 tid: SV_DispatchThreadID)
{
    if (tid.x >= camera.width || tid.y >= camera.height)
        return;

    // ---- Camera ray generation (mirrors PerspectiveCamera::sample_ray) ----
    float px = float(tid.x) + 0.5;
    float py = float(tid.y) + 0.5;
    float inv_w = 1.0 / float(camera.width);
    float inv_h = 1.0 / float(camera.height);

    float3 near_p = transform_point(camera.s2c0, camera.s2c1, camera.s2c2, camera.s2c3,
                                    float3(px * inv_w, py * inv_h, 0.0));
    float3 d = normalize(near_p);
    float inv_z = 1.0 / d.z;

    float3 ro = transform_point(camera.c2w0, camera.c2w1, camera.c2w2, camera.c2w3, float3(0.0, 0.0, 0.0));
    float3 rd = transform_vector(camera.c2w0, camera.c2w1, camera.c2w2, camera.c2w3, d);

    float t_min = camera.near_clip * inv_z;
    float t_max = camera.far_clip * inv_z;

    // ---- BVH traversal: explicit stack, near-child-first (mirrors bvh.cpp) ----
    float3 rd_inv = float3(1.0, 1.0, 1.0) / rd;

    uint stack[64];
    int sp = 0;
    uint node_idx = 0;

    bool hit = false;
    float best_t = t_max;
    float3 best_normal = float3(0.0, 0.0, 0.0);

    for (;;)
    {
        GPUBVHNodeGPU node = bvh_nodes[node_idx];
        float near_t, far_t;
        bool box_hit = ray_aabb(ro, rd_inv, node.bmin, node.bmax, t_min, best_t, near_t, far_t);

        if (box_hit)
        {
            uint count = node.meta & BVH_META_COUNT_MASK;
            if (count > 0)
            {
                uint base = node.offset;
                for (uint i = 0; i < count; ++i)
                {
                    uint tri = bvh_primitives[base + i];
                    float3 p0 = positions[indices[tri * 3 + 0]].xyz;
                    float3 p1 = positions[indices[tri * 3 + 1]].xyz;
                    float3 p2 = positions[indices[tri * 3 + 2]].xyz;

                    float t, u, v;
                    if (ray_triangle(ro, rd, p0, p1, p2, t, u, v) && t < best_t && t >= t_min)
                    {
                        best_t = t;
                        best_normal = normalize(cross(p1 - p0, p2 - p0));
                        hit = true;
                    }
                }

                if (sp == 0)
                    break;
                node_idx = stack[--sp];
            }
            else
            {
                uint axis = node.meta >> BVH_META_AXIS_SHIFT;
                bool neg = rd[axis] < 0.0;
                uint near_child = neg ? (node_idx + 1) : node.offset;
                uint far_child  = neg ? node.offset : (node_idx + 1);
                stack[sp++] = far_child;
                node_idx = near_child;
            }
        }
        else
        {
            if (sp == 0)
                break;
            node_idx = stack[--sp];
        }
    }

    float3 color;
    if (hit)
    {
        // Map normal components from [-1,1] to [0,1] for visualization.
        color = best_normal * 0.5 + float3(0.5, 0.5, 0.5);
    }
    else
    {
        color = float3(0.05, 0.05, 0.08); // background
    }

    uint idx = tid.y * camera.width + tid.x;
    output_image[idx] = float4(color, 1.0);
}
