#pragma once

// Denormalised triangle positions for BVH traversal.
//
// Deliberately its own header with no other declarations: the diagnostic tools
// (src/gpu/*_main.cpp) keep private copies of the packed scene structs, and
// including gpu_scene_upload.h from them collides on names like GPUTextureGPU.
// This file carries the one function they all need and nothing else, so the
// packing logic still has a single definition.

// core/gpu_scene.h only pulls in core/common.h, which merely FORWARD-DECLARES
// TArray/TRGBSpectrum/TTransform; the definitions must be named explicitly for
// this header to stand on its own (same fix as in gpu_scene_upload.h).
#include <core/array.h>
#include <core/spectrum.h>
#include <core/transform.h>

#include <core/gpu_scene.h>

#include <cstdint>
#include <stdexcept>
#include <vector>

namespace tiny_renderer::gpu {

// BVH traversal is the hottest loop in the renderer. Doing
// positions[indices[tri*3+k]] there costs three dependent gathers (an index
// load, then a vertex load that depends on it) per triangle test; pre-expanding
// the vertices turns that into one contiguous read the hardware can coalesce.
// Duplicated vertices cost memory but no correctness - it is the same data, laid
// out for the access pattern.
//
// Layout: 3 x float4 per triangle (12 floats), .xyz = position, .w unused.
inline std::vector<float> build_packed_triangles(const GPUScene &gs) {
    const size_t tri_count = gs.indices.size() / 3;
    std::vector<float> out(tri_count * 12, 0.0f);
    for (size_t t = 0; t < tri_count; ++t) {
        for (size_t k = 0; k < 3; ++k) {
            const uint32_t vi = gs.indices[t * 3 + k];
            // Bounds-checked: an out-of-range index would otherwise read past
            // the end of vertex_positions, which shows up as stray triangles
            // elsewhere in the image rather than as a crash.
            if (vi >= gs.vertex_positions.size())
                throw std::runtime_error("build_packed_triangles: triangle index out of range");
            const auto &v = gs.vertex_positions[vi];
            float *dst = &out[(t * 3 + k) * 4];
            dst[0] = v.x();
            dst[1] = v.y();
            dst[2] = v.z();
            dst[3] = 0.0f;
        }
    }
    return out;
}

} // namespace tiny_renderer::gpu
