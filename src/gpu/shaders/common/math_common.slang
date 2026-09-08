// Small shared math helpers that declare NO bindings.
//
// This exists because common/scene_common.slang - the file every path-tracing
// kernel includes - also declares the whole flattened scene on bindings 1..15.
// That is exactly what those kernels want, but the DENOISING passes
// (svgf_temporal / svgf_atrous) need a couple of these helpers while binding
// their own screen-sized textures on the same slots; including scene_common
// there produces overlapping-binding warnings and, worse, silently wrong
// descriptor wiring.
//
// So: anything with no shader-resource dependency lives here and is included by
// both sides. Keep this file free of [[vk::binding]].

#ifndef M_MATH_COMMON_SLANG
#define M_MATH_COMMON_SLANG

// Transform a point / vector by a row-major 4x4 passed as four float4 rows.
// Point gets the translation column and the homogeneous divide; vector does not.
float3 transform_point(float4 r0, float4 r1, float4 r2, float4 r3, float3 p)
{
    float4 hp = float4(p, 1.0);
    float w = dot(r3, hp);
    return float3(dot(r0, hp), dot(r1, hp), dot(r2, hp)) / w;
}
float3 transform_vector(float4 r0, float4 r1, float4 r2, float4 r3, float3 v)
{
    float4 hv = float4(v, 0.0);
    return float3(dot(r0, hv), dot(r1, hv), dot(r2, hv));
}

// Rec.709 luminance, matching Color3f::luminance() in core/spectrum.h.
float luminance(float3 c) { return dot(c, float3(0.212671, 0.715160, 0.072169)); }

// Octahedral normal encoding (Cigolle et al. 2014, "A Survey of Efficient
// Representations for Independent Unit Vectors"). Cheaper than storing xyz and
// free of the sign ambiguity that would otherwise make the temporal pass reject
// valid history on back faces. rt_path_trace encodes, the denoiser passes decode.
float2 oct_encode(float3 n)
{
    n /= (abs(n.x) + abs(n.y) + abs(n.z));
    // select() rather than `?:` - the conditional operator on vectors is
    // deprecated in Slang (it does not short-circuit, which is a real trap when
    // one branch divides). Both arms here are safe to evaluate.
    float2 folded = (1.0 - abs(n.yx)) * select(n.xy >= 0.0, float2(1.0, 1.0), float2(-1.0, -1.0));
    float2 p = select(n.z >= 0.0, n.xy, folded);
    return p * 0.5 + 0.5;
}

float3 oct_decode(float2 e)
{
    float2 f = e * 2.0 - 1.0;
    float3 n = float3(f.x, f.y, 1.0 - abs(f.x) - abs(f.y));
    float t = max(-n.z, 0.0);
    n.x += n.x >= 0.0 ? -t : t;
    n.y += n.y >= 0.0 ? -t : t;
    return normalize(n);
}

#endif // M_MATH_COMMON_SLANG
