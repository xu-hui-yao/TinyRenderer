// Stage 2, wavefront kernel 5a/5: prepare the shadow dispatch's indirect args.
//
// A single-thread (numthreads(1,1,1), dispatched as (1,1,1)) housekeeping
// kernel. Must run strictly AFTER wavefront_shade (which just finished
// atomically building the shadow-ray queue into counters[2]) and strictly
// BEFORE wavefront_shadow (which needs indirect_args[3..5] to know how many
// groups to dispatch). Deliberately does NOT reset counters[2] to 0 here -
// wavefront_shadow's own bounds check (`slot >= counters[2]`) still needs
// to read the value this kernel is reading right now. The reset happens
// later, in wavefront_advance, only after wavefront_shadow has consumed it.
#include "common/wavefront_common.slang"

[shader("compute")]
[numthreads(1, 1, 1)]
void main(uint3 tid: SV_DispatchThreadID)
{
    uint shadow_count = counters[2];
    uint groups = (shadow_count + GROUP_SIZE - 1u) / GROUP_SIZE;
    indirect_args[3] = groups;
    indirect_args[4] = 1u;
    indirect_args[5] = 1u;
}
