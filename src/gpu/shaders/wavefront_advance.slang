// Stage 2, wavefront kernel 5b/5: advance to the next bounce.
//
// A single-thread housekeeping kernel, dispatched after wavefront_shadow
// has finished consuming counters[2] (see wavefront_prepare_shadow_indirect
// for why that ordering matters). Snapshots the next-active queue's size
// (built by wavefront_shade into counters[1]) into indirect_args[0..2] for
// bounce D+1's extend/shade dispatches, promotes it into counters[0] (the
// bound those dispatches will check against), and resets counters[1]/[2] to
// 0 so bounce D+1's wavefront_shade can atomically build fresh queues into
// them. The host flips pc.read_is_a between bounces D and D+1 so that
// bounce D+1's extend/shade read from whichever active_indices_* buffer
// this bounce's shade just finished writing into.
#include "common/wavefront_common.slang"

[shader("compute")]
[numthreads(1, 1, 1)]
void main(uint3 tid: SV_DispatchThreadID)
{
    uint next_count = counters[1];
    uint groups = (next_count + GROUP_SIZE - 1u) / GROUP_SIZE;
    indirect_args[0] = groups;
    indirect_args[1] = 1u;
    indirect_args[2] = 1u;

    counters[0] = next_count;
    counters[1] = 0u;
    counters[2] = 0u;
}
