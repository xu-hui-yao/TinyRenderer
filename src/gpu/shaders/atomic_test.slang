// Minimal atomic-counter correctness test: dispatch N threads, each does
// InterlockedAdd(counter, 1). If atomics work correctly, the final counter
// value must equal N exactly, regardless of how many threads race on it.
// If InterlockedAdd silently compiles to a non-atomic read-modify-write on
// this target, concurrent threads will overwrite each other's increments
// and the final value will be LESS than N.
[[vk::binding(0, 0)]] RWStructuredBuffer<uint> counter;

[shader("compute")]
[numthreads(256, 1, 1)]
void main(uint3 tid: SV_DispatchThreadID)
{
    uint dummy;
    InterlockedAdd(counter[0], 1u, dummy);
}
