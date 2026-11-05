# Two-coset producer/consumer schedule

The exact-work schedule is a measured negative result. On the matched
K20 Toeplitz comparison, streaming scratch stores increase complete time
by 6.35%; cached scratch stores increase it by 55.01%. All tested outputs
agree with the existing encoder. Keep the experiment default-off and do
not extend it into another tile-size search.

## Exact map and changed schedule

The experiment retains the grouped GF256 outer, alternating pair route,
independent pair swaps, h4 inner, refresh keys, and native output order.
Fix all setup values and an input. The old schedule computes every outer
coordinate before running the inner over regions 0 through 255.

The new schedule uses the existing `coset` primitive:

1. Generate regions 0 through 127 across every source group.
2. Fence streaming stores and fill these regions' wrap padding.
3. Consume these regions with the original inner recurrence.
4. Generate regions 128 through 255 across every source group.
5. Fence streaming stores and fill the remaining wrap padding.
6. Consume the remaining regions with the same carried state.

For outer group width `Group`, the two coset offsets are 0 and
`128/Group`. Both calls retain the full scratch base. Their column stores
already use absolute region indices. No pointer is rebased before its
allocation.

At K20, each half performs 2048 four-group coset calls. The combined
schedule executes the same 4096 coset calls as the old schedule. It adds
no FFT butterflies, column mixers, or diagonal products. The 32768 inner
packets also retain their original order. One zero-initialized `State`
spans both halves; it is not reset at region 128. Thus the output map is
identical for each fixed setup, not merely equal in distribution. Existing
distance claims for that setup distribution are unchanged.

The research flags are:

    HDGROUP_COSET_PIPELINE=0          # default: original schedule
    HDGROUP_COSET_PIPELINE=1          # enable the two-coset schedule
    HDGROUP_COSET_PIPELINE_STREAM=1   # default producer stores: streaming
    HDGROUP_COSET_PIPELINE_STREAM=0   # producer stores: cached

The output-store selector remains independent. Unaligned scratch uses
cached producer stores; unaligned output uses cached output stores.
`LONG_RS` uses the original full path. When the pipeline flag is zero,
the new scheduling helpers are excluded and the original hot path remains
intact. The rebuilt default binary has the same SHA256 as the prior control.

Input/output overlap requires the full path. First-half output could
otherwise overwrite input still needed by the second outer coset. The
dispatcher checks complete ranges using unsigned address differences,
without forming potentially overflowing end addresses or byte lengths.
Equal pointers and partial overlap in either direction take the fallback.

Each inner range has compile-time region bounds. Its prefetch stops at
the ready range boundary. In particular, region 127 does not prefetch the
unproduced region 128. The producer fence precedes reading each plane's
first record for odd-pair wrap padding.

## Memory consequence and interpretation

At K20, the input contains 16 MiB and the output contains 32 MiB of
128-bit records. Each outer region contains 128 KiB of payload. Full
scratch remains 32 MiB plus 16 KiB of padding. One produced half contains
16 MiB plus 8 KiB of padding; it is not a private-cache-sized buffer.

Both schedules write and read the same 32 MiB scratch payload. Their
logical producer/consumer roundtrip remains 64 MiB. Both also write the
32 MiB final output. These are payload counts, not measured DRAM traffic.

The potential benefit was a shorter interval between producing and
consuming each half. The main competing effect is input reuse. The old
64-group tile revisits a 128 KiB input tile across its two cosets. The new
schedule revisits the 16 MiB input only after consuming half the output.
The number of input loads is unchanged, but their cache history differs.
Cached producer stores also change allocation and writeback behavior.

These effects are hypotheses for the regression. The timing experiment
did not isolate cache misses, write allocation, memory bandwidth, or
instruction scheduling. It does not establish which effect caused the
slowdown, and isolated phase times must not be added to explain it.

## Validation and matched result

The default, streaming-pipeline, and cached-pipeline builds pass `--check`
and `--large-check`. The driver compares the original all-outer/full-inner
schedule and the new schedule against the same independent dense oracle.
Checks cover packet input/state basis directions, pair-swap masks,
wrapping pairs, natural lengths, K20, guards, padding, unaligned pointers,
full input/output aliasing, and partial overlap in both directions.
Additional integer tests check overlap boundaries near the address limit.

The streaming-pipeline ASAN/UBSAN build also passes `--check`. Its linked
static library is not instrumented. These tests validate implementation
consistency; the exact-map argument above explains certificate transfer.

The coordinating native run reported the following three matched serial
run medians on 2026-10-09. All builds use the Toeplitz outer, unchanged
interpolation refresh, independent pair swaps, and streaming final output.

| Schedule | Three run medians, ms | Median, ms | Change from control |
|---|---|---:|---:|
| Original | 3.805333, 3.727739, 3.784825 | 3.784825 | — |
| Pipeline, streaming scratch | 4.113258, 4.025123, 4.011598 | 4.025123 | +6.35% |
| Pipeline, cached scratch | 5.866731, 5.878313, 5.854729 | 5.866731 | +55.01% |

All three K20 output checksums are `7701295210980077712`. No XOR2 timing
extension was run after this negative result.

Use the common build definitions in
[the native refresh receipt](HD_XOR_STATE_NATIVE.md#reproduction-and-pins),
select `HDGROUP_TOEPLITZ=1`, and leave `HDGROUP_XOR_REFRESH` disabled.
Add the pipeline definitions above for each candidate. Run `--check` and
`--large-check` before the existing full-call command:

    --bench 501 full 20 nt 4 4

Benchmarks must run serially. Setup and allocation remain outside encoding.

Frozen source SHA256 values:

- `HdGroupedOuter.h`: `62ea6a44a9d3067f2b08bc9bb935501d304b1a4af5a93cb4d2089825ca216efe`.
- `hd_grouped_forward.cpp`: `32d02ba2a134f6d3938dfda5fa1a33bf37eac3d9ce7d8d669d56f01488c82c48`.

| Build | Binary SHA256 |
|---|---|
| Original | `0b31e23d66948665135485da064868ad769d007bb5cc8a2c4c4570698f1638bc` |
| Pipeline, streaming scratch | `30f23ce9a0fc7c8842d8a8ec203afd6f28487c9f67d52ca0d850dbb93a3ee82d` |
| Pipeline, cached scratch | `c18919fc385e2edc23ba7835eebffe1d0d52679df87a4766a173927304acc5b7` |

## Constraints on finer fusion

An independently routed region cannot generally consume its first packet
after the producer finishes a small source-group tile. Its requested pairs
can lie throughout the source array. All 256 regional routes are different.
Reordering the stateful inner to follow producer completion would change
the map, not merely its storage layout.

Recomputing smaller FFT subtrees trades scratch for arithmetic. For an
outer constituent of dimension D, let L divide D and be a power of two.
Recomputing the shared ancestors separately for every L-output tile uses

    (D/L)(D-L) + (D/2) log2(L)

butterfly slots per constituent coset, instead of `(D/2) log2(D)`.
For Group=4, D=32, and a 16-region tile, L=4 and this ratio is 3.2.
Constant-zero or constant-one butterflies have different instruction
costs; this count is not a runtime estimate.

Saving shared prefixes can preserve the butterfly count. A complete
one-coset prefix frontier, however, contains 128 records per source group:
16 MiB at K20. Spilling and rereading that frontier for both cosets adds
another 64 MiB of logical traffic. Such a design exchanges one global
intermediate for another; it does not automatically delete the roundtrip.

A construction-level group-tile stream has a separate distance obstacle.
Fix a linear setup. Suppose all inputs of t complete outer groups are
processed inside one local output window, and only s field symbols carry
their influence beyond that window. The source space has dimension 128t.
The terminal-state map therefore has a kernel of dimension at least
128t-s. If 128t>s, a nonzero source produces zero terminal state. With all
other inputs zero, its nonzero output is confined to at most 256t positions.

For t=64 and s=16, this confines a nonzero word to 16384 of the 2097152
output positions, far below 35%. Keeping several independent chains only
replaces s by their total carried dimension. Rotating chain labels does
not remove this rank argument. An invertible source premix also cannot
repair a low-distance image code.

Consequently, a new cache-local construction must retain long-range source
influence or introduce enough nonlocal coupling to defeat this cut. That
requires a new proof. Preserve the exact-work negative and prioritize the
cheaper XOR2 outer's unfinished whole-code coverage before another native
dataflow experiment.
