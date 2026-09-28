# Joint tuning of packet routing and the inner

## Objective and checkpoint

Minimize precomputed, single-threaded transposed encoding time for K=2^20
and 128-bit elements, subject to a whole-code distance certificate.
Compare candidates at a common distance target and at least 40 bits of
setup-failure margin. The completed two-bit point gives a fallback at
9.25% distance and more than 49.11 bits; see
[two_bit/FIRST_CLOSURE.md](two_bit/FIRST_CLOSURE.md).
The incomplete 9.5% search is paused, with its local witnesses preserved.

The inherited inner parameters are not constraints. A more expensive
inner may permit cheaper routing and reduce total time. Conversely,
smaller batches may improve mixing enough to reduce the number of updates.
Neither tradeoff follows from inner-only timing or a selected proof case.

## Initial comparison

Keep BCH[256,128], t=128, s=19, and the existing expansion and feedback
maps. Compare packet sizes g in {2,4} and update counts r in {2,3}.
Each row has an independent coordinate permutation. Each of 256 regions
has an independent packet permutation, and each packet has an independent
lane permutation. The inner emits before applying its r transvections
and adding feedback. Its initial state is zero; there is no flush.

For g=2,r=2, this is the ensemble covered by the completed proof.
Every changed configuration needs its own proof coverage. In particular,
the three-update two-bit result at selected 10% comparison compositions
does not certify every message. Previously proved two-update sparse
bounds cannot be imported by assuming monotonicity in r.

The benchmark uses seeded samplers to instantiate these setup operations.
Timing does not validate the ideal-distribution distance theorem for a
particular seed. No heuristic bank or shared row shuffle is substituted.

## Implementation and checks

`joint.cpp` retains the fixed-width unrolled emission and zeta circuits.
Update count and packet size are template parameters; the hot path has
no virtual dispatch, type-erased callbacks, or allocation. The r=3 path
adds one update to the reversed product. The r=2 path is checked against
the existing two-update kernel as well as a dense matrix-row reference.

Packets occupy 32 or 64 bytes. Both cached and non-temporal stores are
tested because a half-cache-line store may have different memory costs.
The encoder routes into padded four-row tiles. It copies each tile
sequentially to local storage, restores each row's coordinate order, and
uses the unchanged four-row GFNI BCH transpose.

Before timing, the harness checks route bijectivity, regional occupancy,
packet membership, and local offset completeness. It checks the inner
against a dense reference and a forward/transpose identity using two
different input vectors. It then compares the full output against a
materialized route plus the original BCH circuit, including the untouched
suffix. These are correctness checks, not distance experiments.

Setup, allocation, reference checks, and checksums are excluded from timing.
Each process performs three warmups, followed by in-place encodings without
input reset. Separate phase means include warmups. They need not sum to
the median full-encoder time.

## Reproduction

Build the workstream's `joint` target in Release mode with GCC and
`SPIN_TUNE=znver4`. The isolated executable requires AVX-512F/VL/BW and GFNI.
For sanitizer checks, use a separate build with
`-O1 -fsanitize=address,undefined -fno-omit-frame-pointer`.

```sh
bash research/workstreams/permutation_locality/run_joint.sh ROOT check
bash research/workstreams/permutation_locality/run_joint.sh ROOT sanitize
bash research/workstreams/permutation_locality/run_joint.sh ROOT screen
bash research/workstreams/permutation_locality/run_joint.sh ROOT confirm
```

Here ROOT contains `build` and `build-sanitize`. The script acquires all
three shared benchmark locks and pins timing runs to CPU 15. Run these
commands serially; do not compile during measurement. `screen` uses two
seeds and 31 calls per configuration. `confirm` uses three repetitions
and 101 calls for the faster store policy at each packet size, reversing
the configuration order in the middle repetition.
Raw samples and build outputs remain outside Git.

For the four-bit proof screen, run:

```sh
python -B research/workstreams/permutation_locality/independent_rows/candidates/round_screen.py --rounds 2 3 --tilts .056 --groups 64 --supports 192 200
```

This regenerates the actual local operators for each update count. It
tests selected homogeneous support events at 10% distance, not every
support vector or occupancy. Its binary64 scores are proposals, not
outward certificates. Regression tests check that the requested update
count reaches every round-dependent refinement, including the density
record's provenance. The original driver still defaults to two updates.

## Decision after the screen

First identify the fastest store policy for each packet size. Compare the
cost of a third update against its improvement at the known proof
bottlenecks. Then consider alternative batch widths and state sizes for
the competitive routes. Changed maps require new expansion spectra,
feedback character counts, and sparse operators; the current proof's
fixed geometry must not be reused silently.

Keep the earlier shared-shuffle and four-bit proof work intact. A full
certificate is attempted only after correctness, a useful timing result,
and promising proof diagnostics identify a candidate.

## Initial measurements (2026-09-28)

Peach Ryzen 7950X, GCC 15.2, Release/znver4, CPU 15, K=2^20,
128-bit elements. All runs were serial and protected by the three shared
locks. No builds ran during measurement. The table reports the mean and
range of six process medians: two route/mask seeds, three repetitions,
and 101 timed calls after three warmups. Ranges are not confidence intervals.

| Packet bits | Updates | Store policy | Mean median (ms) | Range (ms) |
|---:|---:|---|---:|---:|
| 2 | 2 | Cached | 13.384 | 13.240--13.523 |
| 2 | 3 | Cached | 13.969 | 13.824--14.179 |
| 4 | 2 | Streaming | 6.442 | 6.406--6.473 |
| 4 | 3 | Streaming | 6.727 | 6.698--6.763 |

The third update costs about 0.285 ms (4.4%) for the four-bit route.
That route is about twice as fast as the current two-bit implementation.
These measurements are not lower bounds on either construction's cost.
In particular, no optimized multistage routing scheme for two-bit packets
was implemented in this initial comparison.

The two-seed, 31-call screen found streaming stores worse for two-bit
packets (14.71 ms with two updates), whereas cached stores were worse
for four-bit packets (12.30 ms). Phase timing puts the two-bit cached
inner-plus-route near 9.5 ms and the four-bit streaming equivalent near
2.75 ms. Their local rearrangement-plus-BCH costs are approximately
3.9 and 3.7 ms. This locates the observed gap in the route/inner phase.
The explanation in terms of half-cache-line traffic remains an inference;
no hardware-counter diagnosis was performed.

All sixteen parameter/seed/store cases passed the reference checks in
Release and separately under AddressSanitizer/UndefinedBehaviorSanitizer
at K=2^14. Every timed K=2^20 process repeated the reference checks.
The proof regression suites also pass: 89 two-bit tests and 47 four-bit
candidate tests, including the two new update-count propagation tests.
The latter tests check integration of existing bounds, not a new distance
certificate for three updates.

The measurement executable SHA-256 is
`e6bd184162acefa54a75dee7ae17bf9ea8b7fe8486fbfeaccb7d2b5c6cb55eea`.
The unchanged production library SHA-256 is
`ecbebbfa16d15b4b5460c8811053720d1dd534f2aa81ea5b82757bebd153f0d7`.
The isolated remote root is `/tmp/spin-joint-9n57TT`; raw samples are in
`measurements/joint-screen-t80iRT` and `measurements/joint-confirm-zA6P6C`.
They are not part of the Git commit.

The next proof candidate is the four-bit route with three updates.
Its 6.73 ms measurement is still above the approximate 5 ms objective.
First test harder proof cases at 9.25% and 10%, retaining two updates
as a cheaper alternative if its bound suffices. Then compare batch widths
before investing in a full cover. The completed two-bit theorem remains
the proof fallback, not the current performance winner.

## Initial proof sensitivity

The four-bit screen completed at output-weight cutoff 209715 (10%),
64 active groups, output tilt .056, and all-one penalty .9. It regenerated
the local operators for each update count using the six-window feedback
census, eight-window output moments, three-window joint cancellation,
and translated-density refinement. Direct empty/single-window component
checks passed for both update counts.

Here u is the common union size of the shuffled row supports in each
active group. The numbers are binary64 proposals for the log2 of an
upper bound on the selected class's expected bad-message count.

| u | Two updates | Three updates |
|---:|---:|---:|
| 192 | -114.36 | -558.55 |
| 200 | -124.71 | -571.03 |

The extra update improves these selected bounds by approximately
444--446 bits for 4.4% measured runtime overhead. Both two-update cases
already pass this diagnostic; this is evidence of added slack, not
closure of an unresolved event. Nor are these whole-code margins.
There is no outward replay or complete support/occupancy cover for the
three-update construction. Next include the harder 80-group cases,
then the intermediate and dense ranges, before choosing a full-proof run.
