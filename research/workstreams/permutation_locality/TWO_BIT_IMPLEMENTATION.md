# Tiled implementation of the proved two-bit construction

The initial 13.38 ms result measured direct scatter, not an optimized
two-bit encoder. The older tiled experiment already showed comparable
one-bit and two-bit timings with a one-update inner. Packet size alone
therefore did not explain the regression.

This implementation retains the two-bit construction in
[FIRST_CLOSURE.md](two_bit/FIRST_CLOSURE.md): BCH[256,128], independent row
shuffles, independent regional packet permutations and lane swaps, and
IMT(128,19) with two updates. At K=2^20, its existing ideal-setup theorem
gives distance greater than 9.25%, except with probability below 2^-49.11.
The implementation changes the execution schedule, not that distribution.
It does not claim that a chosen benchmark seed has been distance-certified.

## Implementation

`TwoBitTiled.h` converts an already sampled inverse route into bucket slots
and local BCH offsets. Each bucket contains a contiguous group of outer
rows. The two rows in a packet always belong to the same bucket.

The transpose executes in two stages:

1. Run the unchanged, unrolled two-update inner backward. Assemble each
   pair of 128-bit values in a 256-bit register and write it to its bucket.
   Slots preserve inner order within each bucket, so each bucket receives
   sequential writes. The lane order is restored in the next stage.
2. Read each bucket sequentially and scatter into its local BCH tile.
   Precomputed offsets restore the exact row coordinates and lane order.
   Apply the existing four-row GFNI BCH circuit and write the compressed
   output in place.

The selected layout uses 16 rows per bucket and one cache-line pad between
buckets. The pad changes addresses, not code coordinates. The local tile
is 64 KiB. At K=2^20, the routing metadata occupies 12 MiB: one 32-bit
slot per pair and one 32-bit local offset per value. Bucket storage needs
32 MiB plus 32 KiB of padding. The diagnostic harness retains additional
reference buffers; these are not requirements of the encoding kernel.

All allocation and setup occur before encoding. Packet assembly, the
inner recurrence, and BCH evaluation use fixed-width SIMD. The hot path
has no dynamic dispatch, allocation, or type-erased callback.

## Measurement protocol

Measurements use Peach's Ryzen 7950X, GCC 15.2, Release/znver4, and CPU 15.
The workload is a precomputed, single-threaded, in-place transpose at
K=2^20 with 128-bit elements. Setup, allocation, reference checks, and
checksums are excluded. Three warmups precede timing; inputs are not
reset between calls. Three shared locks prevent concurrent benchmarks.
Compilation is also excluded from measurement periods.

The 31-call, two-seed screen found approximately 8.5 ms with an unpadded
512-row layout. Padded 16-row buckets reduced this to 7.51--7.55 ms.
Padding is not uniformly beneficial: the padded 128-row candidate was
slower. Tile geometry must be measured rather than inferred from packet
width. These screens select the candidate; they are not confidence intervals.

The confirmation compares both old and new implementations,
plus scalar-store and original-BCH controls. It uses two setup seeds,
three processes per configuration and seed, and 101 timed calls per
process. The second repetition reverses the configuration order.
The production one-bit baseline keeps its original one-update inner;
only the two-bit controls have an identical route and two-update inner.

## Confirmed results (2026-09-28)

The table gives the mean and range of six process medians. The ranges
describe these measurements; they are not confidence intervals.
All tiled two-bit controls below use 16 rows per bucket.

| Implementation | Updates | Mean median (ms) | Range (ms) |
|---|---:|---:|---:|
| Production one-bit route and original BCH | 1 | 9.646 | 9.583--9.689 |
| Two-bit direct scatter and GFNI BCH | 2 | 13.424 | 13.342--13.487 |
| Two-bit tiled scalar stores and GFNI BCH, unpadded | 2 | 9.207 | 9.111--9.350 |
| Two-bit tiled pair stores and original BCH, unpadded | 2 | 9.637 | 9.583--9.686 |
| Two-bit tiled pair stores and GFNI BCH, unpadded | 2 | 8.844 | 8.816--8.868 |
| **Two-bit tiled pair stores and GFNI BCH, padded** | **2** | **7.422** | **7.386--7.465** |
| Same padded two-bit implementation, diagnostic only | 1 | 7.226 | 7.180--7.281 |

The selected implementation takes 23.1% less time than the production
one-bit encoder, a 1.30x throughput improvement. It takes 44.7% less time
than the direct two-bit implementation. This is an implementation comparison,
not evidence that the two distributions have identical distance guarantees.

The controls separate three gains: pair stores, GFNI BCH, and bucket
padding. With the selected padded layout, removing the second update
saves only 0.196 ms, about 2.6% of the two-update time. That diagnostic
does not inherit the two-update distance certificate. The second update
was not the cause of the original 13 ms result.

For each setup seed, all two-update implementations have identical
checksums after the full sequence of 104 encodings. The one-update
control and production one-bit baseline intentionally implement different
maps and do not share those checksums.

The initial padded-layout phase screen measured roughly 3.4 ms for
inner-plus-bucket writes and 4.2 ms for local restoration-plus-BCH.
Those phase means include warmups and timing overhead; they are not a
decomposition of the final 7.422 ms median.

## Checks and reproduction

The harness checks route bijectivity, regional occupancy, and packet
membership. It compares the unrolled inner against a dense matrix-row
reference, checks the forward/transpose identity, and compares complete
encoding against the original BCH circuit. It checks the untouched suffix.
Check-only runs replay two dense inputs and a sparse input. Additional
layout tests replay coordinate destinations and reject malformed packets,
duplicate coordinates, invalid sizes, and out-of-range coordinates.

All 114 cases in `run_two_bit.sh check` passed in Release and separately
under AddressSanitizer/UndefinedBehaviorSanitizer. The matrix includes
108 small cases and six larger cases through K=2^20. Each case checks
three input patterns. The earlier sixteen joint-screen cases also pass
in both builds. Every timing process runs reference checks before timing.

Build the existing research targets `joint` and `locality`. The latter
provides the unchanged production baseline. For sanitizers, build `joint`
separately with AddressSanitizer and UndefinedBehaviorSanitizer enabled.
Then run these commands serially:

```sh
bash research/workstreams/permutation_locality/run_two_bit.sh ROOT check
bash research/workstreams/permutation_locality/run_two_bit.sh ROOT sanitize
bash research/workstreams/permutation_locality/run_two_bit.sh ROOT screen
bash research/workstreams/permutation_locality/run_two_bit.sh ROOT padded
bash research/workstreams/permutation_locality/run_two_bit.sh ROOT confirm 16
```

`ROOT` contains `build` and `build-sanitize`. A nonzero lock-wait limit can
be supplied through `SPIN_LOCK_WAIT`. Raw results go under `ROOT/measurements`,
not into Git. The single-configuration command is:

```sh
ROOT/build/joint 20 2 2 5 17 101 0 16
```

Its arguments are message exponent, packet size, update count, execution
mode, seed, timed calls, phase-timing flag, and tile rows. Mode 5 selects
padded tiled pair stores; zero calls select correctness checks only.

The remote root is `/tmp/spin-joint-9n57TT`. Confirmation logs are in
`measurements/two-bit-confirm-xcGA4S`; Release and sanitizer checks are in
`measurements/two-bit-check-dem8lT` and `measurements/two-bit-sanitize-Eyn7D4`.
The tile screens are in `measurements/two-bit-screen-4igSur` and
`measurements/two-bit-padded-TbtH8L`. These raw files remain outside Git.

Confirmed `joint` executable SHA-256:
`2e5106abb779b92472272b8b0ca2772d1b65ba371ae97b84a5f88bcbbc85de68`.
Production-control `locality` executable SHA-256:
`c1a215d6f527c731133adaacdf719bb3584a1e64b5f4debffb79cf2fc23e626a`.
Unchanged production library SHA-256:
`ecbebbfa16d15b4b5460c8811053720d1dd534f2aa81ea5b82757bebd153f0d7`.

This remains a research implementation. Public library defaults and the
paper are unchanged. The approximately 5 ms objective is not yet met.
The next optimization should target the remaining routing/BCH work while
retaining this proved construction as the measured control.
