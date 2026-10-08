# Speed-first GF256 SPIN: 83.33% rate with a smaller outer and inner

Historical first checkpoint. The [exact-map follow-up](GF8_SPEED_ROUTED.md)
now selects a sixteen-packet pipeline with relocated diagonals and a fused
outer, measuring about3.2ms. The construction and certificate below remain
unchanged; this report preserves the earlier measurements.

The new recommendation is GRS[64,40] with a four-input/four-state inner,
retaining three outputs per step. Its complete K20 certificate gives 10%
symbol distance and 57-bit setup-failure margin. The selected implementation
encodes in about 4.1 ms, versus 4.9 ms for the previous 84% checkpoint.
Production and the earlier certified implementations are unchanged.

The [construction and proof](GF8_SPEED_DESIGN.md) are separate from this
implementation report. These timings measure precomputed transposed encoding
of 128-bit records, with sixteen parallel GF256 lanes. Setup, allocation,
payload restoration, and output checksums are excluded; all encoding and
stores are included. This is not a forward-encoder or application benchmark.

## Matched K20 results

Both profiles encode requested dimension 1,048,576. The smaller design has
1,259,520 retained records, versus 1,254,400 for the old 84% profile: only
0.4082% more. Its nominal rate is 83.3333%; padding lowers requested rate
to 83.2520%. The prior profile's requested rate is 83.5918%.

Three separate launches used 101 full calls per candidate, rotating candidate
order and restoring the payload outside timing before each batch of five.
All candidates shared input and scratch addresses within a launch.

| Configuration | Full-call launch medians | Representative p90 | Status |
|---|---:|---:|---|
| Previous 84% implementation | 4.821--4.921 ms | 4.948 ms | Preserved control |
| New 83.33%, affine maps and packed diagonals, NT regional stores | 4.279--4.433 ms | 4.525 ms | Stable alternative |
| New 83.33%, scatter | 3.980--4.042 ms | 4.970 ms | Lower median, worse tail |
| **New 83.33%, cached regional stores, eight-packet batching** | **4.064--4.131 ms** | **4.340 ms** | **Recommended K20 schedule** |

The representative p90 column comes from the first of these three full-call
launches, not the best p90 across runs. Exact launch medians, in microseconds:

```text
old84:                 4821.110, 4921.367, 4874.630
new40-affine-packed:    4433.256, 4425.401, 4278.827
new40-scatter:          4042.217, 4036.315, 3979.720
new40-cached-eight:     4130.751, 4109.823, 4064.397
```

The selected schedule reduces time by about 14--17% in paired launches.
The scatter median alone would suggest a larger win, but its tail and
payload-restoration sensitivity make it a weaker default.

### Restore the input before every call

An additional 101-call experiment restores the same original payload outside
the timer before every call. This changes cache pressure, not the code or
the timed work. It is not a claim to measure cold memory or a full producer.

| Configuration | Median | p90 |
|---|---:|---:|
| Previous 84% | 4.848 ms | 4.907 ms |
| New affine/packed, NT regional | 4.340 ms | 4.410 ms |
| New scatter | 4.336 ms | 5.110 ms |
| **New cached regional, eight-packet** | **4.189 ms** | **4.269 ms** |

The recommended path retains a 13.6% improvement under this protocol.

### Final entry-point check

After selecting the schedule, the full-call and restore-every-call drivers
were changed to invoke `rs4TransposeCached` directly. Release correctness
checks were repeated before a final pair of serial 101-call measurements.
This checks the selected entry point, not just separately invoked phases.

| Protocol | Previous 84% median | Selected entry-point median | Selected p90 |
|---|---:|---:|---:|
| Restore before each five-call batch | 4.922 ms | 4.095 ms | 4.315 ms |
| Restore before every call | 5.095 ms | 4.241 ms | 4.328 ms |

The selected entry point preserves the measured improvement. Across the
restore-every-call launches, its median is 4.189--4.241 ms.

## What changed

The construction change halves the outer length and the inner width/state.
The outer uses the first 64-point subtree of the existing additive FFT;
only the first 40 source coefficients are returned. Its final butterfly is
fused with stores. Source coordinates remain four-row-interleaved normalized
additive-polynomial coefficients, with a shortened requested-dimension suffix.

The normalized four-point inner needs five constant products per circulant,
instead of nineteen for the eight-point maps. Counting the state update,
feedthrough scale, and input diagonal gives 22 GFNI sites per sixteen raw
records, versus 31 for the old inner. This is a static circuit count, not a
prediction of whole-encoder speedup. The new code processes 4.59% more raw
records, and routing still costs work.

The selected exact-map execution schedule uses:

- GFNI affine instructions for constant field multiplications; variable
  diagonal entries still use GFNI byte multiplication.
- One compact diagonal load per four-packet tile, followed by register
  broadcasts, rather than separate loads for each four-coefficient group.
- Eight-packet batches, exposing independent emission work before the
  reverse state chain. The recursive state continues across batches and regions.
- Folded local layouts, with the fixed record permutation absorbed into
  precomputed diagonal and routing tables.
- Cached regional stores, then non-temporal final output stores.

Each item after the construction change preserves the sampled linear map.
The selected entry point is `rs4TransposeCached` in `RateSpeed4Code.h`.
It preserves the unused suffix when input and output alias exactly. Scratch
is caller-owned and disjoint from input/output; the hot path allocates nothing.
The timing driver is `spin_char2_rate_speed4`.

## Rejected or nondecisive experiments

The old 84% code supplied useful controls. Tiling the routed planes into
16, 64, 256, or 1,024-row tiles did not improve full encoding; some tiles
roughly doubled outer time. A small routing prefetch was neutral. Fusing
the outer's right leaves with parent stages reduced function calls and
scratch traffic but did not improve matched full time. Constant affine maps
alone saved only about 2% on that older code.

In the new design, eight-packet scheduling with non-temporal regional stores
was not itself a win. It became useful with cached regional stores. This
interaction is why the recommendation names a complete schedule, not a list
of independently assumed speedups. A phased cached-eight run measured about
2.694 ms for inner/routing and 1.396 ms for the outer. The preceding
NT affine/packed path measured about 2.535 and 1.848 ms in the same launch.
Cached stores make the first phase dearer but the complete call cheaper.

The [64,39] fallback was implemented as a control. It expands more and was
slower; no whole-code certificate is claimed for it in this report.

## Validation and reproduction

Release validation passes for both experimental outer dimensions. The
independent scalar oracle constructs and normalizes the raw Cauchy matrix
with Gauss--Jordan elimination, then evaluates the outer polynomials directly.
It does not use optimized inner constants or FFT circuits.

Checks cover all byte values in each local input/state coordinate, dense
packet adjoints, natural-size boundaries, partial four-record tails, full
scalar transpose, forward/transpose identities, setup folding, compact and
wide inverse indices, input preservation, exact aliases, suffix preservation,
unaligned offsets, guards, and tight/page-multiple scratch strides.
Large checks cover K16/K18/K20 and all four/eight-packet and routing controls.
Only the stated K20 [64,40] configuration has the new distance certificate.
The final ASAN/UBSAN replay passes the independent small checks and full
K16/K18/K20 checks for both new drivers, including the selected cached wrapper.

Measurements use Peach's Ryzen 9 7950X, GCC 15.2, AVX512/GFNI, CPU2, and boost
disabled. Buffers at least 16 MiB request huge pages. The scratch column
stride is a 4 KiB multiple. All native benchmarks run serially under a lock.

```text
cmake -S spin/experiments/char2_field_spin -B build/field -DCMAKE_BUILD_TYPE=Release
cmake --build build/field --target spin_char2_rate_speed4 -j2
build/field/spin_char2_rate_speed4 --check
build/field/spin_char2_rate_speed4 --large-check
taskset -c 2 build/field/spin_char2_rate_speed4 1048576 101 full
taskset -c 2 build/field/spin_char2_rate_speed4 1048576 101 fresh
taskset -c 2 build/field/spin_char2_rate_speed4 1048576 101 phases
```

The proof replay and twelve proof tests are documented in
[the design note](GF8_SPEED_DESIGN.md). Their complete coverage and exact
map identities also received an independent review.

Next, optimize the coupled inner/routing path of this certified design.
It still accounts for roughly two thirds of full time. Keep the cached-eight
implementation and restore-every-call comparison as controls; do not trade
the distance guarantee for an unproved routing heuristic.
