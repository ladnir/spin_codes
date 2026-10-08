# Native comparison near the GF256 rate frontier

This is the preserved rate-climb checkpoint. The later
[speed-first comparison](GF8_SPEED_NATIVE.md) selects a smaller 83.3333%
construction at about4.1ms; the measurements and certificates below are unchanged.

The 84% periodic profile is the current engineering recommendation at K20.
It keeps the frozen eight-input/eight-state map, retains outputs in the period
`(6,6,6,7)`, and uses GRS[128,84]. The complete certificate gives 10% symbol
distance and 115 bits of setup-failure margin. The slightly higher-rate
`(5,6)` profile remains a useful proof checkpoint, not the measured speed winner.

This report concerns precomputed transposed encoding of 128-bit application
records. Each record carries sixteen GF256 lanes of the same scalar code.
It does not measure setup, forward encoding, or a complete application.
The new research driver is `spin_char2_rate_climb8`; production is unchanged.

## Geometry and repeated full calls

Requested dimension is `K=2^20=1,048,576` for every row. All new configurations
have complete [K20 certificates](GF8_RATE_FRONTIER.md). The code still has one
regional routing stage. The retained-output length, not raw inner length,
determines rate and the 10% distance threshold.

The following three-launch ranges use 4 KiB-multiple column strides, folded
routing, fused outer loads, and a fused final outer butterfly/store. Each
launch contains 101 full calls per candidate. The 75% frozen control uses
its unchanged implementation and gather schedule.

| Profile | Nominal rate | Requested rate | Raw symbols | Retained records | Certified margin | Full-call launch medians |
|---|---:|---:|---:|---:|---:|---:|
| Frozen seven-output control | 75% | 74.6356% | 1,605,632 | 1,404,928 | 116 bits | 5.270--5.308 ms |
| Fixed six outputs, outer80 | 83.3333% | 83.2520% | 1,679,360 | 1,259,520 | 264 bits | 4.903--4.918 ms |
| Period `(6,6,6,7)`, outer84 | 84% | 83.5918% | 1,605,632 | 1,254,400 | 115 bits | **4.807--4.849 ms** |
| Period `(5,6)`, outer74 | 84.0909% | 83.8657% | 1,818,624 | 1,250,304 | 183 bits | 5.318--5.350 ms |

The 84% profile saves 10.714% of retained records versus the 75% checkpoint,
with about 9% lower encoding time in this matched comparison. The 84.09%
profile saves only another 0.3265% of records, while processing 13.27% more raw
symbols and taking about 10.4% longer. Thus a higher final rate need not give
a faster transposed encoder. The 84% profile is also faster and higher-rate
than fixed-six, while both exceed the 40-bit margin target.

Exact 101-call medians in microseconds, in launch order:

```text
frozen75:       5269.727, 5304.673, 5307.638
fixed-six:      4917.630, 4904.576, 4902.613
periodic84:     4807.294, 4849.363, 4821.140
periodic74:     5324.920, 5350.087, 5318.048
```

The refactored seven-output control agrees bit-for-bit with the frozen encoder,
including seeded route/diagonal tables, shortening, and aliases. In a separate
three-launch comparison with folded routing and the materialized outer stage,
its medians were 5023.117, 5013.929, and 5051.038 microseconds; the frozen
control measured 5282.812, 5273.714, and 5307.638. This isolates an exact-map
layout improvement from the change in code rate.

## Implementation and A/B results

The new code is isolated in `RateClimb8Code.h` and `RateClimb8Inner.h`.
Profile types select the outer dimension and puncturing period at compile time.
The hot kernel preserves the frozen unrolled GFNI circulants and backward
state recurrence. It makes no allocations or indirect callback calls.

Retained records are compact within four chronological packets, in coordinate
order and then phase order. Fixed-six uses 24 records per tile. The 84% profile
has 24 common records followed by coordinate six of phase three. The 84.09%
profile has 20 common records followed by coordinate five of phases one and
three. Partial loads never read into a discarded or nonexistent record.
This serialized layout is explicit; a conversion from packet-major inputs
is neither performed nor hidden outside the measured call.

### Folding the local layout into setup

The original packet hook transposes its final coordinate-major vectors back
into packet-major order before routing. The new variant leaves them in place
and precomputes the corresponding route and input-diagonal order.
For original local index `i=32*t+8*phase+coordinate`, the packed slot is
`32*t+4*coordinate+phase`. State diagonals remain chronological.
Forward and inverse route tables are transformed consistently; the sampled
linear map and setup distribution do not change.

This removes sixteen record-shuffle instructions per four-packet batch.
An isolated phased comparison reduced the 84% inner-plus-route time from
3.216 ms to 2.996 ms. Repeated full-call results above use this variant.
The original tables and unfurled schedule remain available as controls.
The research plan currently prepares both sets of tables; reducing setup
storage after selecting a schedule is separate from the measured encode gain.

### Column stride and outer final stores

Column stride was a substantial performance variable, not merely bookkeeping.
At K20, the outer84 profile has 12,544 rows: its tight stride is already
49 times 4 KiB. Outer80 and outer74 have 13,120 and 14,208 rows, whose tight
strides end 1,024 and 2,048 bytes into a page. Rounding their scratch strides
to 13,312 and 14,336 records aligns column starts to 4 KiB boundaries with
the page-aligned K20 scratch allocation. The `page` selector guarantees the
stride multiple, not page alignment of smaller allocations' base addresses.
This changes neither the code nor the number of live symbols and adds at
most 1.47% to these profiles' column scratch storage.

With tight columns, folded routing, and materialized final stores, the best
measured fixed-six and periodic74 schedules took 6.168--6.224 and
6.371--6.402 ms. Page-aligned strides brought fused-load schedules to
4.998--5.038 and 5.424--5.446 ms. The outer84 stride is unchanged;
its corresponding range was 4.859--4.929 ms. These are measured layout
effects; the underlying cache/TLB mechanism has not been isolated.

Padding every column by only four records was a negative control. It did
not repair the slower profiles and increased the outer84 phased time to
about 6.3--6.5 ms. Arbitrary padding is not interchangeable with the selected
page-aligned stride.

`RateClimb8Outer.h` reuses the frozen FFT stages but fuses the last butterfly
with output stores. It agrees with the dense oracle, including shortened tails.
This does not remove field arithmetic: outer dimensions 74, 80, and 84 all
depend on the same 321 GFNI / 769 XOR outer circuit per four rows.

With tight strides its benefit was inconsistent. With page-aligned strides
and fused outer loads, three further full-call launches gave the winning
ranges in the main table. The 84% profile improves by about another 2%.
Gather schedules still regress, so this is a recommendation for the measured
combined schedule, not a universal improvement to every outer path.

At K20, scatter routing lost to the regional non-temporal schedule. A folded
scatter launch gave about 6.57 ms for the 84% profile with fused outer loads,
versus approximately 4.94 ms for regional routing. This is a size-specific
choice, not a recommendation to replace the smaller-size scatter schedule.

## Correctness and measurement protocol

The independent scalar oracle directly evaluates the outer polynomials and
the full dense 16-by-16 packet map. It does not call the optimized FFT or
circulant circuit. Its retained-coordinate indexing is independently enumerated.
It checks the folded tables through explicit tile/coordinate/phase loops.

Release checks pass for:

- Every byte value in each packet input/state coordinate, both packet layouts,
  random cases, and both returned records and dual state.
- Natural requested sizes around one row, four rows, and a 64-row padding unit.
- Independent full forward/transpose bilinear identities on shortened inputs.
- Both outer schedules, both final-store choices, all routing schedules,
  compact and full-width inverse indices, unaligned offsets, and guard records.
- Full schedule agreement and exact input/output aliases at K16/K18/K20.

The 4 KiB-multiple strides also pass the K16/K18/K20 schedule and alias checks.
The final ASAN/UBSAN build passes both the independent small checks and the
full K16/K18/K20 replay, including all column-stride variants.
The K16/K18 checks establish implementation correctness, not new distance
certificates for those sizes.

Measurements use Peach's Ryzen 9 7950X, GCC 15.2, Release, CPU2, boost disabled.
All benchmarks are serial and guarded by a benchmark lock. Candidates share
payload and scratch addresses within each launch; candidate order rotates.
A fresh payload is restored outside timing before each batch of five calls.
Setup, allocation, and page preparation are excluded. Buffers at least 16 MiB
request huge pages. The whole call includes the inner, routing, outer, and
stores. Equal-map schedules are checked through matching output checksums.
Different codes are not expected to share checksums.
Timing uses a fixed reproducible setup; the margin column refers to the
uniform-setup ensemble covered by the frozen certificates.

## Reproduce

From the repository root on the measured ISA:

```text
cmake -S spin/experiments/char2_field_spin -B build/field -DCMAKE_BUILD_TYPE=Release
cmake --build build/field --target spin_char2_rate_climb8 -j2
build/field/spin_char2_rate_climb8 --check
build/field/spin_char2_rate_climb8 --large-check
taskset -c 2 build/field/spin_char2_rate_climb8 1048576 101 full auto page folded final-store
taskset -c 2 build/field/spin_char2_rate_climb8 1048576 51 phases auto page folded final-store
```

The final selectors choose routing, column layout, route folding, and the
outer final-store experiment. Run alternatives serially. The frozen 75%
control remains unchanged by those new-profile selectors.

The recommendation is the 84% profile with folded regional routing,
page-aligned column strides, fused outer loads, and fused final stores.
Retain the 84.09% profile as the highest-rate certified checkpoint. The next
useful engineering step is a smaller-size comparison of this same design,
with separate certificate replays before making K16/K18 distance claims.
Production promotion still requires separate user approval.
