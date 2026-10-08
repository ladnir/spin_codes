# Where the remaining GF256 encoding time goes

The [3.2 ms routed-sixteen implementation](GF8_SPEED_ROUTED.md) remains the
selected K20 checkpoint. Sixteen exact-map schedules were compared in this
round. None supplies a material, repeatable improvement over that checkpoint.
No slower experiment was promoted, and no production code changed.

The construction remains GRS[64,40] with four input and state coordinates,
three retained outputs, 83.3333% nominal rate, and 10% symbol distance.
The existing K20 certificate still gives 57-bit whole-code setup-failure
margin. These experiments change execution, not the sampled code.

## Measured decomposition

The diagnostic driver measures each region's computation followed by its
routing. It carries the actual recursive state across all regions. Thus the
route sees the producer's temporary in its normal cache state; it is not
measured on unrelated prefilled scratch. The outer then runs normally.

| Phase | Representative median | Approximate share |
|---|---:|---:|
| Inner computation, including temporary stores | 1.258 ms | 39% |
| Routing, input diagonal, and plane stores | 1.161 ms | 36% |
| Outer and final output stores | 0.843 ms | 26% |

The instrumented complete call measured 3.260 ms. Medians need not add
exactly, and the displayed shares are rounded. Clock calls are compiled out
of the ordinary benchmark variants. Their complete call times decide the
winner; phase timings guide diagnosis rather than establish a speedup.

The important result is the balance: routing is almost as expensive as the
inner arithmetic. Optimizing the recurrence alone leaves substantial work.

## Final matched comparison

All rows below come from the final 101-call launches, after correcting an
outlined helper in the two-buffer experiment. Every candidate sees the same
input and scratch addresses. Candidate order rotates in batches of five.
The first protocol restores input before each batch; the second restores
input before every call. Restoration and setup are outside the timer.

| Schedule | Full median | Restore-every-call median |
|---|---:|---:|
| Selected routed-sixteen | **3.205 ms** | **3.380 ms** |
| Equivalent profiler control, clocks disabled | 3.205 ms | 3.385 ms |
| Balanced assembly of four loaded records | 3.212 ms | 3.390 ms |
| Dword coefficient broadcast | 3.226 ms | 3.405 ms |
| Compact constant-map circuits | 3.238 ms | 3.433 ms |
| Continuous two-tile inner pipeline | 3.221 ms | 3.401 ms |
| Compact maps plus continuous pipeline | 3.200 ms | 3.382 ms |
| Precomputed 32-bit byte offsets | 3.562 ms | 3.736 ms |
| Precomputed 16-bit qword offsets | 3.322 ms | 3.420 ms |
| Explicit 128-bit record loads | 3.398 ms | 3.597 ms |
| Two-region producer/route pipeline | 5.538 ms | 5.918 ms |

The compact/continuous combination is effectively tied with the selected
implementation. Its 0.005 ms repeated-call difference does not justify a
new recommendation; it loses that difference under input restoration.
The remaining variants combine balanced record assembly, coefficient
broadcasts, and these alternatives. None changes that conclusion.

Measurements use Peach's Ryzen 9 7950X, GCC 15.2, AVX512/GFNI, CPU 2,
and disabled boost. Large benchmark buffers request huge pages. Benchmarks
were strictly serial under the existing host lock. Absolute timings vary
slightly across launches; comparisons use matched controls within a launch.

## Why the alternatives were worth testing

The compact inner exploits the missing fourth emission input and normalized
map identities. Its source circuit removes 24 Boolean operations per
sixteen packets while retaining 72 GFNI sites. This count did not predict
a native speedup. The continuous pipeline also retains two live tiles;
extending lookahead across loop boundaries was essentially neutral.

Byte offsets remove index scaling but double the compact route metadata.
Qword offsets instead retain 16-bit entries while permitting x86's
8-byte addressing scale. They are available when a region has at most
32,768 records; larger regions use the existing route representation.
Neither improved the complete encoder. The explicit-load control also
showed that fewer source-level loads do not guarantee a faster gather.

The two-region experiment uses two 410 KiB temporaries at K20. While it
computes a 64-record interval of the next region, it routes a 64-record
interval of the previous completed region. The recursive chronology and
sampled map are unchanged. The larger active streams raise a cache-capacity
concern on this host's 1 MiB L2, but no hardware-counter experiment isolates
that contribution here.

An assembly check first found a separate implementation defect: GCC outlined
the route's unroll wrapper. The resulting helper repeatedly loaded captured
pointers through a closure. Sixteen explicit forced-inline calls removed
that avoidable overhead. The rebuilt compact- and wide-route kernels have
no calls or vector stack spills. The table reports the corrected implementation;
it still loses substantially. The original slowdown must not be attributed
entirely to the algorithm or cache capacity.

## Validation and retained artifacts

Release and ASAN/UBSAN checks pass for all sixteen alternatives, with leak
detection enabled. Small cases use the independent dense-Cauchy and direct
polynomial oracle. Large cases compare the exact sampled map at K16, K18,
and K20. The checks include dimensions 39 and 40, shortened tails, all four
record-alignment offsets, input preservation, exact aliases, unchanged
suffixes, guards, the instrumented schedule, and forced wide-route fallback.
Additional local checks compare arbitrary incoming states and final states
for the compact and continuous inner kernels at several region lengths.

The selected headers `RateSpeed4RouteNext.h`, `RateSpeed4InnerNext.h`, and
`RateSpeed4OuterNext.h` are unchanged from the previous checkpoint.
New experimental code is isolated in `RateSpeed4Profile.h`,
`RateSpeed4InnerMore.h`, `RateSpeed4Pipeline.h`, and `rate_speed4_profile.cpp`.

```text
cmake -S spin/experiments/char2_field_spin -B build/field -DCMAKE_BUILD_TYPE=Release
cmake --build build/field --target spin_char2_rate_speed4_profile -j2
build/field/spin_char2_rate_speed4_profile --check
build/field/spin_char2_rate_speed4_profile --large-check
taskset -c 2 build/field/spin_char2_rate_speed4_profile 1048576 101 phases
taskset -c 2 build/field/spin_char2_rate_speed4_profile 1048576 101 full
taskset -c 2 build/field/spin_char2_rate_speed4_profile 1048576 101 fresh
```

Keep the selected implementation as the control. A larger improvement now
needs a credible reduction in actual arithmetic or random-access work,
not merely fewer source statements or a more complicated execution schedule.
This is an empirical conclusion about the tested variants, not an optimality
claim for the construction.

## Next construction experiment: route pairs

A concrete next candidate is to permute pairs of adjacent outer rows within
each region. Each route entry would address two contiguous 128-bit records,
allowing a 256-bit load and halving the number of indices. Independent input
and state diagonals would remain. This changes the permutation distribution;
neither the current certificate nor a speedup has been established for it.

The local inner bounds already allow related input forms. The new difficulty
is correlated support placement: two equal outer words remain active together
in a pair. Independent nonzero diagonal coefficients randomize their values
but do not separate their supports. Consequently, the scalar occupancy and
independent-support bounds cannot simply be reused.

Before building a full certificate, test the sparse bounds for pairs with
both rows active, starting with total active-row counts 2, 32, 64, and 96.
Then include singly active pairs and derive a joint two-row source polynomial
for the middle and dense ranges. Dimension 39, giving 81.25% nominal rate,
is a fallback within the target rate range if dimension 40 loses too much
margin. A serial implementation comparison can independently determine
whether the contiguous loads justify this proof work.
