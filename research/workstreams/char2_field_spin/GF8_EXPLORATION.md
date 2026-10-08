# GF256 exploration: cheaper encoding at 10% distance

Subsequent work: the [fixed-inner outer-rate sweep](GF8_OUTER_RATE_SWEEP.md)
closes a separate dimension88 K18 certificate. Its size improvement is
not included in the dimension87 timings recorded in this note.

The goal is a useful higher-rate GF256 code on 128-bit records, with one
regional routing stage and competitive full encoding cost. Each record
contains sixteen components sharing the same scalar code. The existing
binary half-rate/10% implementation remains the performance baseline.
This work does not replace production or the completed certificates.

## Proof work: retain state rank

The [rank-history transfer](GF8_RANK_HISTORY_SCREEN.md) and its
[independent audit](GF8_RANK_HISTORY_AUDIT.md) retain the incoming state-form
rank instead of maximizing it afresh in every packet. The generalized
positive transfer preserves the required continuation cone. Its first
K20 GRS[128,86], six-state screens improve the .10 and .14 input-density
classes by about 335 and 389 bits, but leave negative margins.

An independent exact enumeration of the six-state return map
I+2AC proves output-weight minima 15 and 13 for one and two active inputs.
The earlier generic minima were 10 and 9. Despite this large local
improvement, inserting the exact minima helps the .10 bound by only 1.4 bits.
This identifies the continuation/source part of the middle-class bound,
not sparse returns, as the current numerical bottleneck.

The completed source check removes impossible positive outer weights
below the MDS minimum and uses the sharper modal-conditioning factor.
Homogenizing the truncated row polynomial needs an additional condition:
beta=tau*lambda>=2 makes its discrete second differences nonnegative.
Outside that range, retain the old bound. Do not insert the truncated
homogeneous polynomial without this justification.

Modal conditioning removes another 827.27 bits of loss. Feasible source
supports remove less than .001 additional bits: the tilted occupied-row
term already favors weights well above the MDS minimum. The final sampled
margins are -17150.56 bits at density .10 and +123.88 bits at density .14.
These floating point screens are not whole-code certificates. The .10
class still blocks the six-state candidate; no native implementation of
that unclosed candidate is justified by these results. Thirteen transfer/cost
tests and three independent return-geometry tests pass.

## Implementation opportunity: remove systematic interpolation

The certified outer is the set of evaluations of degree-at-most86
polynomials at the128 GF256 labels0 through127. The old interface takes
the first87 evaluations as source coordinates. Instead, take coefficients
in the normalized additive-polynomial basis of degrees0 through86.

Let E_s be the old128-by87 systematic generator. Let E_c evaluate the
coefficient basis at all128 labels, and let B be its first87 rows.
These basis polynomials have distinct degrees and nonzero leading
coefficients. They form a basis of the polynomials of degree below87.
Evaluation at87 distinct points is invertible, so B is invertible and

    E_c = E_s B.

For each fixed route and each fixed set of inner/input diagonals, the two
full padded encoders therefore have exactly the same codeword image.
All existing whole-code distance certificates transfer pointwise, with
no probability loss. Requested-prefix zero padding can select different
subcodes in the two bases, but both remain subcodes of this same certified
full code. The K16/K18/K20 rate, distance, and margin claims are unchanged.

The transposed output changes by B^T on each source row. It is not
byte-identical to the systematic output, and the change must not be hidden
from callers who depend on the original source coordinates. This is a new
research interface, not a transparent replacement of an existing API.

The [generated outer](../../../spin/experiments/char2_field_spin/Cauchy128OuterCoefficients.h)
uses 321 GFNI and 769 XOR operations per four rows, versus 516 and 1275 for
the complete old outer transpose. It uses no lane shuffles. The
[generator](../../../spin/experiments/char2_field_spin/generate_cauchy128_outer_coefficients.py)
passes 4147 exact checks, including direct polynomial evaluation, an
independent interpretation of the emitted C++ statements, adjoints, and
the generator identity above. These counts do not establish a timing win.

The first native implementation passes Release and ASAN/UBSAN checks,
including the direct polynomial oracle, basis inputs, SIMD lanes, adjoints,
alignment, aliases, output guards, requested tails, and the full padded
old-to-new basis transform. Large checks cover K16/K18/K20. The inner,
setup distribution, routing, and source-record ordering stay unchanged.
The full comparison includes all 128 outer-input loads, changed output
coordinates, scratch, padding, and final stores.

## First matched full-encoder measurements

Peach Ryzen 7950X, CPU2, GCC15.2, Release -O3, AVX512/GFNI,
`-mtune=znver4 -fno-lto`. All runs are serial. Each launch uses five warmups
and 101 timed calls, alternating implementations in rotating five-call
batches. The candidates share input, route, planes and scratch addresses.
Input restoration, setup and allocation are outside timing. Ordinary pages
are used at K16; K18/K20 request huge pages. The binary reference is the
existing public half-rate encoder in the same executable.

The ranges contain medians from seeds 113 and 7919, not confidence intervals.
Every element holds 128 payload bits. Timings are milliseconds.

| Requested K | Output N | Systematic outer | Coefficient outer | Binary reference |
|---:|---:|---:|---:|---:|
| 2^16 | 98304 | .287--.289 | .259--.260 | .0983--.0986 |
| 2^18 | 393216 | 1.125--1.138 | 1.043--1.045 | .463--.470 |
| 2^20 | 1548288 | 5.690--5.750 | 6.058--6.069 | 3.149--3.173 |

The source-basis change saves about 10% at K16 and 7--8% at K18, but loses
5.6--6.5% at K20. All three retain 10% distance and their existing
63/79/67-bit certificates. The requested rates remain 2/3, 2/3, and 128/189.
This is a gain within the higher-rate family, not a speed win over binary
SPIN. No old interface is silently changed.

A separate K20 phase diagnostic gives 3.151/2.532 ms for the systematic
inner-route/outer and 3.143/2.920 ms for the coefficient version. The
regression is in the outer, not routing or inner computation. Fewer field
products do not guarantee lower runtime: the first coefficient wrapper
gathers all 128 planes into scratch before its first transform stage.
A subsequent fused-load control tests that scratch traffic directly.

Initial timed Release executable SHA256:
`158641f011bebcba0000e54245aee43693352d0ac3f5af464577dc65a310693d`.
Its registered coefficient-encoder CTest passes.

### Fused source loads: select only at K16/K18

The additional exact-map control loads sixteen planes directly into the
registers for the first transform stage. It removes 128 scratch stores
and 128 scratch reloads per four-row tile (16 KiB), without changing the
321 GFNI/769 XOR arithmetic count. Higher transform stages are unchanged.

The final executable compares all three field paths and the binary control
at the same addresses and with the same method as above. Two independent
launches use seeds 113 and 7919. Timings are milliseconds.

| Requested K | Systematic | Coefficient gather | Coefficient fused | Binary reference |
|---:|---:|---:|---:|---:|
| 2^16 | .2860--.2866 | .2573--.2579 | .2430--.2437 | .09845--.09857 |
| 2^18 | 1.1254--1.1318 | 1.0373--1.0405 | .9555--.9557 | .4626--.4628 |
| 2^20 | 5.7582--5.8294 | 6.0631--6.0830 | 6.7741--6.8736 | 3.1780--3.1823 |

The fused path reduces full encoding time by 14.8--15.2% at K16 and
15.1--15.6% at K18 relative to the systematic path. At K20 it is worse than
both other field paths. A final phase diagnostic measures inner-route
3.131/3.117/3.102 ms and outer 2.663/2.951/3.625 ms for systematic/gather/fused.
Thus the larger-size regression remains in the outer. Reduced local traffic
alone does not compensate for its source-load schedule at this working set;
the exact cache/latency mechanism is not established by these timings.

Keep the fused coefficient interface as the K16/K18 research winner and
retain the systematic K20 reference. They are different source bases, not
interchangeable byte-level implementations of one fixed interface. No
runtime size-dependent basis switch or production promotion is introduced.

Final Release executable SHA256:
`0830e2faa5894a920afb6160f9c03eb364db2dda211f267277daaf5d619e5f84`.
Both coefficient paths pass Release, ASAN/UBSAN, registered CTest and large
K16/K18/K20 checks. Strided tests additionally cover strides 64, 83 and 192
bytes, byte-unaligned inputs, source preservation, and output guards.
The final generator performs 4147 exact checks. The matched coefficient
checksums agree; systematic checksums deliberately differ.

```sh
python -B spin/experiments/char2_field_spin/generate_cauchy128_outer_coefficients.py --check
cmake -S spin/experiments/char2_field_spin -B build-gf8 -DCMAKE_BUILD_TYPE=Release -DSPIN_CHAR2_PUBLIC_ROOT=/absolute/path/to/spin
cmake --build build-gf8 --target spin_char2_cauchy128_coefficients
build-gf8/spin_char2_cauchy128_coefficients --check
build-gf8/spin_char2_cauchy128_coefficients --large-check
taskset -c 2 build-gf8/spin_char2_cauchy128_coefficients 65536 113 101 full all
taskset -c 2 build-gf8/spin_char2_cauchy128_coefficients 262144 113 101 full all
taskset -c 2 build-gf8/spin_char2_cauchy128_coefficients 1048576 113 101 full all
```

Repeat serially with seed 7919. `phases` replaces `full` for the diagnostic;
`gather` or `fused` replaces `all` to preserve either paired comparison.

## Lower-rate six-state check

A bounded follow-up changes the outer dimension from 86 to 80, for padded
rate 5/8. At K20, all six sampled densities .06/.08/.10/.12/.14/.18 have
positive margins under the refined six-state transfer. This is evidence
for a plausible lower-rate proof, not coverage of all weights or rows.

The cost gate is less encouraging. Six state symbols still occupy two
packed SIMD vectors. With that layout, the dense inner drops from 70 to
62 GFNI operations per packet, not to 50. The audited combined systematic
outer/inner count falls only 4.84% on the same requested K20, while the
output length and route traffic increase 8.47%. No full encoder timing
for this construction is claimed. The [rank-history note](GF8_RANK_HISTORY_SCREEN.md)
records the sampled values and cost derivation. A new native or full-proof
campaign needs a stronger cost case first.

## Completed exploration milestone and next step

1. Rank-history transfer generalized and independently audited; sampled
   six-state gap remains. Do not start directed certification from these points.
2. Sparse return geometry sharpened exactly; its measured effect on the
   failing bound is small.
3. Source-support truncation is valid in its stated domain but contributes
   negligible numerical improvement. Sharper conditioning helps, without
   closing the six-state candidate.
4. Coefficient-basis outer preserves the proved code image and removes
   interpolation work. Fused loads bring the K16/K18 full-encoding gain
   to about 15%, but both new schedules lose at K20. Preserve that negative
   result rather than promoting a uniform replacement.

The milestone delivers a checked implementation improvement at K16/K18
and preserves the negative K20 and six-state findings. It does not meet
the broader goal of beating binary encoding cost at higher rate. Next,
target an inner geometry or packing that reduces actual SIMD work;
reducing scalar state dimension alone is not enough here.
All native benchmarks ran serially. No commit, push, deletion, or promotion
was part of this exploration request.
