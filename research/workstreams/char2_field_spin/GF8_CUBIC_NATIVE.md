# Cubic width-eight encoder: native checkpoint

The certified GRS[128,83] / width8-state4 construction has a validated native
transposed encoder. The first packed implementation was slower than the old
field control. An exact-map two-packet schedule removes an artificial state
dependency and cuts its time by about 14--16%, reaching .237 / .94 / 5.71--5.86
ms at K16/K18/K20. It still has a lower rate than the previous field profile
and remains slower than binary half-rate SPIN. Keep it as a research control,
not an automatic production replacement.

The accompanying [certificate](GF8_PACKET_CERTIFICATE.md) covers K20 at
10% symbol distance and a 65-bit whole-code setup-failure margin. This
implementation round establishes no additional size or rate certificate.
The next construction target remains at least 3/4 rate at 10% distance.

## Construction and implementation

Each 128-bit record contains sixteen GF256 lanes. The scalar code uses an
MDS[128,83] outer in additive-polynomial coefficient coordinates, independent
uniform regional permutations in the proof ensemble, and the repaired
cubic maps A=(u+v)^3+1, C=(u+v)^3+3 with eta=2. Independent nonzero input
and state diagonals are fixed at setup. The implementation's seeded setup
is the research test fixture used by the other native controls.

The first 83 coefficients of the full length-128 additive FFT span exactly
the required outer code. Four source rows are interleaved, with source
index `(row/4)*332 + 4*coefficient + row%4`. Padding fixes the unused source
suffix to zero. There is one regional routing stage, and the four-record
state persists across regions. Encoding uses caller-owned scratch and
writes exactly the requested transposed prefix, preserving the input tail.

The paired map uses nine GFNI512 instructions and sixteen XOR/XOR3
instructions per eight records, including the diagonals and eta scaling.
It also uses two record transposes, at most sixteen lane shuffles. These
counts exclude loads, stores, coefficient expansion, and spills. In
particular, a lower multiplication count does not bound recurrence latency.

Two exact-map outer schedules are retained: gather then transform, and
source loads fused into FFT leaves. Both use the same route and matrices.

## First matched serial measurements

Host: Peach, Ryzen 7950X, boost disabled, one thread pinned to CPU2.
Compiler: GCC15.2, Release, AVX512/GFNI, `-mtune=znver4`, no LTO.
Each row uses 101 calls, five warmups, rotating batches of five calls, and
seeds 113 and 7919 in separate launches. The field controls share one input
address and scratch buffers. Setup, allocation, and input restoration are
outside timing; the entire encoder, including routing and output stores,
is inside. Public SPIN supplies the binary control and allocation policy.
K16 uses normal pages; K18/K20 request huge pages where supported.
No benchmarks ran concurrently.

Median times in milliseconds, showing both launch results:

| Encoder | K16 | K18 | K20 |
|---|---:|---:|---:|
| Previous [128,87], systematic outer | .3130 / .3158 | 1.2518 / 1.2480 | 6.0592 / 7.8048 |
| Previous [128,87], coefficient outer | .2723 / .2745 | 1.0870 / 1.0851 | 6.0406 / 6.2542 |
| New [128,83], gather outer | .2788 / .2795 | 1.1395 / 1.1364 | 6.8457 / 6.8749 |
| New [128,83], fused outer | .2770 / .2776 | 1.0923 / 1.0845 | 7.4878 / 7.4383 |
| Binary RS8 public library, rate 1/2 | .09846 / .09851 | .46233 / .46290 | 3.1940 / 3.1847 |

The previous coefficient control uses fused loads at K16/K18 and gather
at K20, as selected in its earlier experiments. These are different codes
with different rates, not equal-output-size comparisons.

| New profile | Requested K | Padded K | Output N | Requested K/N |
|---|---:|---:|---:|---:|
| K16 | 65536 | 69056 | 106496 | 61.538% |
| K18 | 262144 | 265600 | 409600 | 64.000% |
| K20 | 1048576 | 1051776 | 1622016 | 64.646% |

The nominal unshortened rate is 83/128=64.844%. Only the last row currently
has the complete certificate for this inner. The previous field profile
emits 98304 / 393216 / 1548288 records. Thus the new profile emits
8.333% / 4.167% / 4.762% more records on these requests.

The systematic control's K20 launch variation is retained, not discarded.
It also occurred in earlier experiments and remains unexplained. This
build's small-size controls are slower than their earlier separate
executable; the table is the current matched comparison, not an overwrite
of that earlier receipt.

## Phase evidence

A separate seed113 phase run gives these median milliseconds:

| Profile and phase | K16 | K20 |
|---|---:|---:|
| Previous coefficient: inner plus route | .1706 | 3.3470 |
| New cubic: inner plus route | .1972 | 4.1894 |
| Previous coefficient: outer | .1015 | 2.8487 |
| New cubic: outer (fused K16, gather K20) | .0794 | 2.4924 |

Phase medians need not sum exactly to the median total. The observed
regression is in the inner-plus-route phase, not the outer. The new inner
performs twice as many state updates per output count, with more record
shuffles per output. This is a plausible explanation, not an isolated
measurement of shuffle or recurrence cost. No multiplication-throughput
roofline claim follows from these counts.

## Exact-map paired schedule

The original packed SIMD circuit places input-only A-transpose moments and
state-dependent C-transpose emissions in different lanes of the same ZMM
instructions. This creates a register dependency on the incoming state for
work that is mathematically independent of it.

The paired variant computes the input moments of two consecutive reverse
packets together, advances their two state recurrences, and then batches
their emissions. It retains the original packet order, diagonals, route,
maps, and outputs. No new setup, global route, or scratch allocation is
needed. The arithmetic counts are unchanged; the loop-carried state path
no longer includes the packed circulant circuit.

The final `paired` comparison uses the same outer for both cubic variants:
fused at K16/K18 and gather at K20. It retains five candidates and the same
serial benchmark protocol. Median milliseconds for seeds113 / 7919:

| Encoder | K16 | K18 | K20 |
|---|---:|---:|---:|
| Previous [128,87], systematic | .3151 / .3124 | 1.2897 / 1.2434 | 5.7831 / 7.5632 |
| Previous [128,87], coefficient | .2720 / .2688 | 1.1373 / 1.0773 | 6.0816 / 6.2335 |
| Cubic8, original packing | .2787 / .2767 | 1.0891 / 1.0952 | 6.7676 / 6.9630 |
| Cubic8, paired schedule | .2372 / .2364 | .9396 / .9438 | 5.7123 / 5.8587 |
| Binary RS8 public control | .09871 / .09915 | .46007 / .46145 | 3.1854 / 3.1807 |

The exact-map reduction is about14--16%. In a separate K20 phase run,
inner plus routing falls from4.1907 to3.1343ms; outer times remain close
at2.7614 and2.7392ms. This supports the scheduling diagnosis. It does not
separately measure shuffle throughput or establish a machine roofline.

The paired path passes all original scalar, alias, schedule, and sanitizer
checks. Its output checksums match the original in every timing launch.
The final executable SHA256 is
`c910a5372d740d75e0a4fd15b0eed33637b760e44b3f5638ef443faba820ab86`.

## Validation and reproduction

Release and AddressSanitizer/UndefinedBehaviorSanitizer both pass:

- exact dense/scalar packet comparison across phases, basis vectors,
  and all byte values;
- direct GRS[128,83] evaluation and strided transpose references;
- complete scalar transpose and shortened forward/transpose identities;
- gather/fused schedules, offsets 0..3, exact aliasing, untouched tails,
  guards, and uint32 route fallback;
- K16/K18/K20 schedule and alias comparisons.

Configure the research CMake project with optional `SPIN_CHAR2_PUBLIC_ROOT`
pointing at the library source, then build `spin_char2_cubic8`.

```sh
./spin_char2_cubic8 --check
./spin_char2_cubic8 --large-check
taskset -c 2 ./spin_char2_cubic8 1048576 113 101 full all
taskset -c 2 ./spin_char2_cubic8 1048576 7919 101 full all
taskset -c 2 ./spin_char2_cubic8 1048576 113 101 phases all
taskset -c 2 ./spin_char2_cubic8 1048576 113 101 full paired
taskset -c 2 ./spin_char2_cubic8 1048576 7919 101 full paired
taskset -c 2 ./spin_char2_cubic8 1048576 113 101 phases paired
```

Repeat the timing commands with K=65536 and 262144, serially. The first-round
executable SHA256 was
`ee6b50951c44f6a9ffb063ec3e8e3e35f3a30b4eec9146ba3fe477169857ac13`.
Its full-call gather/fused output checksums agree within each launch.

The implementation lives in `spin/experiments/char2_field_spin` and remains
separate from production. No certificate, production profile, or paper
performance number was replaced. The useful next step is a higher-rate
construction/proof change, retaining paired scheduling as an implementation
lesson. See [larger-outer screens](GF8_HIGH_RATE_OUTER.md) and the
[inner/return analysis](GF8_HIGH_RATE_INNER.md).
