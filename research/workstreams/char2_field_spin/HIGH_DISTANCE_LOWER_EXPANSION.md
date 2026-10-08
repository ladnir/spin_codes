# Lower internal expansion for high-distance SPIN

This round prioritizes complete encoding cost over closing the earlier
8K internal-expansion control. K=2^20 and final rate 1/2 remain fixed.
Application records contain 128 bits; distance counts field symbols.
The frozen 3.2 ms high-rate winner and its certificate are unchanged.

## Best proof direction so far

GRS[256,64], four input symbols and four state symbols per step, and two
retained outputs gives the desired flow K -> 4K -> 2K. Both diagonals
use independent uniform nonzero GF(2^16) coefficients. The fixed outer
and Cauchy maps use embedded GF256 coefficients. There is one routing
stage: an independent uniform row permutation inside each of 256 regions.

Compared with the earlier GRS[64,8]/keep1 control, this halves the number
of inner steps. At K20, it also shrinks the region buffer from 2 MiB to
256 KiB. The larger outer needs more work per raw symbol, so neither
change alone establishes a speedup. The full prototype now measures about
15.6 ms, versus about 29 ms for the earlier 16-bit 8K control. This halves
the observed time but is not yet a compelling low-latency implementation.

At **40% distance**, the floating bounds give:

| Class | Sampled margin in bits |
|---|---:|
| Active rows A=1 / 16 / 64 / 128 | 1415 / 1205 / 828 / 873 |
| A=256 / 512 / 651 | 1397 / 3555 / 5272 |
| Middle raw density .03 / .10 / .30 / .60 | 1611 / 34102 / 270898 / 788917 |
| Middle .95 / .99 | 1724014 / 1470898 |
| Dense .99 / 1 | 1148441 / 1315647 |

These are individual point screens, not a whole-code certificate. The
minimum nonzero outer-row weight is 193. A raw weight below .03 times
4K therefore has at most 651 active rows. This suggests joins between
sparse, middle, and dense bounds, but coverage between samples and
directed replay are still required. Negative low-middle bounds at .003
and .01 need not be used if the sparse argument covers those words.

The new `theory/high_distance_expansion_screen.py` also supports periodic
retention, composing the actual chronological operators rather than
averaging their output counts. A 3K control with GRS[192,64] and retention
pattern (3,3,2) has positive tested points around 30%, but fails tested
bounds at 32%. Narrowing the processing width or shrinking state to one
or two symbols did not repair the tested 40% gaps.

We also tested direct-2K GRS[512,256] and [1024,512] over GF(2^16),
explicitly changing the fixed outer coefficient field. At 30%, their
sampled A=64 bounds are -13,268 and -25,240 bits; their middle .03 bounds
are -42,387 and -44,558 bits. These failed bounds do not establish actual
code obstructions. They give no reason yet to implement those larger
outer transforms. The helper records the outer field separately and
checks its embedding into the diagonal field.

## A real obstruction, not just proof slack

The simplest direct-2K construction, GRS[64,32] with four state elements,
cannot meet the desired high-distance/setup-margin target under uniform
regional routing. A space of messages can be supported inside short
windows in each region. If its dimension exceeds the number of state-reset
constraints, some nonzero message produces output only inside those windows.

An exact expectation argument gives failure probability at least
7.33e-5 at 35% distance, and at least .00615 at 40%. This applies regardless
of the diagonal alphabet. It is not an observed frequency or a failed
certificate search. The [proof and precise assumptions](HIGH_DISTANCE_WINDOW_OBSTRUCTION.md)
are independently audited and tested against exhaustive small examples.
This obstruction does not rule out the larger-outer 4K candidate.

## Full cost controls

`HighDistanceWide.h` and `high_distance_wide.cpp` implement precomputed
transposed cost controls, with genuine GF(2^16) or GF(2^32) diagonals,
the frozen explicit 16-packet inner schedule, inverse regional routing,
and an FFT outer. These are new research files, not production changes.
They are not fully tuned encoder implementations.

The first native launch measured these **GRS length-64 controls**, not
the promising length-256 outer:

| Outer / retention | Internal symbols | GF(2^16), ms | GF(2^32), ms |
|---|---:|---:|---:|
| [64,32] / keep4 | 2K | 6.69 | 8.22 |
| [64,16] / keep2 | 4K | 14.16 | 17.20 |
| [64,8] / keep1 | 8K | 28.79 | 36.59 |

The 8K concern is therefore substantiated by a complete implementation,
not just an operation estimate. These first controls use the generic
length-64 outer prefix; they do not establish an optimized latency floor.
They also run out of place with unchanged input, unlike the frozen
high-rate benchmark's repeated in-place calls.

For the 16-bit 4K control, a separate diagnostic run gives approximately
3.81 ms inner, 4.37 ms routing, and 5.86 ms outer. For the 8K control,
these phases take 6.81, 11.42, and 10.76 ms. Individual phase medians
need not sum to the full-call median. The extra work costs memory traffic
and outer evaluation as well as inner arithmetic.

Measurement host: Peach, Ryzen 9 7950X, GCC 15.2, CPU 2, AVX512/GFNI.
The initial table uses one launch, three warmups, and 31 timed calls per
configuration. All benchmarks run serially under the host lock. Setup,
allocation, and input generation are outside timing. The benchmark
allocator requests huge pages for payload buffers at least 16 MiB; the
OS need not grant the request. Precomputed coefficient arrays use ordinary
vector storage. Output checksums are consumed outside timing.

### Actual GRS[256,64] candidate

The new length-256 outer uses a pruned additive FFT. The forward map has
only its first 64 additive-basis coefficients, so its transpose combines
four length-64 cosets. It does not evaluate a dense matrix or perform a
systematic-basis conversion. All 16,384 point/coefficient values are
independently validated, including each of four packed record lanes.

Two serial launches per mode, each with 101 timed calls, give:

| Mode, GF(2^16), K20 | Median time across launches |
|---|---:|
| Out of place, unchanged input | 15.616–15.709 ms |
| Exact alias, no restore between calls | 15.621–15.632 ms |

The alias mode does not allocate the separate full-size output buffer.
Checksums agree between repeated launches of the same mode. A separate
31-call alias diagnostic gives 3.994 ms inner, 4.115 ms routing, and
7.539 ms outer. Alias allocation is therefore not the source of the
remaining cost. The smaller regional buffer helps the working-set size,
but does not make outer evaluation free.

The timing command waited for the shared benchmark lock before running;
it did not overlap the other task's benchmark. The two modes have the
same construction, but different input evolution between calls.

## Validation and reproduction

The length-64 controls pass independent scalar outer/inner comparisons,
extension-field dot-product adjoints, shortening, offset, guard, exact
alias, and full-width route-index checks. Release and ASAN/UBSAN with
leak detection pass. The new length-256 variant also passes its Release
and ASAN/UBSAN outer-basis and complete small scalar/adjoint checks.
The registered full CTest passes. The final Release binary SHA256 is
`cd573e3fc51145b88c7df4817b7e77ddf26176c3d7774305568b726099af7963`.
Twenty-eight
high-distance Python tests pass,
including exact window probabilities, chronological retention, old-model
agreement, and selected 40% witness replays. These checks do not prove
distance between sampled classes.

```text
python -B -m unittest discover -s research/workstreams/char2_field_spin/theory -p "high_distance*_tests.py" -v
python -B research/workstreams/char2_field_spin/theory/high_distance_fields.py --bits 16 --outer 256 64 --outputs 2 --distance 2/5 --active 1 16 64 128 256 512 651 --density .03 .1 .3 .6 .95 .99 --iterations 60
python -B research/workstreams/char2_field_spin/theory/high_distance_dense.py --bits 16 --outer 256 64 --outputs 2 --distance 2/5 --density .99 1
python -B research/workstreams/char2_field_spin/theory/high_distance_window_obstruction.py --distance 7/20 2/5
cmake --build build/char2-high-distance --target spin_char2_high_distance_wide -j 4
build/char2-high-distance/spin_char2_high_distance_wide --check
flock /tmp/spin-char2-benchmark.lock taskset -c 2 build/char2-high-distance/spin_char2_high_distance_wide --bench 16 2 31 full
flock /tmp/spin-char2-benchmark.lock taskset -c 2 build/char2-high-distance/spin_char2_high_distance_wide --bench 16 2 101 alias 256
```

Next test an exact-map outer-schedule improvement before changing routing.
The current outer materializes length-16 leaves, length-32/64 stages, and
coset accumulation through scratch arrays. Fuse those stages into the
accumulator and final stores, following the frozen fused-outer design.
Then test a tiled plane layout separately. The length-256 outer performs
about 2.89 million GFNI instructions at K20, versus 2.11 million for the
length-64/dimension-16 control: fewer rows did not mean less arithmetic.
Retain the
4K/40% construction as the best current proof candidate, not a promoted
performance winner. Revisit its certificate once the end-to-end cost is
acceptable. Do not return to an 8K certificate campaign without a concrete
cost advantage.
