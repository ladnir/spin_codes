# Smaller higher-rate outer: implementation and measurements

The GRS[128,87] research profile now has complete 10% distance certificates
at K16/K18/K20 with 63/79/67 bits of setup-failure margin. See the
[smaller-size certificates](FIELD128_SIZE_CERTIFICATES.md) and the preserved
[K20 certificate](FIELD128_K20_CERTIFICATE.md). Requested rates are 2/3,
2/3, and 128/189. No production profile changes.

The first implementation improved matched K20 time from 7.25--7.28 ms to
6.01--6.10 ms. The latest exact-map pass finds a modest compact-route win
in direct comparisons, but K20 process-level timings vary from 5.73 to
7.46 ms. Binary SPIN remains about 3.17 ms. Do not present the best launch
as a stable speedup; the latest results and rejected inner are recorded below.

## Construction and source layout

The inner and setup distributions are unchanged: width sixteen, eight
state symbols, the symmetric Cauchy phases, independent nonzero GF256
input/state diagonals, and one independently shuffled region per outer
coordinate. The outer changes from GRS[256,176] to systematic GRS[128,87].

The implementation accepts positive requested K within its uint32 routing
limit and pads to a multiple of 87*64. It computes exactly the requested
transposed prefix and preserves the remainder of an in-place N-record buffer.

Source coordinates are explicitly ordered in groups of four outer rows:

`source_index(row,j) = (row/4)*348 + 4*j + row%4`.

Thus each output vector already contains four adjacent source records.
The outer adds the systematic terms and stores this vector directly,
without the previous four-by-four record transpose. This is a fixed
source-coordinate permutation, not an additional randomized routing stage.
The full code is unchanged by this relabeling; arbitrary zero padding
restricts its message space and preserves the distance lower bound.
The scalar oracle uses the explicit new ordering rather than silently
comparing against the old row-major map.

## Kernels and validation

The [outer circuit](CHEAPER_OUTER_CIRCUIT.md) interpolates the source prefix
using 87=64+23 and 23=16+7. Its parity-only transpose maps 41 vectors to
87 vectors in 516 GFNI512 and 1188 XOR operations per four rows.
Systematic terms are added at the final stores. Caller and nested explicit
scratch total 10.5 KiB; compiler frames may add more. The dense packed
inner remains the fastest checked inner implementation.

Source lives in the isolated experiment:

- `spin/experiments/char2_field_spin/Cauchy128.h`: plan, route, dense inner,
  parity fusion, and four-row-interleaved output.
- `Cauchy128OuterFft.h` and its generator: independently checked fixed circuit.
- `cauchy128.cpp`: independent dense interpolation, scalar adjoints,
  correctness checks, serial phase timing, and matched public-API comparison.

The generator passes 3919 exact arithmetic/map checks and an independent
replay. Native forward, transpose, and parity-only outer hooks agree with
direct Lagrange interpolation. Full scalar tests use 1/4/8/16-byte records.
They check the field adjoint, size boundaries, padding, routing, all three
schedules, four alignments, guards, input preservation, aliases, and repeated
calls. Native K16/K18/K20 checks compare schedules and in-place suffixes.
All pass in Release and ASAN/UBSAN. Both registered Cauchy CTests pass.

The native record representation is a GF256-coordinate basis for the
compatible extension field. Conversion from another application basis is
not included; an integration requiring it must account for that extra cost.
There are no general large-field multiplications or timed heap allocations.

## Preserved initial matched timings

Peach, Ryzen 9 7950X, GCC 15.2, AVX512/GFNI, CPU 2. The comparison uses
the public PacketRs8T64 baseline and both higher-rate encoders. All receive
the same requested K, payload address, and initial bytes, with separate
workspaces. Five-call batches rotate candidate order after five warmups;
101 calls are timed per candidate. Restoration, setup, allocation, and
initialization are untimed. Benchmarks run serially.

Both field profiles use scatter at K16/K18 and regional non-temporal
routing/output at K20. All candidates use normal pages at K16 and prefer
huge pages at K18/K20. A page request is not an OS guarantee.

| Requested K | Seed | Binary SPIN (ms) | GRS256/176 (ms) | GRS128/87 (ms) |
|---:|---:|---:|---:|---:|
| 2^16 | 113 | 0.09918 | 0.27780 | 0.28822 |
| 2^16 | 7919 | 0.09845 | 0.27881 | 0.28780 |
| 2^18 | 113 | 0.46081 | 1.28689 | 1.12763 |
| 2^18 | 7919 | 0.46164 | 1.23559 | 1.17355 |
| 2^20 | 113 | 3.17621 | 7.28017 | 6.09705 |
| 2^20 | 7919 | 3.17104 | 7.24690 | 6.01213 |

The smaller outer regresses K16 by about 3--4%, improves K18 by about
5--12%, and improves K20 by 16--17%. K18 has substantial tail noise:
the new candidate's p90 is 1.59--1.70 ms. None of these times supports
promotion over binary SPIN on encoding speed alone.

At K16/K18 the new code pads to 66816/267264 source records, with
98304/393216 output records. At K20 it pads to 1052352 source records
and produces 1548288 outputs, 26.171875% fewer than half-rate SPIN.
The older higher-rate profile has 1540096 outputs, so the new K20 timing
does not benefit from a shorter output than its immediate reference.

Executable SHA256:
`ea8a427edb5213fd3722abf2804911d36fec14dff2f0d5b195ed441d53e7e5b6`.

Separate 51-call phase diagnostics, seed113, give:

| Requested K | Schedule | Full (ms) | Inner + route (ms) | Outer (ms) |
|---:|---|---:|---:|---:|
| 2^16 | scatter | 0.28798 | 0.14408 | 0.14397 |
| 2^18 | scatter | 1.13279 | 0.59048 | 0.54358 |
| 2^20 | region-nt-all | 5.62723 | 3.02103 | 2.60516 |

Phase diagnostics use a different allocation/cache history from the
matched comparison. Do not substitute their 5.63 ms result for the
6.01--6.10 ms matched result. K20 scatter costs 7.41 ms; the regional
schedule is important at that size. Medians of phases need not sum to
the median of complete calls.

## Reproduction and next decision

Use the configuration in [the preceding native checkpoint](FIELD256_NATIVE.md)
and build the new target:

```sh
cmake --build out/char2-field --target spin_char2_cauchy128 -j4
out/char2-field/spin_char2_cauchy128 --check
out/char2-field/spin_char2_cauchy128 --check-large
taskset -c 2 out/char2-field/spin_char2_cauchy128 --compare 65536 113 101 scatter full compact-batch
taskset -c 2 out/char2-field/spin_char2_cauchy128 --compare 262144 113 101 scatter full compact-batch
taskset -c 2 out/char2-field/spin_char2_cauchy128 --compare 1048576 113 101 region-nt-all full compact-batch
```

Repeat with seed 7919. For phase diagnostics, omit `--compare` and replace
`full` with `phases`. Run the two check modes in the separate sanitizer
build as described in the preceding checkpoint; do not time sanitized code.

The remaining gap is material in both stages. The next useful task is to
reduce inner-plus-route work or obtain a substantial packing/cache win,
not increase the already sufficient failure margin. Reduced-state and
cheap-mixing candidates must clear their known proof obstructions before
replacing this reference. The high-distance track remains open.

## Exact-map pass: compact indices and diagonal batches

Six compile-time variants now share the same plan and encoder map:
`baseline`, `compact`, `batch`, `compact-batch`, `lane-fft`, and
`compact-lane-fft`. There is one variant dispatch per encoding; inner
loops remain specialized. The default stays `baseline` for reproducibility.

The compact regional schedule stores local row indices in uint16 when
`L<=65536`. It avoids the global-offset subtraction and halves index
traffic. Larger regions fall back once per encode to the uint32 path.
Scatter uses the original route. The batch option loads sixteen input
diagonal bytes together and expands four coefficient vectors from them.
These changes leave the sampled code and its proof unchanged.

All controls retain both global route tables, plus the compact inverse
table where supported. At K20 this adds 3,096,576 setup bytes; total
experimental setup storage is 17,805,312 bytes. This is a measurement
control, not a claim of smaller total setup storage. No timed allocations
or application-basis conversions are introduced.

The first four-way executable had SHA256
`e63d901bf8a099957c0c1e16c0151c3b3713c64c4ad11c8af575162e8e5b9a3d`.
Its 101-call rotating, same-buffer K20 medians were:

| Seed | Baseline (ms) | Compact (ms) | Batch (ms) | Compact + batch (ms) |
|---:|---:|---:|---:|---:|
| 113 | 5.98187 | 5.73467 | 6.03479 | 5.68943 |
| 7919 | 6.01854 | 5.88966 | 6.12831 | 5.74112 |

That build's separate public three-way comparisons gave 5.85492 and
5.83749 ms for compact + batch. K16/K18 gains were below one percent.
These results are preserved, not substituted for the final build below.

## Rejected lane-local FFT inner

The new `Cauchy128InnerLaneFft.h` retains the certified residue coordinates.
It transposes each vector's four-by-four matrix of 32-bit chunks, so FFT
butterflies use lane-local shuffles. State remains in that representation
across packets and regions; its diagonal is transformed consistently.
This is a representation change, not a different state-randomization law.

The full packet needs 34 GFNI operations instead of the dense kernel's
70, but still needs 46 local shuffles and eight boundary permutations.
The reduced arithmetic count does not win. In the final six-way rotating
comparison, lane-FFT is about 5% slower at K16/K18 and also loses at K20.
The dense inner remains selected for the measured candidate.

All six variants pass independent scalar/full-map, alignment, boundary,
guard, alias, and K16/K18/K20 differential checks. New packet checks cover
every input/state byte-basis position, nonzero states, four phases, both
directions, and exact input/output aliasing. Release and ASAN/UBSAN pass;
both registered Cauchy CTests pass. A forced-empty compact table exercises
the uint32 fallback without a multi-million-record test allocation.

## Final-build timings and stopping point

Final executable SHA256:
`5b6aebec44b3c0fb48dbe6cb09e567cdc824bd612d948faa58346b7806898991`.
The six-way same-buffer diagnostic gives:

| K | Seed | Dense baseline (ms) | Compact + batch (ms) | Lane-FFT (ms) | Compact lane-FFT (ms) |
|---:|---:|---:|---:|---:|---:|
| 2^16 | 113 | 0.29099 | 0.28909 | 0.30641 | 0.30702 |
| 2^16 | 7919 | 0.28801 | 0.28637 | 0.30395 | 0.30374 |
| 2^18 | 113 | 1.12941 | 1.12901 | 1.18671 | 1.18771 |
| 2^18 | 7919 | 1.13057 | 1.13255 | 1.18869 | 1.18858 |
| 2^20 | 113 | 6.47963 | 6.34349 | 6.80189 | 6.66005 |
| 2^20 | 7919 | 5.87707 | 5.61440 | 6.17299 | 6.01993 |

`--variants` uses the diagnostic allocator, which requests huge pages
only for allocations at least 16 MiB. The public `--compare` flow instead
uses the common public memory policy described above. Both restore outside
timing and rotate five-call batches. They must not be pooled as one sample.

The final public three-way comparisons, with `compact-batch`, give:

| K | Seed / launch | Binary SPIN (ms) | GRS256/176 (ms) | GRS128/87 (ms) |
|---:|---|---:|---:|---:|
| 2^16 | 113 | 0.09873 | 0.27908 | 0.28715 |
| 2^16 | 7919 | 0.09925 | 0.28347 | 0.28846 |
| 2^18 | 113 | 0.46142 | 1.23074 | 1.12628 |
| 2^18 | 7919 | 0.46041 | 1.23768 | 1.12862 |
| 2^20 | 113 / first | 3.16768 | 7.25372 | 5.75997 |
| 2^20 | 7919 / first | 3.16597 | 7.25216 | 7.46337 |
| 2^20 | 7919 / repeat | 3.17046 | 7.27615 | 7.14506 |
| 2^20 | 113 / repeat | 3.19221 | 7.28443 | 5.73477 |

The slow K20 launches are real recorded results, not discarded outliers.
Their cause is unresolved. The same seed was faster in the four-way build
and in the isolated six-way diagnostic, so these observations do not by
themselves identify a bad route or a construction issue. They also do not
establish a dependable 5.7 ms production profile. Keep the previous checkpoint
and all controls; nothing is promoted.

The round ends here for the requested status discussion. The mathematical
gain is a certified higher-rate family at all three sizes. The implementation
still costs substantially more than binary SPIN. Further work should be chosen
after discussing that tradeoff, not justified by a larger failure margin.
