# High-distance exploration: first larger-field results

Date: 2026-10-05. Requested dimension K=2^20; final rate 1/2;
128-bit application records; one regional routing stage. The eventual
acceptance target is at least 40 bits of whole-code setup-failure margin.
The [3.2 ms high-rate winner](GF8_HIGH_RATE_WINNER.md) remains frozen.

## Outcome and next decision

Two candidates pass the sparse, middle, and dense **points tested so far**:
42% symbol distance over GF(2^16), and 45% over GF(2^32). Neither has a
complete distance certificate. Coverage between samples, joins between
the three bounds, and directed numerical replay remain open.

The arithmetic results make both worth pursuing. Genuine GF(2^16)
diagonal multiplication is close to byte multiplication in the measured
pass. A precomputed GF(2^32) multiplication matrix is about 2.7 times
faster than the old tower kernel at K16/K18, and reduces its K20 pass time
by about 25%. These are kernel measurements, not full encoder timings.

The main cost concern is the construction: the current positive control
processes eight raw symbols per message symbol. First try to close its
45% certificate, retaining the cheaper 42% control. Then reduce the raw
expansion and measure the complete encoder. Do not claim a fast high-distance
encoder from the multiplier timings alone. The tested 49% bounds still
have a substantial middle-class gap; increasing the alphabet alone does
not resolve it.

## Concrete construction being screened

Let F=GF(2^m), for m=16 or 32, with an embedded AES GF256. Encode each
eight-symbol message row with the fixed GRS[64,8] outer over F. The fixed
outer coefficients remain in GF256. Route each of its 64 coordinates
to its own region and apply an independent uniform permutation of the rows
within every region.
Process the resulting stream with the frozen four-input/four-state
full-Cauchy inner, now over F. Retain only output coordinate zero from
each four-input step. State continues across region boundaries.

Input and state diagonal coefficients are mutually independent uniform
nonzero elements of **the full field F**, sampled at setup. The Cauchy
maps remain fixed over the embedded GF256. This is a genuinely different
sampled code, not scalar extension of a sampled GF256 code. A restricted
coefficient family would need its own cancellation and membership bounds.

| Quantity at K20 | Value |
|---|---:|
| Message symbols K | 1,048,576 |
| Outer rows | 131,072 |
| Raw inner inputs | 8,388,608 = 8K |
| Retained outputs | 2,097,152 = 2K |
| Final rate | 1/2 |
| Width / state / retained outputs per step | 4 / 4 / 1 |
| Bytes in one full region of 128-bit records | 2 MiB |

The last quantity exceeds the measurement host's 1 MiB per-core L2.
The frozen high-rate winner uses about 410 KiB per region. Routing and
scratch traffic therefore need new measurements. No complete native
encoder for these high-distance candidates has been implemented or timed.

A 128-bit record holds eight parallel GF(2^16) lanes or four parallel
GF(2^32) lanes, with the same scalar code applied to all lanes. Distance
here counts nonzero field symbols; it is not a bit-distance claim.

## Proof evidence and its limits

The new `theory/high_distance_fields.py` applies the existing full-Cauchy
rank-potential and regional bounds with the actual field cardinality.
Its sparse classes use active outer-row count A. Its middle classes use
raw-input weight divided by raw length. The values below are the negative
base-two logarithms of the floating upper bounds for individual sampled
classes; they are not whole-code margins.

| Candidate | Sparse A=1 | A=16 | A=64 | A=128 | A=256 |
|---|---:|---:|---:|---:|---:|
| GF(2^16), distance 42% | 446.43 | 798.62 | 1,530.49 | 2,800.18 | 5,524.75 |
| GF(2^32), distance 45% | 851.89 | 681.88 | 1,148.65 | 2,044.43 | 4,019.95 |

At A=512, the GF(2^32) bound gives 8,282.61 bits using the A=256
witness parameters without reoptimizing them.

| Candidate | Middle density .003 | .01 | .03 | .95 | .99 | .995 |
|---|---:|---:|---:|---:|---:|---:|
| GF(2^16), 42% | 867 | 26,904 | 108,568 | 1,426,227 | 782,079 | not sampled |
| GF(2^32), 45% | 1,678 | 24,180 | 79,626 | 2,000,984 | 1,568,173 | 1,353,472 |

These rounded results come from numerical witness searches. Optimizer
convergence is not a proof of optimality; any valid evaluated witness
would suffice for an upper bound once rigorously replayed.

The separate `theory/high_distance_dense.py` uses exact MDS support counts
and the established positive-excess/Maclaurin bound, evaluated in floating
logarithms. It includes the output-ball bound and the completion factor
q^(raw length - retained length) for punctured outputs. Density one is
handled here, not by the middle routine.

| Candidate | Dense density .99 | .995 | .999 | 1 |
|---|---:|---:|---:|---:|
| GF(2^16), 42% | 117,094 | 412,162 | not sampled | 625,028 |
| GF(2^32), 45% | -409,302 | 472,965 | 1,131,675 | 1,272,356 |

The GF(2^32) middle and dense bounds are both positive at .995, providing
a candidate join location. This is not interval coverage. In particular,
untested active-row counts and the sparse-to-middle join remain to be checked.

### What did not solve the problem

- Alphabet alone: for GRS[64,16], width 4/state 4, keep 2, distance 36%, and
  A=64, the same witness gives 149.316 bits over GF256 and 149.469 bits
  over either GF(2^16) or GF(2^32). Its effective rank tilt is about 15,
  below the byte-field limit. More coefficient values add almost nothing.
- More state alone: width 4/state 8, keep 2 at 40% still has negative sampled
  sparse and low-middle bounds. State 16 at 45% improves A=64 but leaves
  negative bounds at A=128 and 256. These are bound failures, not codewords.
- A minimum-source-weight Chernoff mark enforcing raw weight at least
  49A for the tested GRS[64,16] controls chose its neutral value in all
  three cases. This particular refinement gave no
  gain; it does not rule out a stronger row-spectrum argument.
- GF(2^128) at 49%, using the new [64,8]/keep 1 geometry, has positive
  sampled sparse points and dense points sufficiently near density one.
  However, middle densities .95 and .98 give about -595,047 and -887,649
  bits. Expanding the optimizer domains to the full-field scale and adding
  entropy/full-tilt starts did not repair them. The .95 witness is interior;
  the .98 witness uses the full rank tilt. This remains a gap in the
  tested upper bound, not a counterexample or a proof that 49% is impossible.

## Native multiplication results

`WideDiagonal16.h` uses the quadratic tower u^2=u+0x20 over AES GF256.
Setup stores four coefficient bytes. Four 128-bit records take two GFNI
products, coefficient shuffles, a byte-pair swap, and XOR.

`WideDiagonal32.h` preserves the exact field and map of `MediumField32.h`.
Setup represents multiplication by a scalar as a 4-by-4 GF256 matrix, stored
as four cyclic diagonals. Four records then take four GFNI products,
three within-word rotations, coefficient shuffles, and XOR. This exchanges
arithmetic for coefficient traffic: sixteen stored bytes per scalar,
versus four in the old tower implementation. It is not a restricted or
heuristic randomizer.

The table reports **microseconds per multiplication pass**. Ranges are
medians from two launches, not confidence intervals.

| 128-bit records | GF256 | New GF(2^16) | Old GF(2^32) tower | New GF(2^32) matrix |
|---|---:|---:|---:|---:|
| 4,096 | 1.019–1.021 | 1.131–1.132 | 5.497–5.500 | 1.905–1.908 |
| 65,536 | 21.246–21.403 | 21.574–21.584 | 85.286–85.780 | 31.208–31.218 |
| 262,144 | 82.043–85.299 | 85.940–87.053 | 342.108–343.632 | 125.154–125.293 |
| 1,048,576 | 631.388–635.617 | 673.598–675.561 | 1,439.107–1,444.647 | 1,070.558–1,072.502 |

Host: Peach, Ryzen 9 7950X, GCC 15.2, CPU 2, AVX512/GFNI. Each launch uses
three warmups and 101 samples, with 31/3/1/1 inner passes at the four sizes;
times are normalized per pass. Seeds are 113 and 7919. Benchmarks ran
strictly serially under a host lock. Each control uses separate ordinary
vector allocations, without an explicit huge-page request. Controls run
in fixed order, not interleaved on a shared address. Data and matching
GF(2^32) coefficients agree; outputs are consumed and their checksums match.
Setup, allocation, and coefficient preparation are excluded. Payload loads
and stores are included. No full-code speedup follows from this table.

The final native binary SHA256 is
`e1c78702396b38fab07d94de656b94c517d30fe5b21136e18c6dd3ce37481935`.

## Validation and replay

Release checks, registered CTest, and ASAN/UBSAN with leak checks pass.
Checks include all 65,536 GF(2^16) coefficients on 16 binary basis inputs,
independent polynomial/tower oracles, GF(2^32) matrix/tower equality,
random SIMD tests, exact aliasing, alignments, tails, guards, and unchanged
inputs/setup. Fourteen high-distance Python regression tests pass. The
independent proof audit caught and repaired overflowing floating dimension
moments at large fields and an incorrectly admitted middle endpoint.
The new dense helper also checks exact rational row agreement.

From the repository root:

```text
python -B -m unittest discover -s research/workstreams/char2_field_spin/theory -p "high_distance*_tests.py" -v
python -B research/workstreams/char2_field_spin/theory/high_distance_fields.py --bits 16 --outer 64 8 --outputs 1 --distance 21/50 --active 1 16 64 128 256 --density .003 .01 .03 .95 .99 --iterations 45
python -B research/workstreams/char2_field_spin/theory/high_distance_fields.py --bits 32 --outer 64 8 --outputs 1 --distance 9/20 --active 1 16 64 128 256 --density .003 .01 .03 .95 .99 .995 --iterations 45
python -B research/workstreams/char2_field_spin/theory/high_distance_dense.py --bits 32 --outer 64 8 --outputs 1 --distance 9/20 --density .99 .995 .999 1
cmake -S spin/experiments/char2_field_spin -B build/char2-high-distance -DCMAKE_BUILD_TYPE=Release
cmake --build build/char2-high-distance --target spin_char2_wide_diagonal16 -j 4
build/char2-high-distance/spin_char2_wide_diagonal16 --check
```

On a compatible Linux host, run the two benchmark launches **sequentially**:

```text
flock /tmp/spin-char2-benchmark.lock taskset -c 2 build/char2-high-distance/spin_char2_wide_diagonal16 --bench all 101 all auto 113
flock /tmp/spin-char2-benchmark.lock taskset -c 2 build/char2-high-distance/spin_char2_wide_diagonal16 --bench all 101 all auto 7919
```

Witness optimization can vary by starting points and numerical libraries.
These commands reproduce screens, not certified receipts. No production
profile, prior certificate, or frozen high-rate source was changed.
