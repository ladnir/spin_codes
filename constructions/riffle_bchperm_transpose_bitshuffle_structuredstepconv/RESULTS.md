# Parameter and kernel results

Fix a message length of \(2^{20}\) bits, an output length of \(2^{21}\) bits,
and a distance threshold of 188743. The numerical calculation uses the
modeled even, complement-symmetric \([256,128,38]\)-shaped spectrum.

## One-active parameter sweep

The 40-bit boundary follows a simple batching schedule.

| Step bits \(t\) | State bits \(s\) | Map width \(d\) | One-active margin |
|---:|---:|---:|---:|
| 4 | 17 | 21 | 40.5495 |
| 8 | 16 | 24 | 40.5539 |
| 16 | 15 | 31 | 40.5634 |
| 32 | 14 | 46 | 40.5911 |

Adding one state bit gives about 52 bits of margin:

| \(t\) | \(s\) | \(d\) | One-active margin |
|---:|---:|---:|---:|
| 4 | 18 | 22 | 51.9477 |
| 8 | 17 | 25 | 51.9530 |
| 16 | 16 | 32 | 51.9637 |
| 32 | 15 | 47 | 51.9851 |

The balanced point \((t,s)=(32,32)\) gives 80.8314 bits. The numerical
Chernoff grid and the real-valued spectrum make every value in this section
a diagnostic rather than a certificate.

At the balanced point, the one-active calculation remains above 40 bits
through relative distance 14.05%:

| Relative distance | One-active margin |
|---:|---:|
| 13.00% | 46.8134 |
| 14.00% | 40.4453 |
| 14.05% | 40.1343 |
| 14.10% | 39.8227 |

This frontier does not imply 14.05% distance for the complete code. The
multi-active spectrum contribution has not been evaluated.

## Bitsliced kernel baseline

The benchmark evaluates an in-place reverse inner chain on \(2^{21}\)
128-bit elements. It excludes setup and input restoration. Each row reports
the median of 15 trials. The process was pinned to one logical processor on
an Intel Core i7-13700H.

The code was compiled with MSVC 19.50 using `/O2 /Ob3 /arch:AVX2`. The hot
kernel uses an AVX2 128-by-128 transpose and 256-bit `VPCLMULQDQ`. Benchmarks
were executed serially.

| \(t\) | \(s\) | One-active margin | FieldMul ms | Toeplitz ms |
|---:|---:|---:|---:|---:|
| 16 | 15 | 40.5634 | 56.7623 | 49.0360 |
| 16 | 16 | 51.9637 | 56.4151 | 48.8910 |
| 32 | 14 | 40.5911 | 29.2352 | 25.8625 |
| 32 | 15 | 51.9851 | 29.4721 | 25.9712 |
| 32 | 16 | 60.5496 | 29.6148 | 26.1506 |
| 32 | 32 | 80.8314 | **27.1398** | **25.8285** |

The final vectorized implementation improved the \((32,15)\) FieldMul
kernel from 90.64 ms to 29.47 ms. It improved Toeplitz from 51.26 ms to
25.97 ms. The changes were:

1. Replace scalar Eklundh transposes with the AVX2 transpose schedule.
2. Use carryless multiplication for field reduction.
3. Process two independent bit lanes per vector carryless multiplication.
4. Reuse the padded transpose workspace across reverse steps.

At \(t=32\), every tested map fits in the same 64-bit bitslice. Increasing
the state from 14 to 32 bits therefore has no measured penalty. FieldMul
stores one 64-bit coefficient per step throughout this range. Toeplitz stores
two 64-bit words per step.

These results are a bitsliced baseline. They do not establish that a
physical transpose is the best evaluator for every map width.

## Direct block kernels

The direct evaluator keeps the 128-bit elements in their native layout. It
builds four-bit XOR-combination tables for the current step and applies the
binary matrix without a physical data transpose. Coefficients remain compact
and matrix rows are derived inside the timed region.

| \(t\) | \(s\) | Map width | Direct FieldMul ms | Direct Toeplitz ms |
|---:|---:|---:|---:|---:|
| 16 | 16 | 32 | 37.3897 | **15.6509** |
| 32 | 14 | 46 | 30.4932 | **17.9484** |
| 32 | 15 | 47 | 30.6888 | **18.4905** |
| 32 | 32 | 64 | 42.7634 | 34.3959 |

The crossover depends on the map width. At width 32, direct Toeplitz is 3.12
times faster than bitsliced Toeplitz. At widths 46 and 47, it is about 30%
faster. At width 64, the bitsliced Toeplitz kernel is faster. Direct FieldMul
wins at width 32, but the bitsliced FieldMul kernel wins at widths 46, 47,
and 64.

Thus, a physical transpose is not part of the preferred design by default.
It is useful only when the 128-way carryless-multiplication batching repays
its fixed transpose cost. The current best measured evaluator near the
40-bit one-active boundary is direct Toeplitz at \((t,s)=(32,14)\), taking
17.9484 ms. The balanced \((32,32)\) point still has more distance room, but
its fastest measured evaluator is bitsliced Toeplitz at 25.8285 ms.

FieldMul uses about half as many coefficient bits as Toeplitz and supports
the nonzero-coefficient invertible variant. Its direct row construction is
also more expensive; the current benchmark derives those rows within every
step rather than storing expanded matrices.

### Block-resident Toeplitz algorithm audit

Three Toeplitz evaluators were compared while keeping every value in
128-bit block form. Four Russians builds four-input XOR tables. The diagonal
kernel enumerates nonzero Toeplitz diagonals and XORs each shifted block
range with AVX2. The Karatsuba kernel recursively multiplies a block-valued
polynomial by the binary Toeplitz polynomial. Its base-case cutoff was tuned
at 8, 16, and 32 coefficients; 32 was fastest at width 46.

| \(t\) | \(s\) | Width | Four Russians ms | Diagonal AVX2 ms | Karatsuba ms |
|---:|---:|---:|---:|---:|---:|
| 16 | 16 | 32 | **14.8980** | 39.0523 | 36.7084 |
| 32 | 14 | 46 | **18.9647** | 33.3514 | 55.0741 |

The width-46 Four-Russians repetition is consistent with the earlier
17.9484 ms measurement. The comparison shows that the earlier time is not
caused by an accidental scalar implementation. At these widths, Karatsuba's
padding and temporary-polynomial work outweigh its reduction in elementary
products. The diagonal method performs more block XORs than Four Russians,
although its accesses are contiguous. Four Russians remains the preferred
block-resident evaluator.

### Four-Russians cost profile

The width-46 chain was split into state movement, Toeplitz-row generation,
four-input table construction, and indexed table application. The recurrent
movement floor touches all 32 input blocks per step and carries the 14-block
state, but replaces the dense map with one block XOR per output. It is an
implementation floor, not a code construction.

| Variant | Median ms |
|---|---:|
| Recurrent movement floor | 5.7071 |
| Compact Four Russians before local tuning | 19.2238 |
| Compact Four Russians after local tuning | 18.5588 |
| Four Russians with rows expanded during setup | **17.4525** |

Pre-expanding all 46 row masks uses about 23 MiB of additional setup data at
65,536 steps. It saves 1.11 ms despite the extra sequential reads. The
compact variant was improved by streaming the new Toeplitz coefficient bits
and initializing every output from its first table selection. Attempting to
specialize the partial final input group changed compiler code generation
and regressed the median to 22.6792 ms; that change was reverted.

A sampled serialized timestamp profile of the compact kernel reports:

| Phase | TSC ticks per step | Share of accounted ticks |
|---|---:|---:|
| State assembly and commit | 53.2 | 6.2% |
| Toeplitz-row generation | 104.9 | 12.3% |
| Four-input table construction | 147.2 | 17.3% |
| Indexed table application | 546.6 | 64.2% |

The profile samples one of every 256 steps and subtracts the minimum measured
timestamp overhead. Its ticks are for attribution, not an independent wall-
clock benchmark. The conclusion is stable: framework overhead and state
movement do not explain the 17--19 ms time. Most work is the dense block-XOR
map itself.

The benchmark measures only the inner chain. It excludes the outer encoder,
both permutation layers, and their transpose implementations.

## Validation

The C++ self-test compares Four Russians, diagonal AVX2, Karatsuba, and the
bitsliced implementation against dense reference maps at widths 31, 32, 46,
47, 48, and 64. It checks complete reverse chains at every benchmarked
\((t,s)\) point. It also compares optimized field reduction against
schoolbook polynomial arithmetic.

The algebraic checker applies the Rabin irreducibility test to every field
modulus. For each width, it also verifies full one-vector rank for all basis
vectors and sixteen deterministic random vectors.

## Next proof step

Extend the spectrum calculation beyond one active outer block. Condition on
the outer weights, retain the without-replacement row supports, and sum the
resulting region occupancies. This step determines whether the large margin
at \((32,32)\) survives the complete spectrum.
