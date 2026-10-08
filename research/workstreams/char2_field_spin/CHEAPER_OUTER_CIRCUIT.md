# Direct additive-FFT outer for GRS[128,87]

The isolated GRS[128,87] outer has an exact four-row packed circuit with 516 GFNI512 multiplications. Its parity-only transpose uses 1,188 XOR instructions and 2.5 KiB of nested scratch. The generated header is complete and passes independent interpolation checks. No native build or benchmark was run for this note.

This outer has rate `87/128`, above two thirds. It is a separate higher-rate profile, not a change to the archived GRS[256,176] map. Whole-code probability bounds are maintained separately. The arithmetic uses only AES GF256 constants, applied independently to each record byte; it introduces neither general extension-field multiplication nor an additional routing stage.

## Map and interface

The source evaluations are at byte labels `0,...,86`. Output evaluations are at `0,...,127`. Let `E` be this systematic interpolation map, and write `E=[I;P]`, where `P` is `41 x 87`.

[Cauchy128OuterFft.h](../../../spin/experiments/char2_field_spin/Cauchy128OuterFft.h) exports three ordinary noinline functions in `spin::research::char2`:

```cpp
c128OuterForwardFft(const __m512i* input87, __m512i* output128);
c128OuterTransposeFft(const __m512i* input128, __m512i* output87);
c128OuterParityTransposeFft(const __m512i* input41, __m512i* output87);
```

Each vector holds one 128-bit record from each of four independent outer rows. The parity hook computes `P^T`, not `E^T`: its input points are `87,...,127`, and the frontend must add the 87 systematic contributions. Public input and output ranges must not overlap. Unaligned pointers are accepted. There is no heap allocation or caller-provided scratch.

The implementation uses fixed-index, fully unrolled arithmetic. FFT leaves hold at most sixteen data vectors; larger butterfly stages process eight pairs at a time. Scratch belongs to ordinary functions, not coroutine frames.

## Prefix decomposition

For a power of two `n`, let `T_n` interpolate on `0,...,n-1` and evaluate on `n,...,2n-1`. Its matrix is a symmetric XOR-circulant, so `T_n^T=T_n`. The circuit is an inverse normalized additive FFT followed by a forward FFT on the other coset.

Split `87=64+23` and `23=16+7`. The last completion, from seven points to eight, needs no multiplication: the eighth evaluation of a degree-below-seven polynomial is the XOR of the first seven. This follows from the constant derivative of the eight-point subspace polynomial.

Let `P_23` be the parity block of systematic GRS[32,23]. For its transpose, split the nine inputs as `(v,w)`, with lengths one and eight. Compute

\[
q=T_8w,\qquad q_7\leftarrow q_7+v,
\qquad u_i=q_i+q_7\quad(0\le i<7).
\]

Then

\[
P_{23}^T(v,w)=\bigl(T_{16}(u\mathbin\Vert v\mathbin\Vert w),\ u\bigr).
\]

For the desired parity transpose, split its 41 inputs into `v` of length nine and `w` of length 32. Compute

\[
q=T_{32}w,\qquad q_{23:32}\leftarrow q_{23:32}+v,
\qquad u=P_{23}^Tq_{23:32}+q_{0:23}.
\]

The result is

\[
\boxed{P^T(v,w)=\bigl(T_{64}(u\mathbin\Vert v\mathbin\Vert w),\ u\bigr).}
\]

The header saves `u` in output coordinates `64,...,86`. It then uses output coordinates `0,...,63` as the workspace for `T_64`. Only `q[32]` and a nested `q[8]` require separate arrays.

There is a second exact description useful for validation. The complete 128-point evaluation set has a constant nonzero vanishing derivative. Complementary Lagrange interpolation therefore identifies `P^T` with interpolation from the 41 parity points to the 87 source points, without diagonal scalings. Translating every point by byte 127 maps these sets to reversed prefixes `0,...,40` and `41,...,127`. The checker compares all 3,567 matrix entries with that description as well.

## Arithmetic and memory counts

All counts below are per four independent outer rows. XOR counts are binary XORs; they exclude moves, loads, stores, coefficient materialization, and compiler spills.

| Operation | GRS[256,176] FFT | GRS[128,87] FFT |
|---|---:|---:|
| GFNI512, either direction | 1,203 | 516 |
| Forward XOR | 2,771 | 1,234 |
| Full transpose XOR | 2,851 | 1,275 |
| Parity-only transpose XOR | 2,675 | 1,188 |
| Parity hook's explicit nested scratch | 4 KiB | 2.5 KiB |
| Parity input + output + nested scratch | 20 KiB | 10.5 KiB |
| Forward explicit nested scratch | 13 KiB | 7.5 KiB |

The smaller outer requires about twice as many four-row calls at a fixed source length. Consequently, the per-call reduction is not the end-to-end reduction. With the current dense state-eight inner, a sixteen-record packet costs 70 GFNI512 instructions, including its input diagonal.

| Per padded source record | GRS[256,176] | GRS[128,87] |
|---|---:|---:|
| Outer GFNI512 | `1203/704 = 1.70881` | `43/29 = 1.48276` |
| Inner GFNI512 | `70/11 = 6.36364` | `560/87 = 6.43678` |
| Combined GFNI512 | `5683/704 = 8.07244` | `689/87 = 7.91954` |
| Combined Boolean operations | `5283/704 = 7.50426` | `2491/348 = 7.15805` |

These totals include the systematic XOR additions but exclude output-layout conversion. The GFNI volume falls about 1.9%; the Boolean volume falls about 4.6%. The inner expands slightly more data because the new rate is lower than `11/16`. No twofold speedup follows from these counts.

At requested K20, padding to a multiple of 64 rows gives:

| Profile | Rows | Padded source records | Encoded records |
|---|---:|---:|---:|
| GRS[256,176] | 6,016 | 1,058,816 | 1,540,096 |
| GRS[128,87] | 12,096 | 1,052,352 | 1,548,288 |

An implementation may omit outer calls whose entire output lies beyond the requested source length. The preceding normalized counts deliberately use complete padded rows. Input-diagonal multiplication remains part of measured encoding work.

## Four-row-interleaved output order

The new profile deliberately serializes source coordinates in four-row-interleaved order. Within a row group, the records are ordered by `(coordinate, row-lane)`. Each of the 87 output vectors is therefore already one contiguous aligned four-record block. The frontend adds its systematic plane and stores the vector directly.

This avoids the old four-by-four record transposes, which require eight lane shuffles per sixteen source records. It also avoids an alignment trap: a row-major stride of 87 records has byte stride 1,392. Its four row starts have offsets `0,48,32,16` modulo 64. An aligned non-temporal ZMM store is invalid at three of those starts.

Interleaving is a fixed permutation of source coordinates, not another random routing stage or a postprocessing pass. The full generated code is unchanged by source-coordinate permutation. Restricting any prefix of the new source order is a subcode, so an established full-code distance bound applies to that shortening as well.

The frontend must encode this order explicitly in its scalar oracle. It must mask or copy the final partial four-record vector when the requested length is not divisible by four. It must not write beyond the requested output, even when the outer operates on a complete four-row group.

## Why not extend to 88 first?

One can reconstruct source evaluation 87 from the first 87 values and use an 88-point circuit. The reconstruction coefficients are constant on the first ten eight-point blocks and the final seven-point block:

```text
101, 211, 250, 76, 255, 72, 109, 218, 12, 13, 1.
```

Their GF2 span has basis `(1,12,72,76,101,211)`. Thus the reconstruction transpose needs five nontrivial fixed products, thirteen binary XORs to form its distinct coefficients, and the final output accumulations.

A direct 88-point prefix FFT uses 491 products. Adding this constraint gives 496 products and 1,216 parity-transpose XORs. It saves twenty products against the chosen direct-87 circuit, but adds 28 XORs and an 88-value intermediate. The older coefficient-shared 128/88 circuit gives a separate 237-product upper count after this constraint; its higher Boolean fan-out and bounded-register schedule remain unimplemented.

The generated candidate uses the direct-87 decomposition because it has a simple bounded working set and independent complete validation. The alternative counts are not evidence of a faster native schedule.

## Validation and next gate

Run:

```sh
python -B spin/experiments/char2_field_spin/generate_cauchy128_outer_fft.py --check
```

The replay passes 3,919 exact small-map, coordinate-basis, complementary-interpolation, random dense-map, and adjoint checks. It also checks the operation counts and requires the saved header to match regeneration exactly. An independent read-only audit checked the wrapper ranges, scratch disjointness, identity prefix, and output-head reuse.

Native compilation, full-profile scalar-oracle tests, requested-length boundary tests, and serial timing remain separate gates. The next comparison should retain the existing dense inner and compare the complete isolated 128/87 profile against the frozen 256/176 profile. The arithmetic gain is modest; reduced cache working set and eliminated output shuffles are the stronger implementation hypotheses.
