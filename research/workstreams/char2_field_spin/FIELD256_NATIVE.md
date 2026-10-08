# Native implementation of the higher-rate checkpoint

The GRS[256,176] checkpoint has a complete K20 certificate for 10% symbol
distance at rate 11/16. Its first native implementation is correct but
substantially slower than the binary production control. A shorter output
does not by itself compensate for the field outer and inner arithmetic.
The [certificate](FIELD256_K20_CERTIFICATE.md) remains a rate improvement;
it is not a performance result or a new half-rate/10% guarantee.

## Implemented map and validation

The isolated experiment is `spin/experiments/char2_field_spin/cauchy256.cpp`.
It uses the certified systematic GRS map, width-sixteen inner, eight state
symbols, coordinate diagonals, and one regional routing stage. No binary
production implementation or earlier reference construction is replaced.

Setup pads the requested source to a multiple of 176*64 symbols. Encoding
uses caller-owned scratch and performs no allocations. Transposed encoding
returns the requested source prefix; its in-place form preserves the input
suffix beyond that prefix. State carries across all region boundaries.

The independent scalar oracle computes the outer by direct Lagrange
interpolation, rather than by the generated native factorization. It checks
the field adjoint and complete encoders with 1/4/8/16-byte records. Native
128-bit tests check all three route schedules, four input alignments,
in-place and out-of-place operation, guards, boundary sizes, repeated calls,
and requested-prefix padding. Separate K16/K18/K20 checks compare full
native outputs and in-place suffixes. All five exact-map variants pass
these checks in Release and ASAN/UBSAN. The registered CTest also passes.
Smaller-size performance tests do not supply their missing distance certificates.

The native records use a GF256-coordinate basis for the compatible extension
field. Timings do not include conversion from an application's different
polynomial basis. Such conversion, if required by an integration, is an
additional cost rather than an implicit benefit of these measurements.

## Initial timings

Measurements use Peach's Ryzen 9 7950X, GCC 15.2, AVX512/GFNI, CPU 2,
seed 113, five warmups, and 31 timed calls. All benchmarks run serially.
Times exclude setup, allocation, initialization, and payload restoration.
The initial executable SHA256 is
`6ff855e5b8fdc4a42b0db37a9ca95f55cfd4a6c57593c6307ce48bf6d5635dcd`.

The phase diagnostic uses scatter at K16/K18 and regional non-temporal
routing and output at K20. Medians of separate phases need not sum to the
median of complete calls.

| Requested K | Padded K | N | Full (ms) | Inner and route (ms) | Outer (ms) |
|---:|---:|---:|---:|---:|---:|
| 2^16 | 67,584 | 98,304 | 0.610 | 0.269 | 0.338 |
| 2^18 | 270,336 | 393,216 | 2.496 | 1.081 | 1.415 |
| 2^20 | 1,058,816 | 1,540,096 | 11.847 | 4.911 | 6.932 |

A separate matched comparison uses the public `spin::Code` PacketRs8T64
control at the same requested K. Both encoders receive the same payload
address and initial bytes. Five-call batches alternate order; restoration
occurs outside timing. Both use normal pages at K16 and prefer huge pages
at K18/K20. The field candidate uses regional non-temporal routing and
output throughout this comparison, which is not its best small-size route.

| Requested K | Binary control (ms) | Initial field candidate (ms) |
|---:|---:|---:|
| 2^16 | 0.0990 | 0.7225 |
| 2^18 | 0.4655 | 3.2623 |
| 2^20 | 3.1618 | 11.8847 |

The initial outer minimizes field products through XOR-circulant
factorization, but executes 8,322 Boolean operations and substantial scratch
traffic per four rows. Product count alone therefore selected a poor
native implementation.

## Exact-map optimization results

Three changes preserve the certified map and all setup distributions:

1. The dense packed inner reuses record broadcasts, folds eta into fixed
   coefficients, and uses ternary XOR accumulation. Its 66 GFNI512 products
   per packet exceed the original count, but replace substantial packing
   and XOR work. Scatter inner time falls from 270 to 144 microseconds at
   K16 and from 1081 to 586 microseconds at K18.
2. The additive-FFT outer uses 1203 GFNI512 products and 2851 XOR operations
   per four rows, versus 503 and 8322. Independent replay checks all 176
   forward and 256 transpose basis vectors, both translation maps, and
   random adjoints: 672 exact checks pass.
3. The parity-only FFT computes the 80-coordinate parity contribution and
   adds systematic coordinates during final output stores. It stages 80
   rather than 256 vectors and reuses its output as FFT workspace. Additional
   wrapper scratch falls to 4 KiB. Its 128 exact scalar checks pass. Native
   tests check the parity map independently and include partial output tiles.

The final selector is `dense-fft-fused`. Original `basis`, inner-only
`dense`, outer-only `fft`, and combined `dense-fft` remain available as
exact-map controls. The default remains `basis`; timed commands select
the optimized variant explicitly. There is one variant dispatch per full
encode, with compile-time specialization inside the scans.

The following final matched measurements use 101 timed calls and otherwise
the same five-call alternating protocol described above. The route is
scatter at K16/K18 and regional non-temporal routing/output at K20.

| Requested K | Binary, seed 113 (ms) | Field, seed 113 (ms) | Binary, seed 7919 (ms) | Field, seed 7919 (ms) |
|---:|---:|---:|---:|---:|
| 2^16 | 0.09835 | 0.27962 | 0.09868 | 0.27926 |
| 2^18 | 0.46078 | 1.23439 | 0.46170 | 1.22676 |
| 2^20 | 3.17438 | 7.32216 | 3.15411 | 7.38317 |

The final executable SHA256 is
`8130c53329d3a8ec99282e731bf3c9bff757d624b382d6fd4d67ab7f8726d818`.
The matched K18 batches have visible tail noise: field p90 is 1.50--1.52 ms
and binary p90 is 0.939--0.941 ms. The gap is much larger than that variation.

A separate final 51-call phase diagnostic gives inner/outer medians of
0.143/0.136 ms at K16, 0.585/0.607 ms at K18, and 3.018/3.911 ms at K20.
These use the diagnostic allocator and are not substituted for the matched
full-call results. Both arithmetic stages remain material costs.

Thus the implementation changes are substantial wins over the first field
kernel, but the resulting code remains about 2.3--2.8 times slower than
binary SPIN in the matched runs. The shorter output and higher certified
rate are real; competitive encoding performance has not been established.
No production profile is promoted or changed.

## Further exact-map inner experiment

The sixth variant, `fft-inner-fused`, replaces the dense inner with an
additive-FFT factorization in the original proved state coordinates.
It does not move the fresh diagonal into a different state basis.
Including input diagonals, its packet needs 34 rather than 70 GFNI512
products, but also 46 lane shuffles and substantial sequential dependence.
Independent arithmetic checks, native checks, and ASAN/UBSAN pass for
all six variants. The measured result rejects this optimization:

| K | Dense full (ms) | FFT-inner full (ms) | Dense inner + route (ms) | FFT inner + route (ms) |
|---:|---:|---:|---:|---:|
| 2^16 | 0.27885 | 0.30300 | 0.14363 | 0.16808 |
| 2^18 | 1.20633 | 1.29788 | 0.58645 | 0.68175 |
| 2^20 | 6.97117 | 7.33452 | 3.02330 | 3.39510 |

These are separate 51-call phase diagnostics, seed113, with scatter at
K16/K18 and region-nt-all at K20, not the preceding same-buffer figures.
Seed7919 with reversed variant order confirms the regression. Checksums
match between exact-map variants. `dense-fft-fused` stays the winner.
The executable SHA256 for this comparison is
`7586a3729d50b1bdf9ebe43692944ce896fade9db08aa5ee66612193e98905ec`.

Moving the random diagonal into a cheaper coefficient basis was also
rejected: [explicit rank witnesses](CHEAPER_INNER_MAP.md) invalidate the
proof properties used here. This is not evidence that an FFT operation
count alone gives a faster encoder.

The subsequent [GRS128/87 checkpoint](FIELD128_NATIVE.md) changes the
outer and source layout, closes its own certificate, and improves K20
matched time to 6.01--6.10 ms. The current page retains the earlier map
and its measurements as a reproducible control.

## Reproduction and next decision

From the repository root, configure the isolated experiment with the
standalone library as the optional control:

```sh
cmake -S spin/experiments/char2_field_spin -B out/char2-field \
  -DCMAKE_CXX_COMPILER=g++-15 -DCMAKE_BUILD_TYPE=Release \
  -DSPIN_CHAR2_PUBLIC_ROOT="$PWD/spin" -DSPIN_TUNE=znver4 \
  -DSPIN_BUILD_TESTS=OFF -DSPIN_BUILD_BENCHMARKS=OFF
cmake --build out/char2-field --target spin_char2_cauchy256 -j4
out/char2-field/spin_char2_cauchy256 --check
out/char2-field/spin_char2_cauchy256 --check-large
taskset -c 2 out/char2-field/spin_char2_cauchy256 --compare 65536 113 101 scatter full dense-fft-fused
taskset -c 2 out/char2-field/spin_char2_cauchy256 --compare 262144 113 101 scatter full dense-fft-fused
taskset -c 2 out/char2-field/spin_char2_cauchy256 --compare 1048576 113 101 region-nt-all full dense-fft-fused
```

For ASAN/UBSAN, use a separate build without the optional production control
and add `-DCMAKE_CXX_FLAGS=-fsanitize=address,undefined` and
`-DCMAKE_EXE_LINKER_FLAGS=-fsanitize=address,undefined`. Run both check modes,
not benchmarks. The standalone phase syntax replaces `--compare K` with
`K` and `full` with `phases`.

The next design pass should seek a cheaper outer or packet map while
preserving the higher-rate objective and controlling the known regional
window obstruction. Current measurements do not justify inflating the
119-bit margin or promoting this profile. Preserve this complete proof
and its exact-map native implementation as the reference for that pass.
