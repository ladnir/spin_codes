# Exact-Map Fusion of Four IMT Updates

Follow-on: [packed GL32/BCH tuning](packed_bch_tune_REPORT.md) retains this
inner and reaches about 5.85--5.90 ms with wide output stores and an 8 MiB
coefficient table. The measurements below remain the original matched
comparison of sequential and fused inner updates.

2026-09-30. The fused scalar-table implementation takes **6.212096 ms**,
versus **6.753111 ms** for sequential R4 in the same executable: a saving
of 0.541016 ms (**8.01%**). Both compute exactly the same sampled encoder.
These are full precomputed transposed encodes at K=2^20, rate one-half,
with 128-bit XOR elements. Setup and allocation are excluded.

This optimization does not weaken or otherwise change the construction in
the [R4 proof record](packed_mixing/R4_CLOSURE.md). The retained sequential
R4 implementation and [its earlier comparison](packed_updates_REPORT.md)
are unchanged. Production defaults and paper claims are unchanged.

## The Change

The original inner applies four dependent transvections to its state.
Their product has an identity-plus-rank-at-most-four representation. Setup
precomputes that representation from the original masks, with no new
randomness. The transpose computes four dot products against the original
state, then combines their contributions in one state update. The feedback
syndrome is added in that same output pass.

Two bounded implementations are retained in [FusedR4.h](FusedR4.h):

- Scalar: build the 16 subset XORs of the four dot products; update the
  19 state words with fixed, unrolled table selections. Metadata occupies
  36 bytes per epoch (576 KiB at this size).
- AVX512: form two four-entry tables in vector registers and use byte
  indices expanded to 64-bit lane indices. Five unrolled groups update the
  state; the last group masks off the nonexistent twentieth word. Metadata
  occupies 96 bytes per epoch (1.5 MiB at this size).

The AVX512 candidate is not measurably better in this small campaign. The
scalar candidate has the smaller setup representation, so it is the current
candidate to carry forward. The benchmark prepares both representations
for every mode, outside timing; the metadata sizes above describe each
representation separately, not the total allocation of the test harness.

The probe includes the retained driver without modifying it. Mode zero
calls the original R4 kernel; modes one and two instantiate the two fused
kernels. Dispatch occurs outside the timed loop. The fixed emitter, route,
packed GL32/BCH kernel object, allocation order through the original buffers,
and streaming-store behavior are retained. There is no hot-path allocation
or runtime callback abstraction.

## Correctness

For seeds 1 and 17 at K=2^14 and K=2^20, both candidates pass:

1. Every epoch's complete 19-dimensional basis map in both forward and
   transpose orientations, against dense application of the original
   sequential masks. An additional SIMD half tests dense input.
2. Inner equality against the original dense and sequential references,
   and the original forward/fused-transpose adjoint identity.
3. Full-encoder equality against the independent scalar GL32 mixer and
   production BCH reference, including the entire untouched suffix.
4. Dense, sparse, and alternate-dense input patterns in check-only runs.

All timed modes also pass validation and have identical final checksums
for the same seed and call count. Independent source review checked the
product/adjoint orientation, syndrome order, table indices, and masked tail.
The GCC build passes; no new sanitizer campaign was run.

## Matched Measurements

Peach Ryzen 7950X, GCC 15.2, pinned CPU 15. Flags and kernel objects match
the earlier driver. Each process performs three warmups followed by 101
timed calls. Setup, allocation, input filling, reference checks, and checksums
are outside timing. There is no timed input reset between calls. All
benchmarks and compilation run serially under the three existing locks.

| Seed | Process order | Sequential R4 (ms) | Fused scalar (ms) | Fused AVX512 (ms) |
|---|---|---:|---:|---:|
| 1 | 0, 1, 2 | 6.768660 | 6.208795 | 6.155115 |
| 1 | 2, 1, 0 | 6.716132 | 6.208945 | 6.267293 |
| 17 | 0, 1, 2 | 6.767017 | 6.215246 | 6.221078 |
| 17 | 2, 1, 0 | 6.739205 | 6.258418 | 6.218843 |
| | Median of process medians | **6.753111** | **6.212096** | **6.219961** |

No samples or processes were discarded. The roughly 0.008 ms difference
between fused candidates is smaller than the observed process variation.
The longer confirmation runs, not the initial 31-call screen or the
instrumented phase runs, supply the headline result.

A separate 31-call instrumentation batch reports these phase means,
including warmups. They are not a decomposition of the headline medians.

| Mode | Inner + route range (ms) | Packed GL32 + BCH range (ms) |
|---|---:|---:|
| Sequential R4 | 3.246569--3.261652 | 3.493011--3.570751 |
| Fused scalar | 2.718425--2.728175 | 3.448741--3.449961 |
| Fused AVX512 | 2.712501--2.748212 | 3.544023--3.569696 |

The inner/route saving is consistently about half a millisecond. The BCH
kernel is unchanged; the smaller variations in its measured phase are not
claimed as an algorithmic improvement. Compiled encoder body sizes are
12,534 bytes for sequential R4, 13,318 for scalar fusion, and 13,577 for
AVX512 fusion. No routing-emitter symbol is present in the symbol listing.

## Reproduction and Provenance

[The runner](packed_fused_r4_run.sh) supplies locked serial `build`, `check`,
`screen`, `confirm`, and `profile` phases. Explicit invocation:

```sh
taskset -c 15 <build>/packed_fused_r4 20 1 101 1 0
```

The arguments are exponent, seed, call count, mode, and optional profiling
flag. Mode zero is the original encoder; zero calls runs correctness checks
for all three implementations. The source is
[packed_fused_r4.cpp](packed_fused_r4.cpp).

Remote build: `/tmp/spin-fused-r4-c7oraT`.
Logs under its `measurements/` directory:

- Correctness: `fused-check-AzGpLv`.
- Initial screen: `fused-screen-jdTOto`.
- Matched confirmation: `fused-confirm-ZM2BjT`.
- Phase instrumentation: `fused-profile-lS7SIu`.

All logs, the symbol listing, and compiler stack-usage output were copied
to ignored `tmp/packed-hill/fused-r4-performance/`. Raw data is not committed.

| File | SHA256 |
|---|---|
| `FusedR4.h` | `a61ae7a633057be5c8c3062cce7555bd9f475a0339aaa193fa16c22557e63781` |
| `packed_fused_r4.cpp` | `7d98ce42439db804c08c57069e2e6f2cd08ea3d4921123d7ed3978968677ecc1` |
| Unchanged `packed_driver.cpp` | `fda1849261c6777d42d13c753452790a3401c33b783174bcc4e07fb29efb3f8b` |
| Executable | `7b8c1b3addf2d83f0126a0bd4dad65940485407c14bbc101330b840cb01ffef2` |
| Unchanged packed-kernel object | `ebe4ef88a815da7ceb8abfb486fc199dbd6d7a6b27d35d0d227a98a694b3a3b5` |

Next: retain this exact-map optimization with the completed R4 certificate.
The packed GL32/BCH phase, roughly 3.5 ms, is now the larger performance
target. Further tuning should use this matched baseline, without changing
the construction unless its new proof benefit justifies the change.
