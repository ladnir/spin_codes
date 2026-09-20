# Bare BCH-256/RM2Sub performance

Measured on Peach on 2026-09-07. All 12 cells pass the independent dense oracle;
the GCC and MSVC output hashes agree. No fanout is included.

## Online encoding time

Times are milliseconds for the transposed block encoder. Each 128-bit block
carries 128 parallel binary instances. Values are medians of three run medians.

| Inner (t,s) | K=2^16 | K=2^18 | K=2^20 | K=2^20 interpretation |
|---|---:|---:|---:|---|
| t64_s20 | 0.625 | 2.577 | 12.059 | full reference, 50.487298 bits |
| t128_s19 | 0.560 | 2.322 | 11.259 | full certificate, 50.448203 bits |
| t64_s16 | 0.584 | 2.400 | 11.547 | shortlisted, not fully certified |
| t256_s14 | 0.517 | 2.143 | 10.545 | exploratory watch |

S19 takes 6.6% less online time than the S20 reference at K=2^20.
The wider-step candidate is faster, but its timing does not resolve its weaker proof evidence.
Use the selected (128,19) map as the main certified performance point; retain (64,20) as a reference.

## Variation and setup costs at K=2^20

| Inner | Range of run medians (ms) | Setup (ms) | Retained setup (MiB) | Workspace (MiB) |
|---|---:|---:|---:|---:|
| t64_s20 | 12.017–12.063 | 46.23 | 14.50 | 40.00 |
| t128_s19 | 11.105–11.276 | 43.91 | 13.19 | 40.00 |
| t64_s16 | 11.517–11.562 | 39.16 | 14.00 | 40.00 |
| t256_s14 | 10.436–10.615 | 42.94 | 12.44 | 40.00 |

Setup is measured separately and includes creating oracle and alternative routing data.
The benchmark discards that data before online timing. The setup column excludes this final
deallocation and workspace construction. Input/output buffers add 48 MiB at K=2^20.
Those buffers are separate from the listed workspace. Setup storage before compaction is
approximately 36–39 MiB; PERFORMANCE.json records exact byte counts and all final observations.

## Measurement method

- Ryzen 9 7950X, Linux, one thread pinned to CPU 15; GCC 15.2.0.
- CMake Release: `-O3 -DNDEBUG -march=znver4 -mavx2 -mpclmul -mvpclmulqdq`.
- Three independent serial processes, three warmups, 101 timed calls per cell per process.
- Identical inputs, routing seeds, and coefficient seeds across repetitions.
- Setup, workspace allocation, correctness checks, and output hashing are outside online timing.
- Every cell checks for other benchmark executables; a process lock prevents concurrent copies.
- The same buffers and workspace are reused. This is steady-state latency, not cold-cache latency.
- The tuning sweep precedes the final runs; final numbers are not minima selected from that sweep.

## Optimizations retained

All configurations share routing, setup, the paired BCH transpose, and the RM2Sub template.
Each selected map supplies compile-time constants and an exactly checked quadratic XOR circuit.
The kernels use pruned SIMD zeta stages, unrolled emission, fixed nibble tables, and
a templated emitter. The hot call performs no allocation or type-erased callback dispatch.
Configuration and layout dispatch occur once per encode.

The bounded tuning effort selected these choices:

- Packed 24-bit routing rather than 32-bit indices at K=2^20.
- 256-row tiles at K=2^16 and K=2^18; 2048-row tiles at K=2^20.
- Optimized nibble grouping for S19 only. Three A/B repetitions showed about 3% lower latency.
  Other maps keep sequential grouping because fewer nominal lookups did not give a reliable win.
- Setup compaction removes the dense route, scalar coefficient list, and unused routing format.

The configuration constructor exposes tile size, and CMake exposes `SPIN_GROUPED_A`.
These are measured defaults on this host, not a global optimum for every processor.

## Construction and correctness

The BCH matrix is not the imported Eq3 matrix. Let p and q denote the length-255 BCH
generator polynomials for designed distances 37 and 39. Let P and Q denote their
parity-extended codes. The implemented code spans Q and the parity extensions of
five representatives `p*x^j`, for `j=0..4`, and uses a systematic basis.
Generation checks rank 128, Q containment, P containment, and inclusion of the all-one word.
This is the intermediate-code family used by the existing bounds.

The generated paired transpose uses 3029 XORs in its scalar circuit and processes two rows
with AVX2. Symbolic checks verify its matrix and transpose. Separate symbolic checks verify
the pruned zeta schedule and the selected inner maps. Runtime tests compare complete outputs
against dense BCH multiplication and a separately implemented scalar RM2Sub recurrence.
Tests also cover both route layouts, tile changes, compaction, multiple setups, boundary
impulses, zero input, and invalid or overlapping buffers.

The known 50.487298-bit label applies to the selected S20 reference at K=2^20.
Subsequent proof work certified the selected S19 map at K=2^16, 2^18, 2^20, 2^22, and 2^24.
See [the paper handoff](../bch_rm2sub_bridge/PAPER_HANDOFF.md) for the exact scope and evidence.
Timings and implementation support remain limited to the sizes described above.
Historical Eq3/fanout timings concern a different map and are not used as the comparison baseline.

## Reproduction and next step

See README.md for build commands and the API. PERFORMANCE.json includes source hashes,
binary hashes, exact memory sizes, quantiles, and the 36 final cell observations.
Local tuning and correctness logs are in `out/bare_spin_receipts/`.

Next, integrate S20 and S19 into the calling application and measure end-to-end throughput.
The S19 distance/setup proof is complete at the listed sizes; retain S16 and T256 as comparison points.
