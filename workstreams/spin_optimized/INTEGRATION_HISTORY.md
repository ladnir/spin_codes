# SPIN integration history

This document preserves the earlier integration notes and measured snapshots.
Some option descriptions refer to the build at that stage, before aligned-length
dispatch was introduced. Use [README.md](README.md) for current defaults and
[ALIGNED_LENGTHS.md](ALIGNED_LENGTHS.md) for the current routing policy.

These targets integrate the four-row BCH kernel into the selected half-rate
SPIN encoder: BCH [256,128], randomized bit-transpose routing, and the
weight-five IMT inner with (t,s)=(128,19). Existing configurations preserve their
binary linear maps. Explicit `T64S12` and `T64S12R2` options select the K16
inners described below. These options do not replace the submitted artifact or its source pins.
Each 128-bit block carries 128 parallel binary instances; the transposed API
maps 2K input blocks to K output blocks.

## Build and use

The transpose target now accepts naturally aligned lengths through
`MessageLength{K}`. It checks structural validity, not certificate coverage.
See [ALIGNED_LENGTHS.md](ALIGNED_LENGTHS.md) for alignment, automatic index width,
and the preserved power-of-two paths. The existing exponent constructor remains.
Forward work is parked; `SPIN_FORWARD_FOUR` is experimental and defaults to OFF.

For other XOR element types, `Spin::generic()` creates an explicitly requested,
owned fallback plan. See [GENERIC_ELEMENTS.md](GENERIC_ELEMENTS.md) for its typed
workspace API. The optimized 128-bit kernels and dispatch are unchanged.

For Hypercat's 256/512-bit **forward** kernels and the equal-work comparison
including row-to-column assembly, see [WIDE_FORWARD.md](WIDE_FORWARD.md).
The optional `spin_wide` target supports both S19 and the selected K16
`T64S12R2` inner, preserving their binary maps. It is forward-only.
See [WIDE_TUNING.md](WIDE_TUNING.md) for the K16 extension and K20 tile tuning.

From the repository root, on Linux x86-64 with GCC, CMake 3.20+, and Python 3:

```sh
cmake -S workstreams/spin_optimized -B out/spin-optimized \
  -DCMAKE_BUILD_TYPE=Release -DSPIN_BCH_AVX512=ON
cmake --build out/spin-optimized -j2
ctest --test-dir out/spin-optimized --output-on-failure -j1
```

Link `spin_half_transpose` and include its generated `Spin.h`. Existing
constructor calls and encoding calls remain source-compatible. Rebuild callers;
this is not an ABI-compatible replacement for previously compiled libraries.
As with the earlier alternatives, link only one library defining
`bare_spin::Spin` into each executable.

`SPIN_BCH_AVX512` defaults to `OFF`. When enabled, setup selects the four-row
kernel only if AVX-512F and AVX-512VL are available and the tile contains at least four rows.
Two-row tiles still work through the original kernel. The selected backend and
its offset representation remain fixed for the object's lifetime.

```cpp
bare_spin::Spin code(bare_spin::Configuration::T128S19, 20);
bare_spin::Spin::Workspace workspace(code);
// Allocate code.codeBlocks() input blocks and code.messageBlocks() output blocks.
// Use encode(), or encodeInplace() with a code.codeBlocks()-sized buffer.
```

An optional final constructor argument selects `BchBackend::Auto` (default),
`Avx2`, or `Avx512`. `bchBackend()` reports the actual selection.
An explicit `Avx512` request throws if the build, CPU, or tile cannot support it;
it never silently changes to another backend. Both Packed24 and Indices32
schedules support validation, compaction, and in-place encoding.

The baseline requires x86-64-v3, PCLMULQDQ, and VPCLMULQDQ; it is not a generic
x86 binary. Only the four-row BCH kernel and its separately dispatched hot
loop enable AVX-512F/VL. Setup and the fallback remain at the baseline ISA.
LTO is disabled across that boundary. `SPIN_TUNE=znver4` selects Zen 4 scheduling
without enabling additional instructions; the default is `generic`.
The current validation host uses GCC 15.2. MSVC and ARM are not supported here.
Clang is accepted by the build but has not been validated in this integration.

## Two-round K16 option: 49.33 bits

Use `Configuration::T64S12R2` at message exponent 16 for the selected
`(t,s,r)=(64,12,2)` inner. It retains the BCH [256,128] outer, routing,
and fixed A/B maps of `T64S12`. Each update applies two independent
transvections before the single feedback addition. Output precedes the
update; the initial state is zero and there is no flush.

The full setup-failure bound is **2^-49.3275773868** at 10% relative distance.
It covers all outer occupancies, not only one active row. See
[the full proof and prototype report](../inner_design/k16_parameters_20260919/FULL49_RESULT.md).
This certificate applies at K=65536. The transposed library accepts other
aligned lengths without asserting a certificate; the caller selects the map.
The parked bidirectional library retains its previous size restrictions.

```cpp
bare_spin::Spin code(bare_spin::Configuration::T64S12R2, 16);
bare_spin::Spin::Workspace workspace(code);
```

Both libraries support this option. The bidirectional library also exposes
`forward`, `forwardReference`, and `forwardBits`, including after compaction.
`wideView()` remains S19-only for compatibility with old consumers.
With `SPIN_BUILD_WIDE=ON`, `wideForwardView()` and `WideWorkspace` also support
`T64S12R2`, with 256-bit and 512-bit logical elements.

`k16_r2.py` adds the configuration after the one-round generator.
`ImtRounds.h` defines explicit map traits and inlined transvection updates.
The reverse kernel applies the second transvection's transpose first.
The forward kernels apply the first transvection first. Neither direction
adds a runtime round-count loop, another permutation, or hot-path allocation.

The benchmark selector is `6412r2`. `r2_check.sh ROOT` checks both libraries
with direct, tiled, AVX2-only, feature-disabled, and sanitizer builds before
running serial transpose benchmarks. The tests include frozen output hashes
from the one-round implementation and the earlier two-round prototype.
Use `test_k16_map.py --proof PATH` with `target49_sub12_r2_full.json` to
check the fixed maps against the full certificate.

The full build matrix passed. Matched integrated transpose timings remain
about 0.34 ms; forward paths were validated but not timed. See
[R2_RESULTS.md](R2_RESULTS.md) for the measurements, memory costs, and checks.

## What changes

`SPIN_K18_DIRECT=ON` selects an additional direct-routing specialization at
`K=2^18` when the AVX-512 BCH backend is active. Set it to `OFF` to keep
the tiled control. K16 keeps its separate specialization; K20 remains tiled.
The K18 change preserves the `(128,19,1)` inner and its full 50.189-bit margin
at 10% relative distance. It does not apply the K16 two-round configuration
at a new length.

The K18 inner writes into an 8 MiB buffer in the four-row BCH layout.
The outer then reads that buffer sequentially, avoiding an intermediate
8 MiB read and 8 MiB write. It uses plain 32-bit destinations without
software prefetching. `k18_screen.sh` records the candidate search;
`k18_confirm.sh` checks the selected build, fallbacks, and sanitizers before
running matched serial benchmarks. `k18_summary.py` checks output equality
and reports the measurements. Detailed results are in [K18_RESULTS.md](K18_RESULTS.md).

K20 keeps tiled routing. A subsequent screen rejected direct routing,
streaming stores, buffered writes, and alternative tile sizes at this length.
`SPIN_K20_EXPERIMENT=0` is the retained default; other values are research
overrides. See [K20_RESULTS.md](K20_RESULTS.md) for timings and replay scripts.
The subsequent [K20 profile](K20_PROFILE.md) finds three comparable phase
costs and identifies BCH spill traffic as the next bounded optimization
target. `SPIN_STAGE_PROFILE=ON` enables a separate diagnostic build; it is
off by default and is not a performance configuration.
The [four-row BCH locality screen](BCH_LOCAL_RESULTS.md) subsequently rejected
local output groups and forced input reloading. `SPIN_BCH_SCHEDULE=global`
remains selected. Other schedule values are research overrides. AVX-512 test
builds now include an exhaustive BCH coordinate/lane check.

The old circuit processes two independent BCH rows per call. The new circuit
processes four, using contiguous 512-bit loads. Setup changes the existing
scratch-offset schedule so each coordinate's four rows occupy adjacent lanes.
The scatter writes directly into this representation; there is no extra
transpose pass or retained buffer. The public routing permutation is unchanged.

The new XOR circuit uses less aggressive sharing to reduce register spills.
Generation checks each output's linear form against the original BCH matrix.
For existing configurations, the inner maps, transvections, setup sampling, and output ordering are unchanged.
Consequently this optimization does not require a different distance certificate.

`generate.py` derives the integration from the pinned implementation and fails
if a patch anchor changes. It copies the required headers into the build tree
to avoid mixing incompatible class definitions. Circuit synthesis reuses
`../inner_design/bch_20260919/generate.py` and
`../../scripts/probe_bch_forward_xor_circuit.py`. No numerical data files are
needed. Generated code belongs to the build directory, not the source tree.

## Validation and measurements

`validate.sh ROOT` builds AVX-512 enabled, disabled, and detection-suppressed
variants, then benchmarks serially. Detection suppression is a test hook, not
a substitute for testing on a physical AVX2-only host. `sanitize.sh ROOT` runs
the enabled build with AddressSanitizer and UndefinedBehaviorSanitizer.

The tests cover dense-reference equality at K=2^16, 2^18, and 2^20; explicit
backend equivalence; two-row fallback; tiles 4/8/16/512; both schedule layouts;
compaction; boundary inputs; and in-place suffix preservation.
`isa_check.py BUILD` inspects the compiled objects for misplaced EVEX instructions.

For benchmarks, add `-DSPIN_BUILD_BENCHMARK=ON`, then run, for example:

```sh
out/spin-optimized/spin_benchmark 20 auto 101 1
out/spin-optimized/spin_benchmark 20 avx2 101 1
```

These are in-place encoding-only measurements: no setup, copying, or input
reset is timed. The benchmark pins CPU 15 and acquires both project benchmark
locks. Never run two benchmarks concurrently. The scripts alternate backend
order and use three processes per configuration, each with three warmups and
101 measured calls. `summary.py` checks output hashes, storage sizes, and medians.
Raw measurements are ignored by Git.

### First integration snapshot (before K16 specialization, 2026-09-19)

On the Ryzen 7950X host, GCC 15.2, Release, `SPIN_TUNE=znver4`, CPU 15:

| K | Forced AVX2 fallback | Auto-selected AVX-512 | Latency reduction |
|---|---:|---:|---:|
| 2^16 | 0.521553 ms | 0.405266 ms | 22.30% |
| 2^18 | 2.191320 ms | 1.687280 ms | 23.00% |
| 2^20 | 10.464042 ms | 9.118250 ms | 12.86% |

Each entry is the median of three process medians under the protocol above.
These comparisons use the same integrated binary and different backend choices.
The earlier native-only four-row experiment measured 9.089597 ms at K=2^20
in separate runs. The integrated result preserves that performance to about 0.3%.
This does not establish a cross-host performance guarantee.

Generic tuning also improves latency: AVX2/AVX-512 pairs are
0.513499/0.412239 ms, 2.122752/1.725231 ms, and 10.733706/9.425354 ms.
At K=2^20, both backends retain 12,713,992 setup bytes after Packed24 compaction
and use a 40 MiB workspace. Setup does not allocate a second offset schedule
for the fast backend.

All eight release test invocations (two tests across four builds) and both
sanitizer tests passed. Generated-source hashes match the measured sources;
all equivalent benchmark outputs have identical hashes. The object-code check
finds EVEX instructions only in the two dispatched AVX-512 objects.

The first integration kept the entire inner/routing loop at AVX2 and lost
part of the experimental gain. Separately compiling the same hot loop for
AVX-512F/VL recovered it. This is an instruction-selection change, not an
additional inner design or a different routing distribution.

## Bidirectional target

Set `SPIN_BIDIRECTIONAL_SOURCE` to the upstream `hypercat/native/spin` directory
to build `spin_half_bidirectional` as well. The author-side source used here is
`C:/Users/peter/repo/hypercat/hypercat/native/spin`. Pass the corresponding Linux
path when building on Linux. The upstream checkout is read, not modified.

`bidirectional.py` generates the patched source in the build directory's
`bidirectional/` subdirectory and records the imported source hashes. Keep
the upstream MIT license when redistributing those derived files. This work
does not regenerate or replace the submitted supplemental archive.

Link either `spin_half_bidirectional` or `spin_half_transpose`, not both.
The bidirectional target preserves `forward`, `forwardBits`, and `wideView`.
It adds the optional backend argument and the checked `encodeInplace` method.
AVX-512 applies to its T128S19 and new T64S12 transposed encoders; other existing
configurations retain the original kernel. One-row tiles are rejected because
the original paired BCH kernel requires at least two rows.

The original row-major offset tables remain intact for forward and wide
consumers. The AVX-512 transpose gets a separate packed offset table. After
Packed24 compaction, that table costs approximately 3N additional bytes at
larger sizes. Tests compare `WideView` arrays byte-for-byte against the fallback.
They also check both dense references, the adjoint identity, in-place encoding,
and packed-bit forward encoding after compaction. No forward speedup is claimed.

## Certified K16 inner option

Use `Configuration::T64S12` with message exponent 16 in either library:

```cpp
bare_spin::Spin code(bare_spin::Configuration::T64S12, 16);
bare_spin::Spin::Workspace workspace(code);
```

This selects BCH [256,128], the same randomized bit-transpose routing, and IMT
with (t,s)=(64,12). The fixed expansion map is a [64,12,24] subcode. Feedback uses
a separate fixed map, and each update uses one sampled transvection.
The recurrence starts at zero, emits before updating, and has no flush.

The full certificate bounds setup failure by 2^-45.38195 at relative distance
10%, K=65536, and N=131072. It covers every outer occupancy q=1,...,512;
this is not merely a one-active-row estimate. See the
[certificate and map-selection record](../inner_design/k16_parameters_20260919/SUBSPACE_RESULT.md).
The transposed library also accepts other aligned lengths, without extrapolating
this certificate. The bidirectional library retains its previous restrictions.
Existing `T128S19` calls keep their parameter choice; selection is not automatic.

`k16_map.py` contains the fixed forward expansion and feedback columns.
It generates symbolically checked XOR circuits without numerical receipts or
scientific Python packages. `k16.py` integrates the option into the generated
sources. It uses a separate shared-XOR reverse kernel, retaining direct K16
routing and four-row BCH processing when AVX-512 is available.
AVX2, two-row tiles, both routing layouts, and compaction remain supported.

The bidirectional target also supports `forward`, `forwardReference`, and
`forwardBits` after compaction. Its dense references and adjoint tests use the
forward convention y=x+Aq, q'=Fq+Bx. The transpose-only adapter explicitly
reverses the historical map names. `wideView()` remains S19-only and rejects
T64S12; its consumers must not interpret a 12-bit state as a 19-bit state.

Run `test_k16_map.py` for the exact expansion spectrum, feedback rank, and
symbolic circuit checks. With local proof receipts, add
`--proof workstreams/inner_design/k16_parameters_20260919/measurements/sub12_full.json`
to compare the fixed columns to the assembled certificate.
An optional `--header BUILD/generated/Map64S12.h` also checks the generated
header byte-for-byte. These checks bind the maps; they do not rerun the bound's
outward-rounding proof producers.

`k16_check.sh ROOT` checks tiled, AVX2-only, detection-suppressed, and sanitizer
builds before running any benchmarks. The last benchmark argument selects the
inner; omitting it retains T128S19:

```sh
out/spin-optimized/spin_benchmark 16 auto 101 1 0 0 6412
out/spin-optimized/spin_bidirectional_benchmark 16 auto 101 1 0 0 6412
```

### Integration measurements (2026-09-19)

The matched integrated builds retain the earlier experiment's improvement:

| Library | Route seed | Existing (128,19) | New (64,12) | Latency reduction |
|---|---:|---:|---:|---:|
| Transpose-only | 1 | 0.341107 ms | 0.330327 ms | 3.16% |
| Transpose-only | 17 | 0.345094 ms | 0.328944 ms | 4.68% |
| Bidirectional, transpose call | 1 | 0.343310 ms | 0.331359 ms | 3.48% |
| Bidirectional, transpose call | 17 | 0.344563 ms | 0.331549 ms | 3.78% |

Peach Ryzen 7950X, CPU 15, GCC 15.2, Release, znver4 tuning, AVX-512 backend.
Each entry is the median of three process medians, each with three warmups
and 101 in-place calls. Setup, copying, and input initialization are excluded.
Configuration order reverses in the middle repetition. All benchmarks run serially.
`k16_summary.py` verifies matching output hashes between the two libraries.
These timings measure transpose calls, not forward performance.

The new option retains 540,672 setup bytes in the transpose-only library and
1,867,784 bytes in the bidirectional library after Packed24 compaction.
The corresponding old values are 532,480 and 1,851,400 bytes.
Both options use a 3 MiB workspace. The bidirectional library retains routing
needed by forward and packed-bit encoding.

All five C++ tests passed in each of five builds: direct, tiled, AVX2-only,
detection-suppressed, and ASan/UBSan. The tests cover both old and new options,
dense forward/transpose references, adjoint equality, both layouts, multiple
tiles, two route seeds, compaction, and in-place suffix preservation.
The exact-map Python suite also passed. The object-code audit found no EVEX
instructions outside the dispatched objects. Raw receipts remain ignored.

## K=2^16 specialization

`SPIN_K16_DIRECT=ON` is the default. It applies only when the selected backend
is AVX-512 and K is exactly 2^16. Larger sizes and AVX2 use their existing path.
Set `SPIN_K16_DIRECT=OFF` to retain the tiled control.

For this size, the whole routed state occupies 2 MiB. A precomputed 32-bit
destination table lets the inner write directly into the four-row BCH layout.
The BCH kernel then reads that state sequentially. This removes the intermediate
bucket-to-tile copy without changing the permutation or inner recurrence.
The small routine is separate from the larger-size hot loop.

Compaction discards transpose tables that the small routine no longer uses.
The transpose-only target retains 532,480 setup bytes, down from 794,632.
The bidirectional target retains 1,851,400 bytes, up from the tiled port's
1,720,332, because it must preserve its forward and wide schedules. Workspace
geometry is unchanged: the existing 3 MiB allocation includes a tile that this
small transposed path does not use.

The tuning compared direct 32-bit and packed24 indices, store-prefetch leads
16 and 64, and tile sizes 32 through 512 with both schedule formats.
Plain 32-bit direct routing was fastest across the tested route seeds.
The numerical comparison is in [K16_RESULTS.md](K16_RESULTS.md).
`small_screen.sh` and `small_confirm.sh` replay the research variants;
`final_check.sh` validates the selected default, fallback, feature masking,
and sanitizer builds before running serial confirmation benchmarks.

Quarter-rate encoding remains unchanged. A useful next step is to adopt the
generated bidirectional port in its upstream application checkout and measure
end-to-end impact; forward-kernel optimization is separate work.
