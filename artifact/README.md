# SPIN core implementation artifact

The artifact covers the core SPIN encoder implementation, its build instructions,
and correctness tests. It does not package the numerical proof searches,
certificate receipts, parameter sweeps, application integrations, comparison
benchmarks, or raw measurements used during the research.

The repository is a broader author workspace. Historical reproduction scripts
and research notes remain available, but they are not artifact deliverables.
Their commands may require local inputs that are not distributed.

## Core source

The selected half-rate encoder uses the BCH [256,128] outer, randomized bit
transpose permutation, and IMT inner with (t,s) = (128,19), weight-five feedback.

- [Selected encoder and inner](../workstreams/inner_design/asymmetric/bch256/weight5/implementation/):
  `Weight5Spin.cpp`, `Weight5Inner.h`, and `AsymmetricMap.h`.
- [Shared interface and routing](../workstreams/bare_bch_rm2sub/):
  shared support code and the generated BCH circuit. The directory name is historical.
- [Generated quarter-rate outer circuits](../workstreams/rate_quarter_bch/implementation/generated/):
  an existing build dependency.
- [Selected finite configurations](../workstreams/inner_design/finite_migration/PAPER_RESULTS.md):
  author-side source mapping for the half-rate and quarter-rate variants.

These paths are retained to avoid duplicating or relocating kernels and their
dependencies. Do not substitute a historical baseline merely because its path
has a shorter name.

## Build and test the half-rate encoder

From the repository root on Linux with CMake 3.20+ and a C++20 GCC-compatible
compiler supporting `-march=znver4`:

```sh
cmake -S workstreams/inner_design/asymmetric/bch256/weight5/implementation -B out/spin-core -DCMAKE_BUILD_TYPE=Release
cmake --build out/spin-core --target sparse_pages sparse_pages_test -j 3
ctest --test-dir out/spin-core -R '^sparse_pages_test$' --output-on-failure
```

The current build targets Zen 4 and requires compatible x86 SIMD instructions;
it is not a portable-binary configuration. The `sparse_pages` target selects
the optimized weight-five encoder. Its correctness test demonstrates setup,
workspace use, and encoding through the shared interface. No numerical receipt
or benchmark log is required to build this target.

A clean standalone package and matching quarter-rate build entry point can be
prepared from these sources; this guide does not claim that such a release has
already been assembled. Never run two benchmarks concurrently.

## Author-side material

[PAPER_MAP.md](PAPER_MAP.md), [REPRODUCING.md](REPRODUCING.md), and the application
measurement notes describe the broader research workspace. They are retained
for provenance and maintenance, not as promises to reproduce all paper results.
Missing experiment data is not a blocker for the core-code artifact.
