# SPIN encoder code

This supplement contains only the core encoders, their dependencies, correctness
tests, and build instructions. It contains no ZK or OT protocols, experimental
data, numerical proof archives, benchmark results, or binaries.

## Implementations

- `spin_half_transpose`: rate-1/2 BCH [256,128] with IMT (t,s)=(128,19)
  and weight-five feedback.
- `spin_quarter_transpose`: rate-1/4 BCH [128,32] with IMT (128,19)
  and weight-three feedback.
<!-- BEGIN FORWARD -->
- `spin_half_bidirectional`: half-rate forward and transposed encoding.

The bidirectional library also provides `forwardBits` for one bit-packed binary
input. There is no quarter-rate forward implementation in this package.
<!-- END FORWARD -->

Each 128-bit block holds 128 parallel binary instances. The libraries are
alternatives with overlapping `bare_spin::Spin` symbols: link **exactly one**
into each executable.

## Build and test

Use Linux x86-64, CMake 3.20+, and a C++20 compiler (GCC 13+ recommended).
The CPU must support AVX2, PCLMULQDQ, and VPCLMULQDQ. No downloads, Python,
Rust, or external library checkout are required.

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j 2
ctest --test-dir build --output-on-failure -j 1
```

The default is `-march=native`; use `-DSPIN_ARCH=znver4` on a compatible
Zen 4 host. ARM and MSVC are not supported by this build. Linux page advice
is best-effort and changes no system-wide settings. Never run two benchmarks
concurrently. These tests check correctness, not timing or distance guarantees.

## Usage

Include `Spin.h`, select `Configuration::T128S19`, construct a `Spin`, and then
construct its `Spin::Workspace`. Reuse setup and workspace between calls.
Use `messageBlocks()` and `codeBlocks()` to size input/output buffers.

- In both transpose libraries, `encode` maps N blocks to K blocks.
  `encodeInplace` replaces the first K blocks and preserves the suffix.
- The quarter-rate constructor selects `Outer::Bch128x32`.
<!-- BEGIN FORWARD -->
- In the bidirectional library, `forward` maps K blocks to 2K blocks;
  `encode` is the transpose. Forward input/output buffers must not overlap.
<!-- END FORWARD -->
- Correctness tests in `tests/` are complete usage examples. The transpose
  suites cover message exponents 16, 18, and 20. Other historical configuration
  enums do not extend the selected IMT operating points.
- Dense-reference checks and `validateSetup` are not performance paths.
  Use them before `compact`, which discards reference schedules.

## Integrity

All encoder kernels and generated circuits retain their original bytes.
Packaging removes original source-tree provenance, not implementation code.
Existing third-party notices are preserved. Any included `licenses/` files
apply to their named components, not automatically to the entire package.

```sh
sha256sum -c SHA256SUMS
```
