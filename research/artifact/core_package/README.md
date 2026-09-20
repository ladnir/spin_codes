# SPIN core implementation

This supplement contains encoder source, its dependencies, and correctness
tests. It contains no numerical proof archives, parameter sweeps, comparison
benchmarks, application protocols, raw measurements, or prebuilt binaries.

## Included implementations

| Build target | Direction | Selected configuration |
|---|---|---|
| `spin_half_transpose` | Transposed, N blocks to K blocks; includes in-place API | Rate 1/2, BCH [256,128], IMT (t,s)=(128,19), weight-five feedback |
| `spin_quarter_transpose` | Transposed, N blocks to K blocks; includes in-place API | Rate 1/4, BCH [128,32], IMT (128,19), weight-three feedback |
| `spin_half_bidirectional` | Forward K to 2K, and transposed 2K to K | Rate 1/2, BCH [256,128], IMT (128,19), weight-five feedback |

Each block holds 128 parallel binary instances. The bidirectional target also
exposes `forwardBits` for a single bit-packed binary input. Its source comes
from the application-side encoder, without importing the application itself.
No quarter-rate forward implementation is included.

These are **alternative libraries**, not three components to link together:
they retain the same `bare_spin::Spin` class name and have overlapping symbols.
Link exactly one target into an executable. The kernels and generated circuits
are copied byte-for-byte; only the packaging build and bidirectional test are new.

## Build and check

Requirements: Linux x86-64, CMake 3.20+, a C++20 compiler (GCC 13+ recommended),
and a CPU with AVX2, PCLMULQDQ, and VPCLMULQDQ. A Linux environment under WSL2
also works on compatible hardware. No Python, Rust, network access, or external
cryptoTools checkout is needed to build the extracted package.

From this directory:

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j 2
ctest --test-dir build --output-on-failure -j 1
```

The default is `-march=native`. On a Zen 4 host, use `-DSPIN_ARCH=znver4`
to select the paper's architecture target. Do not run a Zen 4 binary on hardware
that lacks its instructions. This package does not claim portable ARM or MSVC
support. Linux page advice is best-effort and changes no system-wide settings.

Builds may run in parallel; **never run two benchmarks concurrently**.
The supplied tests check correctness, not performance or distance certificates.

## Using the encoder

The short programs in `tests/` serve as complete usage examples. Include
`Spin.h` and select `Configuration::T128S19`. Construct a `Spin`, then a
`Spin::Workspace` for that object. Allocate buffers using `messageBlocks()`
and `codeBlocks()`. Reuse the setup and workspace for subsequent calls.

- Half-rate transposed: default outer; see `tests/half.cpp`.
- Quarter-rate transposed: pass `Outer::Bch128x32`; see `tests/quarter.cpp`.
- Half-rate forward: call `forward`; `encode` means transpose, not forward.
  See `tests/bidirectional.cpp` for both directions and their adjoint check.
- `encodeInplace` in the transposed libraries replaces the first K blocks of
  an N-block buffer, leaving the suffix unchanged. Ordinary forward encoding
  uses separate, nonoverlapping input and output buffers.
- `validateSetup` and dense-reference routines are correctness tools, not hot
  paths. Run them before `compact`, which discards the reference schedules.

Use message exponents 16, 18, or 20 for the operating sizes exercised by the
transposed tests. Shared headers retain historical configuration enums; their
presence does not extend the selected IMT configuration or its paper claims.
Setup seeds in examples are fixed for repeatable testing.

## Layout and source integrity

`src/half/` and `src/routing/deployment/` contain the selected transpose kernels.
`src/common/`, `src/asymmetric/`, and `vendor/` supply their dependencies.
`src/bidirectional/` is a separate, self-contained half-rate implementation.

`SOURCE_MAP.txt` maps copied files to their source-tree paths. `SHA256SUMS`
authenticates every other archive member. Verify on Linux with:

```sh
sha256sum -c SHA256SUMS
```

The bidirectional source retains its upstream MIT license in
`licenses/bidirectional-MIT.txt`. Existing notices in vendored headers are
preserved. No additional license grant is inferred for other source files.
This archive preserves attribution; it is not claimed to be anonymized.
