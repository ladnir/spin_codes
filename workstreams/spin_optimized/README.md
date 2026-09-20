# Optimized SPIN transpose

The supported target is `spin_half_transpose`: a half-rate encoder with a
BCH [256,128] outer, randomized bit-transpose permutation, and IMT inner.
It maps 2K input elements to K output elements. Use the optimized 128-bit API
by default, or the explicit generic XOR-element fallback for other types.

This workstream does not replace the submitted paper artifact or its source pins.
Forward and wide-element optimization are parked, not part of the default build.

## Build and use

Build from the repository root on Linux x86-64 with GCC, CMake 3.20+, and Python 3:

```sh
cmake -S workstreams/spin_optimized -B out/spin-optimized \
  -DCMAKE_BUILD_TYPE=Release -DSPIN_BCH_AVX512=ON
cmake --build out/spin-optimized -j2
ctest --test-dir out/spin-optimized --output-on-failure -j1
```

Link `spin_half_transpose` and include the generated `Spin.h`:

```cpp
#include "Spin.h"
using namespace bare_spin;

Spin code(Configuration::T128S19, MessageLength{81920});
code.compact();
Spin::Workspace workspace(code);
std::vector<block> input(code.codeBlocks()), output(code.messageBlocks());
// Populate input before encoding.
code.encode(input.data(), input.size(), output.data(), output.size(), workspace);
// Or use encodeInplace() on a 2K-element buffer; the suffix beyond K is unchanged.
```

`Spin(Configuration::T128S19, 20)` still means K=2^20. `MessageLength{K}` takes
an actual length. The caller selects the parameter set and handles certification;
the implementation checks structural validity and does not round or retune K.

| Configuration | (t, s, rounds) | K alignment |
|---|---|---:|
| `T128S19` | (128, 19, 1) | 16384 |
| `T64S12` | (64, 12, 1) | 8192 |
| `T64S12R2` | (64, 12, 2) | 8192 |

Aligned lengths through 2^26 are supported. The supplied maps do not support
K=4096. See [ALIGNED_LENGTHS.md](ALIGNED_LENGTHS.md) for geometry, layout selection,
and range-based routing. The selected K16 two-round point and its certificate
are documented in [R2_RESULTS.md](R2_RESULTS.md).

For other element types, include `GenericSpin.h` and call `code.generic()`.
That owned plan preserves the same realized code, including after compaction.
See [GENERIC_ELEMENTS.md](GENERIC_ELEMENTS.md) for the typed workspace interface.

## Supported build options

| Option | Default | Purpose |
|---|---|---|
| `SPIN_BCH_AVX512` | OFF | Compile runtime-dispatched four-row BCH kernels |
| `SPIN_TUNE` | `generic` | Compiler scheduling target; `znver4` was measured |
| `SPIN_BUILD_BENCHMARK` | OFF | Build the serial encoding-only benchmark |
| `BUILD_TESTING` | ON | Build correctness tests |

The baseline requires x86-64-v3, PCLMULQDQ, and VPCLMULQDQ. Enabling AVX-512
does not require every caller to support it: setup dispatches to the available
backend. An explicitly requested unavailable backend throws. GCC 15.2 was tested;
Clang is accepted but unvalidated here. MSVC and ARM are not supported.

Generated sources stay in the build tree. Source dependencies are tracked in
this repository; the default build needs no proof receipts, measurement files,
external checkout, or scientific Python packages. Rebuild callers after changing
generated headers: this is not an ABI-compatible replacement for older builds.

## Tests and measurements

`ctest` runs dense-reference, backend, aligned-length, generic-element, and
selected-map tests. `generic_check.sh ROOT` adds AVX2-only and sanitizer builds.
`lengths_check.sh ROOT` also exercises the large-index boundary and K=2^26;
those large tests require substantial memory. Both scripts run benchmarks only
after builds and correctness checks finish. `isa_check.py BUILD` audits ISA isolation.
`cleanup_check.sh CLEAN_ROOT [UPSTREAM]` validates a clean source archive and,
when supplied, the optional upstream import. It also checks benchmark metadata.

For an encoding-only benchmark, add `-DSPIN_BUILD_BENCHMARK=ON`, then run:

```sh
out/spin-optimized/spin_benchmark 20 auto 101 1
# m backend trials seed tile layout configuration [actual K override]
out/spin-optimized/spin_benchmark 20 auto 101 1 0 0 12819 81920
```

The harness pins CPU 15 and acquires the shared benchmark locks. Run benchmarks
serially. Setup, copying, and input reset are excluded. The JSON field `K` gives
the actual length; `m` is its base-two exponent, or null for other lengths.
Recent timings are in [ALIGNED_LENGTHS.md](ALIGNED_LENGTHS.md) and
[GENERIC_ELEMENTS.md](GENERIC_ELEMENTS.md). Raw results remain ignored by Git.

## Source layout and parked work

- `generate.py`, `k16*.py`, `k18.py`, `lengths.py`: generate the optimized implementation and dispatch.
- `generic.py`, `GenericSpin.h`: generate and expose the XOR-element fallback.
- `*_test.cpp`, `test_k16_map.py`: correctness tests; `*_check.sh` adds broader build matrices.
- `*_screen.sh`, `*_summary.py`, `*_RESULTS.md`: tuning scripts and measured conclusions, not required dependencies for callers.
- [INTEGRATION_HISTORY.md](INTEGRATION_HISTORY.md): earlier integration notes, options, and measurements, preserved without deleting results.
- `cmake/OptionalForward.cmake`, `bidirectional.py`, `wide.*`: opt-in imports from a caller-supplied upstream checkout. These keep their older length interfaces.
- [WIDE_FORWARD.md](WIDE_FORWARD.md), [WIDE_TUNING.md](WIDE_TUNING.md): parked forward work. `SPIN_FORWARD_FOUR` remains experimental and OFF.

Research switches are advanced CMake options. In particular, disabling an
exact-size direct specialization does not disable the newer range-direct path;
also set `SPIN_DIRECT_MAX_K=0` when reproducing a tiled control. Do not link the
transpose-only and bidirectional libraries into the same executable.
