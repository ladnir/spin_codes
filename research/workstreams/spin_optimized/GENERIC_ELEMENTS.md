# Transposed encoding with XOR elements

Use the optimized `Spin` API for 128-bit elements. For another element type,
include `GenericSpin.h` and create a generic plan from the selected `Spin` object:

```cpp
#include "GenericSpin.h"
using namespace bare_spin;

Spin setup(Configuration::T128S19, MessageLength{81920});
setup.compact();
auto code = setup.generic();
GenericTranspose::Workspace<std::uint64_t> work(code);
std::vector<std::uint64_t> input(code.codeBlocks());  // 2K elements
std::vector<std::uint64_t> output(code.messageBlocks());  // K elements
code.encode<std::uint64_t>(input, output, work);
// Or replace the first K elements of an existing 2K-element buffer:
code.encodeInplace<std::uint64_t>(input, work);
```

The plan copies the realized permutation and IMT samples. It neither resamples
the code nor chooses new parameters. It works before or after `compact()` and
owns its setup, so it remains usable after the original `Spin` object is destroyed.
The caller retains responsibility for certificate coverage at the selected K.

An element represents parallel binary values. Its `operator^` must implement
componentwise XOR, return a value convertible to the element type, and avoid
hidden allocation. The type must support copying and assignment, be trivially
copyable, and be default constructible.
Supported examples include unsigned integers, `std::byte`, `block`, and fixed
arrays wrapped in a type with an XOR operator. `bool` is excluded because
`std::vector<bool>` is not contiguous element storage; use `uint8_t` for individual bits.
This is not arithmetic over an arbitrary odd-characteristic field.

Encoding uses templated XOR circuits and a reusable workspace of 2K elements.
The IMT circuit keeps its compile-time pruned transform and shared XOR terms.
The BCH circuit is checked symbolically against the generator during the build.
There is no allocation, virtual dispatch, type erasure, or lane-packing pass
inside encoding. Element constructors and XOR operators remain the caller's responsibility.

The checked methods reject incorrect geometry and overlap with workspace.
`encode` also rejects overlapping input and output; use `encodeInplace` for that
case. In-place encoding leaves the suffix beyond K unchanged. All input is read
before any output is written. Each concurrent call needs its own workspace;
the encoding plan is immutable during use.

This fallback uses natural row-major routing and one-row BCH circuits. It is
not a replacement for the optimized 128-bit path or a tuned wide-element kernel.
The `Spin` object's layout, encoding methods, and SIMD dispatch are unchanged;
only an explicit call to `generic()` creates the additional setup. Merely
including `Spin.h` does not instantiate the generic circuits.

The current build still requires the existing Linux x86-64 toolchain and baseline
ISA. Generic *element types* do not imply CPU portability. This API belongs to
`spin_half_transpose` with `SPIN_GENERAL_LENGTHS=ON`; the optional bidirectional
library and forward encoder are unchanged.

`generic_test.cpp` compares every byte of several element widths against the
128-bit encoder. It covers all three inner configurations, both compacted index
layouts, both CPU backends, partial tiles, and direct-routing specializations.
It also checks source-object lifetime, in-place suffix preservation, invalid
buffer geometry, workspace aliasing, and absence of encode-time allocations.
`generic_check.sh` builds release, AVX2-only, and ASan/UBSan variants, then runs
the existing 128-bit benchmarks serially to check for regressions.

## Validation result

The release, AVX2-only, and ASan/UBSan suites passed on the Ryzen 9 7950X host
with GCC 15.2. The tested types were `uint8_t`, `std::byte`, `uint32_t`, `uint64_t`,
`block`, and aligned custom 256/512-bit elements. All paired benchmark hashes matched.
ISA inspection found no EVEX instructions in the generic factory or test object.

The existing `Spin.cpp` and `Fast.cpp` object files, and the complete 128-bit
benchmark executable, are byte-for-byte identical to the pre-fallback build.
Three serial process runs, each with 51 timed encodes,
gave these medians for the optimized 128-bit path:

| K | Before (ms) | With generic support available (ms) |
|---:|---:|---:|
| 2^16, T64S12R2 | 0.340566 | 0.342529 |
| 2^18, T128S19 | 1.413428 | 1.418578 |
| 2^20, T128S19 | 9.076624 | 9.150673 |

All differences are below 1%, from runs of identical binaries. These measure the existing optimized path,
not generic-fallback throughput. `generic_summary.py` checks the paired hashes
and summarizes the ignored measurement files.
