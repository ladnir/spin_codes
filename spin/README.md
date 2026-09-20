# SPIN encoding library

This standalone C++20 library provides the optimized half-rate BCH [256,128]
SPIN encoder. It has no build dependency on libOTe, Hypercat, Python, or the
research workstreams. SPIN is [MIT licensed](LICENSE), copyright Peter Rindal.
The checked-in kernels have a [generation record](PROVENANCE.md).

Forward encoding maps K records to 2K records. Transposed encoding maps 2K
records to K records. Both apply the same binary matrix, independently to each
bit of a record. They do not perform extension-field multiplication.

## Build and link

```sh
cmake -S spin -B out/spin -DCMAKE_BUILD_TYPE=Release
cmake --build out/spin -j2
ctest --test-dir out/spin --output-on-failure -j1
cmake --install out/spin --prefix /path/to/install
```

A consumer can use `add_subdirectory(path/to/spin)` or an installed package:

```cmake
find_package(spin CONFIG REQUIRED)
target_link_libraries(my_target PRIVATE spin::spin)
```

The target exports C++20 and its public include directory, not ISA flags or
cryptoTools headers. Public headers contain no SIMD types. The implementation
uses private x86 kernels, with runtime AVX-512 dispatch and an AVX2 baseline.
`SPIN_ENABLE_AVX512=OFF` omits the AVX-512 objects. `SPIN_TUNE=znver4` changes
GCC/Clang instruction scheduling without changing the required ISA.

Linux/GCC and Windows/MSVC are supported. With a Visual Studio generator, use
`--config Release` for building/installing and `-C Release` for CTest.
Clang is accepted but has not been validated in this first package pass.
ARM and CPUs without AVX2 are not supported by `Code` construction.

## Single-stream use

Use your own trivially copyable, 16-byte record type. The buffer addresses must
be 16-byte aligned. No conversion to a library-owned block type is required.

```cpp
#include <spin/Code.h>
#include <array>
#include <vector>

struct alignas(16) Block { std::array<std::uint64_t,2> words; };

spin::Code code({.message_size=81920,
                 .parameters=spin::Parameters::T128S19,
                 .route_seed=1, .inner_seed=2});
auto work=code.make_workspace();
std::vector<Block> message(code.message_size());
std::vector<Block> encoded(code.code_size());

// Populate message before encoding.
code.forward<Block>(message, encoded, work);
code.transpose_inplace<Block>(encoded, work);
// The first K records now contain the transposed result; the suffix is unchanged.
// Transpose is not decoding: applying it does not recover the original message.
```

For separate transpose buffers, use `code.transpose<Block>(input, output, work)`.
Separate input and output ranges must not overlap. The in-place operation requires
exactly 2K records and is supported for 128-bit records only.

`Code` is an immutable shared handle. Copying it shares the plan. Moving or
destroying a handle does not invalidate its workspaces; each workspace retains
ownership of the plan. After moving a handle, use it only for destruction or
assignment. A workspace is move-only and belongs to the plan that created it.
Another independently constructed plan is not interchangeable, even with equal seeds.

Use one workspace per simultaneous call. Encoding does not allocate, resize
scratch, copy whole input buffers, or start threads. The checked interface rejects
bad lengths, overlap, alignment, and mismatched workspaces before encoding.
Invalid arguments throw `std::invalid_argument`; unsupported baseline hardware
throws `std::runtime_error`. Allocation failures propagate from setup/workspace creation.

## Parameters and natural lengths

| Parameter set | Inner (t,s,rounds) | K must be a positive multiple of |
|---|---|---:|
| `T128S19` | (128,19,1) | 16,384 |
| `T64S12` | (64,12,1) | 8,192 |
| `T64S12R2` | (64,12,2) | 8,192 |

`message_alignment` and `valid_message_size` provide allocation-free queries.
The implementation never rounds K or changes the chosen parameter set.
Current 32-bit routing requires K < 2^31; available memory also limits allocations.
There is no additional benchmark-size or certificate-size cap.
The caller remains responsible for certificate coverage at the selected parameters.

`ExecutionOptions` selects a backend and an optional tile size. The default tile
is 256 outer rows, clamped to the code geometry. Use 512 explicitly when desired
for the previously measured large wide workload. These are execution choices,
not changes to the binary code. Compaction and route packing are internal.

The exact-size transpose specializations and the range-direct path through
K=458752 remain available. Full and partial tiles use separate kernel bodies.
Forward direct routing retains its selected threshold of K=2^18.

## Wide records and foreign buffers

For 256- or 512-bit records, create a workspace with `Width::Bits256` or
`Width::Bits512` and call `forward<Record>(...)`. Each record holds two or four
adjacent 128-bit lanes. Input and output use the same coordinate-major layout.
Only 16-byte external alignment is required, including for 512-bit records.

Query `code.supports_forward(width)` before selecting a width. Wide forward
supports T128S19 and T64S12R2. A width request is explicit: dispatch never changes
the application's record layout. 512-bit records require AVX-512F/VL/BW/DQ.

The `forward_bytes`, `transpose_bytes`, and `transpose_inplace_bytes` methods
accept zero-copy byte views. This is the intended boundary for foreign record
types and the future Rust binding. Sizes are in bytes and the workspace fixes
record width. Kernels use compiler-supported aliasing behavior internally;
callers must not reinterpret their records as private kernel types.

`forward_bits` encodes one binary message stored in 64-bit words. Bit j of word
i is coordinate 64*i+j. Its output and caller-owned scratch each contain 2K/64
words; input contains K/64 words. All three ranges must be disjoint.

## Other XOR element types

```cpp
#include <spin/Generic.h>
auto generic=code.generic_transpose();
auto bytes=generic.make_workspace<std::uint8_t>();
std::vector<std::uint8_t> choices(code.code_size());
generic.transpose_inplace<std::uint8_t>(choices, bytes);
```

The generic plan copies the same realized routing and inner samples once during
creation. It owns that data independently of `Code`. Encoding uses compile-time
XOR circuits, without per-element runtime dispatch. Elements must be copyable
and default constructible. By default, `operator^` implements addition.
An optional final, statically typed operation argument can instead supply
characteristic-two addition for types without `operator^`. Copy and addition
must not allocate when allocation-free encoding is required.
Use bytes rather than `vector<bool>`.
Generic forward and wide transpose are not provided in this first release.

## Identity and memory

`descriptor()` returns 40 canonical bytes. A consumer can hash them using its
existing hash implementation. The layout is little-endian:

| Byte offset | Content |
|---:|---|
| 0 | Four-byte tag `SPIN` |
| 4 | 32-bit descriptor/setup version, currently 1 |
| 8 | 32-bit parameter-set identifier |
| 12 | Four reserved zero bytes |
| 16 | 64-bit actual K |
| 24 | 64-bit route seed |
| 32 | 64-bit inner seed |

Version 1 fixes the BCH [256,128] outer, map definitions, setup expansion, and
coordinate conventions in this snapshot. Backend, tile size, and packing do not
affect the descriptor. Any future change to the realized map for these inputs
requires a new version, not a silent implementation update.

`setup_bytes()` and `workspace.bytes()` report allocated vector capacities,
excluding inputs, outputs, allocator overhead, and small object headers.
The initial package prepares both directions and retains routing for bit-packed
forward. Direction-specific setup pruning is a possible later memory optimization.

## Validation and remaining integration

`spin_api_test` checks dense-reference equality, the forward/transpose adjoint
identity, generic byte/64-bit equality, packed-bit equality, and wide lane equality.
It covers natural lengths, partial tiles, exact-size and range-direct paths,
invalid buffers, descriptor fields, workspace ownership, concurrent calls with
separate scratch, and allocation-free successful encoding.
`spin_known_answers` freezes descriptor-v1 forward and transpose outputs for
all three parameter sets, independent of compiler and selected backend.
It also checks that packed-bit encoding agrees with the low-bit projection of
those known-answer outputs. The API test checks the packed BCH lookup separately
from the recursive inner.

The packed encoder's fixed lookup tables must remain `constexpr`. With MSVC
19.51.36256 and `/O2`, the former runtime feedback-table initializer compiled
with an uninitialized loop bound, causing out-of-bounds writes, incorrect output,
or a crash. This was reproduced in an extracted encoder and confirmed from its
faulting instruction and disassembly. Compile-time initialization removes that
code path; it does not disable optimization or change the encoded map.
`tools/probe_packed.py` extracts a small standalone compiler test directly from
the production source. Its optional `--runtime-tables` mode deliberately restores
the failing initialization; the manual `Packed MSVC diagnosis` workflow compares
that mode with the current optimized and unoptimized builds. This diagnostic
does not add a Python dependency to ordinary library builds.

`tests/consumer` is a separate installed-package project. It defines its own
`osuCrypto::block`-named type to check that no imported type leaks into consumers.
Build it with `CMAKE_PREFIX_PATH` pointing at an installed SPIN package.

Validation on 2026-09-20 passed both CTest tests with MSVC 19.50 on Windows,
GCC 15.2 on Linux with AVX-512 BCH and wide kernels enabled, and GCC 13.3
with AVX-512 disabled under AddressSanitizer and UndefinedBehaviorSanitizer.
The separate installed-package consumer also built and ran on Windows and Linux.
These are correctness and integration checks, not performance measurements.

The libOTe adapter now supports coefficient contexts and Silent OT; its build and
tests live in the libOTe checkout. Hypercat's C ABI/Rust migration remains next.
This package pass does not supply new performance measurements
or certificates; previous kernel timing reports remain in the research workstream.
