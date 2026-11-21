# SPIN encoding library

This standalone C++20 library implements Turbo SPIN
(`spin::Parameters::PacketRs8T64`), a half-rate binary linear code with fast
transposed and forward encoding. It has no build dependency on libOTe,
Hypercat, Python, or the research workstreams. SPIN is [MIT licensed](LICENSE),
copyright Peter Rindal. [RS8.md](RS8.md) describes the construction, backends,
measurements, and distance certificates; [PROVENANCE.md](PROVENANCE.md)
records where the kernels came from.

Forward encoding maps K records to 2K records. Transposed encoding maps 2K
records to K records. Both apply the same binary matrix independently to each
bit of a record. The kernels use finite-field instructions internally to
evaluate that binary matrix; callers supply ordinary fixed-width binary records.

SPIN 0.4 supports only Turbo SPIN. The earlier families (paper SPIN with a
BCH[256,128] outer and IMT inner, identifiers 1-3, and the packet profiles
with identifiers 4-7) were removed together with bit-packed forward encoding
(`forward_bits`) and generic XOR element types (`Generic.h`). Their
identifiers are never reused; the last source tree containing them is commit
`584a7fb6`.

## Build and link

```sh
cmake -S spin -B out/spin -DCMAKE_BUILD_TYPE=Release
cmake --build out/spin -j2
ctest --test-dir out/spin --output-on-failure -j1
cmake --install out/spin --prefix /path/to/install
```

A consumer can use `add_subdirectory(path/to/spin)` or an installed package:

```cmake
find_package(spin 0.4 CONFIG REQUIRED)
target_link_libraries(my_target PRIVATE spin::spin)
```

The target exports C++20 and its public include directory, not ISA flags.
Public headers contain no SIMD types. Kernels are private and selected at
runtime: AVX-512, AVX2+GFNI, and AVX2 on x86-64, NEON on AArch64, and an
optimized portable kernel on any 64-bit CPU. `SPIN_ENABLE_AVX512=OFF` omits
the AVX-512 objects. `SPIN_TUNE=znver4` changes GCC/Clang instruction
scheduling without changing the required ISA. Linux/GCC and Windows/MSVC are
supported; with a Visual Studio generator, use `--config Release` for
building/installing and `-C Release` for CTest.

## Single-stream use

Use your own trivially copyable 16-, 32-, or 64-byte record type. Buffer
addresses must be 16-byte aligned. No conversion to a library-owned type is
required.

```cpp
#include <spin/Code.h>
#include <array>
#include <vector>

struct alignas(16) Block { std::array<std::uint64_t,2> words; };

spin::Code code({.message_size=1u<<18, .route_seed=17, .inner_seed=41});
auto work=code.make_workspace();
std::vector<Block> message(code.message_size());
std::vector<Block> encoded(code.code_size());

// Populate message before encoding.
code.forward<Block>(message, encoded, work);
code.transpose_inplace<Block>(encoded, work);
// The first K records now contain the transposed result; the suffix is unchanged.
// Transpose is not decoding: applying it does not recover the original message.
```

`CodeSpec::parameters` defaults to `PacketRs8T64`, the only supported value.
For separate transpose buffers, use `code.transpose<Block>(input, output, work)`;
the ranges must not overlap. In-place transpose requires exactly 2K 128-bit
records.

`code.make_buffer()` returns a zero-initialized, 64-byte-aligned `Buffer` of
exactly `code_size()` records; its `bytes()` view works with the byte-based
methods. Caller-owned storage is equally supported.

`Code` is an immutable shared handle. Copying it shares the plan. Moving or
destroying a handle does not invalidate its workspaces; each workspace retains
ownership of the plan. A workspace is move-only and belongs to the plan that
created it; another independently constructed plan is not interchangeable,
even with equal seeds. For successive maps, call `next.prepare_workspace(work)`
before encoding with `next`: scratch of the same size and width is rebound
without allocation or clearing; other workspaces are recreated.

Use one workspace per simultaneous call. Encoding does not allocate, resize
scratch, copy whole input buffers, or start threads. The checked interface
rejects bad lengths, overlap, alignment, and mismatched workspaces before
encoding. Invalid arguments throw `std::invalid_argument`; an explicitly
requested backend that the CPU or build lacks throws `std::runtime_error`.

## Sizes, backends and capabilities

K may be any positive multiple of 128 up to 132,152,832
(`message_alignment` returns 128; `valid_message_size` is the allocation-free
check). The library never rounds K. Accepting a size is not a distance claim:
only the sizes listed in [RS8.md](RS8.md#certificate-scope) have certificates.

`ExecutionOptions::backend` selects an implementation of the same seeded map:
`Automatic` (default) picks the first available of `Avx512`, `Avx2Gfni`,
`Avx2`, `Neon`, `Portable`; `Reference` is the literal oracle. Backend choice
never changes the map or descriptor. `tile_rows` accepts only 0 or 8.

`capabilities()` reports which optimized kernels this CPU and build can run
(`avx2`, `avx2gfni`, `avx512`). Every record width is supported on every CPU
through the portable kernel; `Code::supports_forward(width)` and
`supports_transpose(width)` answer for a particular code.

| Operation | Record widths |
|---|---|
| Forward (K -> 2K) | 128, 256, 512 bits |
| Transpose (2K -> K), separate or in place | 128 bits |

## Wide records, output stores and foreign buffers

For 256- or 512-bit records, create a workspace with `Width::Bits256` or
`Width::Bits512` and call `forward<Record>(...)`. Each record holds two or
four adjacent 128-bit lanes; input and output use the same coordinate-major
layout. Only 16-byte external alignment is required.

Forward calls take an output-store preference:

```cpp
// Stream a large codeword that will be consumed after a batch of encodings.
code.forward_bytes(message, encoded, workspace, spin::OutputStores::Streaming);
// Keep an immediately consumed result in cache.
code.forward_bytes(message, encoded, workspace, spin::OutputStores::Cached);
```

`Streaming` (the default) requests non-temporal output stores where the
kernel and alignment allow them, and ordinary stores otherwise. `Cached`
forces ordinary stores. The choice is per call and changes neither the
encoded bytes nor the descriptor; streaming writes are fenced before return.
Measure the complete consumer workload, including subsequent reads, when
choosing.

`forward_bytes`, `transpose_bytes`, and `transpose_inplace_bytes` accept
zero-copy byte views; this is the boundary for foreign record types and
language bindings. Sizes are in bytes and the workspace fixes the record width.

## Bit-vector transposed encoding

A 128-bit record transpose is 128 independent copies of one binary map, one
per bit lane. `transpose_bits` applies that map to a single bit per
position, for example a choice-bit vector, at a fraction of the record cost
(one Ryzen 9 7950X core: 0.15 ms against 0.47 ms at K = 2^18, 0.62 ms
against 2.90 ms at K = 2^20):

```cpp
auto bits=code.make_bit_workspace();
std::vector<std::uint64_t> in(code.code_size()/64), out(code.message_size()/64);
code.transpose_bits(in, out, bits); // 2K bits -> K bits
```

Bit i is bit `i%64` of word `i/64`. For every lane j the result equals bit j
of `transpose`. The first bit workspace builds the plan's bit tables (routing
inverse and coefficient lanes, about 0.7 MiB at K = 2^18); later workspaces
share them. `PreparedEncoder` offers the same calls; with `BankedHeuristic`
setup, `setCodeSeed` refreshes the bit tables in place. AVX-512 runs a GFNI
kernel; every other backend runs a portable kernel with the same output.

## Repeated transposed encoding and seed updates

For successive code seeds, use the move-only `PreparedEncoder` and choose the
setup mode when constructing it:

```cpp
#include <spin/PreparedEncoder.h>

spin::PreparedEncoder encoder(
    {.message_size=1u<<18},
    {.mode=spin::SetupMode::BankedHeuristic, .bank_seed=913});
auto work=encoder.make_workspace();

encoder.setCodeSeed({.route=17, .inner=29});
encoder.transpose_inplace<Block>(encoded, work); // Exactly 2K records.
```

| Mode | Seed update | Reused across updates |
|---|---|---|
| `Full` (default) | Rebuild the exact seeded setup | Compatible workspace allocations |
| `BankedHeuristic` | Reorder the bank's routing regions and inner updates; fresh outer scalars | Bank and all workspace allocations |

`setCodeSeed` is a no-op for identical seeds. Banked refresh allocates
nothing and does O(K) copying. `bank_seed` fixes the bank; replace the encoder
to change the bank, mode, or geometry. Both parties must agree on these
settings and the current seeds. After a `Full` refresh,
`encoder.prepare_workspace(work)` rebinds scratch and releases the previous
setup outside the timed call; encoding also does this automatically.

Banked setup samples a correlated family (one bank reused across codes) that
the distance certificates do not cover; use `Full` where the certificate is
required. Refresh requires exclusive access to the encoder; encoding permits
concurrent calls with distinct workspaces.

## Identity and memory

`Code::descriptor()` returns 40 canonical little-endian bytes for hashing:

| Byte offset | Content |
|---:|---|
| 0 | Four-byte tag `SPIN` |
| 4 | 32-bit descriptor version, 1 |
| 8 | 32-bit parameter identifier, 8 |
| 12 | Four reserved zero bytes |
| 16 | 64-bit K |
| 24 | 64-bit route seed |
| 32 | 64-bit inner seed |

`PreparedEncoder::descriptor()` returns 56 bytes: the same fields with version
2 at offset 4 and the setup mode at 12, then the bank seed at 40 and the
bank-family revision (1) at 48; both are zero for `Full`. Backend, tile size,
and store policy never affect a descriptor. A changed realized map requires a
new identifier or version, not a silent implementation update; these
descriptors are unchanged from SPIN 0.3.

`setup_bytes()` and `workspace.bytes()` report owned storage capacities.
Owned buffers and scratch default to `MemoryPolicy::Automatic`: allocations of
at least 16 MiB use 2 MiB-aligned storage with best-effort Linux huge-page
advice; smaller ones keep 64-byte alignment. `Normal` and `PreferHugePages`
override this. Preparation never advises caller-owned buffers.

## Validation

`spin_rs8_test` checks every backend against the literal oracle (CTests
`spin_rs8`, `spin_rs8_avx2`, `spin_rs8_gfni`, `spin_rs8_portable`,
`spin_rs8_neon`; each skips where its ISA is absent): fixed known answers,
natural sizes and overflow bounds, forward/transpose adjoints at every width,
in-place suffix preservation, buffer guards, descriptors, retired identifiers,
workspace ownership and rebinding, allocation-free encoding, and both prepared
modes. `spin_turbo_setup_test` checks that every setup path produces
bit-identical tables. `spin_bits_test` checks `transpose_bits` against bit
lanes 0, 1, 63, 64 and 127 of the record transpose on every backend the CPU
runs, at natural sizes from 128 to 65,536, on dense and sparse inputs, for
`Code` and both prepared modes across reseeds, with allocation-free calls.
`tests/consumer` is a separate installed-package
project that checks no imported type leaks into consumers.

`SPIN_BUILD_BENCHMARKS=ON` builds `spin_turbo_backends_bench`,
`spin_setup_once_bench`, and (with AVX-512) `spin_rs8_stages_bench`. Run
benchmarks one at a time on an otherwise idle core.
