# SPIN encoding library

This standalone C++20 library provides two half-rate binary-code families:
paper SPIN with a BCH [256,128] outer and IMT inner, and the newer RS packet
construction. It has no build dependency on libOTe, Hypercat, Python, or the
research workstreams. SPIN is [MIT licensed](LICENSE), copyright Peter Rindal.
The checked-in kernels have a [generation record](PROVENANCE.md).

The [wide encoding notes](WIDE_ENCODING.md) describe composed routing, the
runtime-dispatched register kernels, regeneration, and portability tests.

Forward encoding maps K records to 2K records. Transposed encoding maps 2K
records to K records. Both apply the same binary matrix, independently to each
bit of a record. The RS packet kernels use finite-field instructions internally
to evaluate that binary matrix; callers still supply ordinary 128-bit records.

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
uses private x86 kernels with runtime AVX-512 dispatch. The IMT families require
AVX2; the packet family also has an SSE2 fallback.
`SPIN_ENABLE_AVX512=OFF` omits the AVX-512 objects. `SPIN_TUNE=znver4` changes
GCC/Clang instruction scheduling without changing the required ISA.

Linux/GCC and Windows/MSVC are supported. With a Visual Studio generator, use
`--config Release` for building/installing and `-C Release` for CTest.
Clang is accepted but has not been validated in this first package pass.
ARM is not supported. On x86-64 CPUs without AVX2, only the packet family is supported.

## Single-stream use

All supplied families use `CodeSpec`, `Code`, and `Workspace`. Changing the
family changes the specification, not the sizing or transposed-encoding calls.

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

Select `Parameters::PacketRsT64S20` for the newer packet construction. It supports
the same transpose calls, but not forward encoding or generic XOR elements yet.
Query `supports_forward(width)`, `supports_transpose(width)`, and
`supports_generic_transpose()` before using an optional operation.
Unsupported operations throw; they never substitute another code family.

For library-owned storage, `code.make_buffer()` returns a zero-initialized,
64-byte-aligned `Buffer` holding exactly `code_size()` records. Its `bytes()` view
works directly with the byte-based encoding methods. Caller-owned storage remains
supported; no copy into library storage is required.

`Code` is an immutable shared handle. Copying it shares the plan. Moving or
destroying a handle does not invalidate its workspaces; each workspace retains
ownership of the plan. After moving a handle, use it only for destruction or
assignment. A workspace is move-only and belongs to the plan that created it.
Another independently constructed plan is not interchangeable, even with equal seeds.

For successive maps, call `next.prepare_workspace(work)` before encoding with
`next`. Compatible 128-bit scratch is rebound without allocation or clearing.
The previous plan can no longer use that workspace. Different geometry and wide
workspaces are recreated when necessary; a moved-from workspace is recreated at
128 bits. Do not rebind while another call uses the workspace. This avoids
repeated scratch allocation when a protocol changes only its code seeds.

Setup and encoding are separate costs. The [setup benchmark report](SETUP_PERFORMANCE.md)
records both, including repeated batches with changed seeds.

Opt-in [Feistel routing experiments](experiments/feistel/README.md) compare
compact permutation families with the exact-shuffle baseline. They are not
public code profiles and do not inherit the existing distance certificates.

Use one workspace per simultaneous call. Encoding does not allocate, resize
scratch, copy whole input buffers, or start threads. The checked interface rejects
bad lengths, overlap, alignment, and mismatched workspaces before encoding.
Invalid arguments throw `std::invalid_argument`; unsupported baseline hardware
throws `std::runtime_error`. Allocation failures propagate from setup/workspace creation.

## Repeated transposed encoding and seed updates

`Code` remains immutable. For successive code seeds, use the move-only
`PreparedEncoder` controller and choose the setup mode when constructing it:

```cpp
#include <spin/PreparedEncoder.h>

spin::PreparedEncoder encoder(
    {.message_size=81920, .parameters=spin::Parameters::T128S19},
    {.mode=spin::SetupMode::BankedHeuristic, .bank_seed=913});
auto work=encoder.make_workspace();

encoder.setCodeSeed({.route=17, .inner=29});
encoder.transpose_inplace<Block>(encoded, work); // Exactly 2K records.
```

The constructor prepares the initial seeds from `CodeSpec`. `setCodeSeed` is a
no-op for identical seeds. Its behavior for changed seeds depends on the fixed mode:

| Mode | Seed update | Reused across updates |
|---|---|---|
| `Full` (default) | Rebuild the selected family's sampled setup | Compatible workspace allocations |
| `BankedHeuristic` | Refresh row transformations, region labels, offsets, and IMT masks | Bank and all workspace allocations |

Banked refresh allocates nothing and never constructs a full per-round route.
The optimized transpose generates addresses for one region at a time. It supports
the three IMT parameter sets and their natural lengths, including non-powers of two.
The packet family supports `Full` only; `BankedHeuristic` is rejected for that family.
`bank_seed` determines the reusable bank; code seeds do not change it. Replace the
controller explicitly to change the bank, mode, or geometry. Both parties must agree
on these settings and the current code seeds.

After a `Full` refresh, `encoder.prepare_workspace(work)` rebinds compatible
scratch and releases its previous setup before the timed encoding call. Encoding
also performs this check automatically; omitting explicit preparation can defer
that deallocation until the first call. `BankedHeuristic` keeps its existing setup.

`generic_transpose()` returns a live view that follows the controller's seed updates.
Its typed workspace remains valid across updates. In contrast, a generic view from
immutable `Code` remains a fixed map. Create views and workspaces during preparation.
Do not refresh or create a generic view concurrently with another use of the controller
or its views. Encoding permits concurrent calls with distinct workspaces.
Workspaces and views retain their required setup; there is no global cache.
Use moved-from controllers only for assignment or destruction.

The prepared interface exposes transpose only, with the same capability queries,
size queries, buffer allocation, and width argument as `Code`. Generic XOR elements
remain available for the IMT families only. The immutable interface retains forward
and wide support for the families listed below.
The banked path uses a region buffer, so `tile_rows` does not affect that path.

Prepared descriptors contain 56 little-endian bytes: magic at 0, version 2 at 4,
parameter set at 8, mode at 12, K at 16, route seed at 24, and inner seed at 32.
The bank seed is at 40 and bank-family revision 1 is at 48; both are zero for `Full`.
The first three numeric fields are 32-bit; K and subsequent fields are 64-bit.
`Code::descriptor()` retains its original 40-byte version-1 format.

`BankedHeuristic` is the generalized `row-rotate1` experiment, not the original
uniform-permutation ensemble. Its distance certificates do not follow from the
original analysis. The K=2^18 regression compares complete routes and IMT masks
against that prototype. The integrated K=2^18 fast path measures approximately
8--10% fresh-code overhead on the tested Ryzen/Linux host. This is not an API-wide
performance guarantee; other lengths retain the generalized routing path.
See the [integration checkpoint](PREPARED_ENCODER.md) for validation and current timings.

## Parameters and natural lengths

| Parameter set | Inner (t,s,rounds) | K must be a positive multiple of |
|---|---|---:|
| `T128S19` | (128,19,1) | 16,384 |
| `T64S12` | (64,12,1) | 8,192 |
| `T64S12R2` | (64,12,2) | 8,192 |
| `PacketRsT64S20` | (64,20), GL20 updates | 256 |

`message_alignment` and `valid_message_size` provide allocation-free queries.
The implementation never rounds K or changes the chosen parameter set.
Current 32-bit routing requires K < 2^31. RS packet routing stores padded groups,
so its representation bound is slightly smaller. Query `valid_message_size`
rather than hard-coding that bound. Available memory also limits allocations.
There is no additional benchmark-size or certificate-size cap.
The caller remains responsible for certificate coverage at the selected parameters.

`ExecutionOptions` selects a backend and an optional tile size. The default tile
for IMT is 256 outer rows, clamped to the code geometry. Use 512 explicitly when desired
for the previously measured large wide workload. These are execution choices,
not changes to the binary code. Compaction and route packing are internal.
RS packet encoding uses eight-row tiles and accepts only `tile_rows=0` or `8`.

| Family | 128-bit transpose | 128-bit forward | Wide forward | Generic XOR transpose |
|---|---|---|---|---|
| `T128S19` | Yes | Yes | 256/512 bits | Yes |
| `T64S12` | Yes | Yes | No | Yes |
| `T64S12R2` | Yes | Yes | 256/512 bits | Yes |
| `PacketRsT64S20` | Yes | No | No | No |

Wide support also depends on the CPU and compiled backends.

For paper SPIN, the exact-size transpose specializations and the range-direct path through
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

For an IMT family, the generic plan copies the same realized routing and inner samples once during
creation. It owns that data independently of `Code`. Encoding uses compile-time
XOR circuits, without per-element runtime dispatch. Elements must be copyable
and default constructible. By default, `operator^` implements addition.
An optional final, statically typed operation argument can instead supply
characteristic-two addition for types without `operator^`. Copy and addition
must not allocate when allocation-free encoding is required.
Use bytes rather than `vector<bool>`.
Generic forward and wide transpose are not provided in this first release.

## Precomputed RS packet transpose

`Parameters::PacketRsT64S20` selects eight parallel GF16 RS [16,8] rows per
group, with 256 binary inputs and 512 binary outputs. Independent GF(2^32)
randomizers act on the aligned 32-bit outer symbols. Four-bit packets then pass
through structured routing and a t64/s20 inner with independent GL20 updates.
The forward construction uses the binary adjoints of the field multipliers;
the transposed kernel evaluates the field multipliers themselves.

This family currently supports transposed encoding of 128-bit records only.
Paper SPIN remains available for forward encoding, wide records, and generic XOR types.

```cpp
#include <spin/Code.h>

spin::Code packet({.message_size=1u<<20,
                   .parameters=spin::Parameters::PacketRsT64S20,
                   .route_seed=1, .inner_seed=1});
auto scratch=packet.make_workspace(); // Prepare once, outside timing.
auto buffer=packet.make_buffer();
// Fill buffer.bytes() with 2K records, then encode repeatedly:
packet.transpose_inplace_bytes(buffer.bytes(), scratch);
```

K must be a positive multiple of 256. Non-powers of two use the same kernel;
construction preserves K without implicit rounding. Transpose consumes exactly
2K records and produces K records. In-place calls preserve the last K records.
Separate buffers must not overlap and need 16-byte alignment. A 64-byte-aligned
output enables the measured non-temporal-store path; other supported alignments
use ordinary stores for the same map. Encoding allocates nothing.

The current numerical certificate covers K=2^20, rate 1/2, relative distance
above 10%, and a 68.103757-bit setup-failure margin under independent ideal setup.
It covers all 4,096 outer-group occupancies. Other accepted lengths are
implementation support, not additional certificate claims.
The [construction and replay instructions](../research/workstreams/k16_design/README.md)
and [standalone replay](../research/workstreams/k16_design/reproduce_rs16_k20.py)
preserve the proof and its frozen source dependencies.

Setup is immutable and sampled once. The implementation retains the research
prototype's deterministic seed expansion. `route_seed` generates the routing
and outer field maps; `inner_seed` generates the GL20 updates.
The seed expansion does not certify every individual seed.

`Backend::Automatic` selects the AVX-512F/VL/BW/DQ, VBMI, and GFNI kernel when
available, otherwise the literal SSE2 fallback. `Portable` forces the fallback;
`Avx512` requires the fast backend. Packet encoding rejects `Avx2`.
`SPIN_ENABLE_AVX512=OFF` supports the fallback without compiling the fast kernel.
The fallback is for compatibility and checking, not comparable performance.
`PreparedEncoder` supports this family in `Full` mode only.

Owned buffers and packet scratch default to `MemoryPolicy::Automatic`: allocations
of at least 16 MiB use 2 MiB-aligned, rounded storage with best-effort Linux huge-page
advice. Smaller allocations and unsupported platforms keep 64-byte alignment.
The threshold applies separately to each allocation, before rounding.
`Normal` keeps 64-byte alignment; `PreferHugePages` requests 2 MiB alignment
at every size. Normal packet scratch retains its best-effort page advice.
Preparation never advises caller-owned buffers or runs inside encoding.
Actual huge pages are not guaranteed.

The RS family has parameter identifier 5. Identifier 4, used by the previous BCH
packet family, is retired and rejected. That implementation is no longer selectable.
The three paper-SPIN identifiers and their seeded maps are unchanged.

`PacketCode({K,seed})` is a convenience wrapper for
`Code({K,Parameters::PacketRsT64S20,seed,seed})`. Its 32-byte `SPKP` descriptor
now has version 2 and family 2; the former version 1/family 1 identified BCH packets.
This is an intentional construction change in library version 0.3, not a
seed-compatible replacement. New consumers should use the common `Code` interface.

`spin_packet_test` checks independent frozen research known answers,
portable/fast equality, natural lengths, buffer guards, and workspace ownership.
It also checks forward/transpose adjoint identities using an independent literal
forward oracle. The promoted package passes all six tests and an installed-package
consumer on Linux/GCC 15.2 and Windows/MSVC, including an AVX-512-disabled Windows
build. Linux exercises the fast packet kernel; the Windows host uses the fallback.
The AVX-512-disabled Linux packet test also passes AddressSanitizer and
UndefinedBehaviorSanitizer with leak detection enabled.
Generate/check the private RS kernels with
`python -B tools/import_rs_packet.py --check` from this directory.
Ordinary builds and installed consumers need neither Python nor research sources.

Build `spin_packet_bench` with `SPIN_BUILD_BENCHMARKS=ON`.
`spin_packet_bench K [seed] [calls] [auto|normal|huge]` measures precomputed
in-place encoding, excluding setup and allocation. Run benchmarks serially.
The frozen research kernel measured 3.259 ms at K=2^20 on Ryzen 7950X,
GCC 15.2, with four seeds, two execution orders, five warmups, and 501 timed calls.
The library imports that fused routing and outer schedule directly.

The promoted public `Code` path was compared against that frozen binary on the
same host and core, using two seeds, both execution orders, and 501 calls:

| K | Frozen reference (ms) | Library (ms) | Buffer policy |
|---|---:|---:|---|
| 2^16 | 0.21290 | 0.21215 | Normal |
| 2^18 | 0.83851 | 0.82970 | Normal |
| 2^20 | 3.29603 | 3.27232 | PreferHugePages |

Entries are medians of four process medians. Every matched output checksum
agrees; all differences are within 1.1%. These measurements check preservation
of the reference performance, not a speedup from packaging. `Automatic` selects
the listed allocation policy for these sizes on Linux.

For future changes within the RS family, retain the accepted executable and use:

```sh
bash spin/tools/compare_packet.sh /path/to/accepted/spin_packet_bench \
  /path/to/candidate/spin_packet_bench /existing/external/log/directory
```

Use matched build settings and memory policy. This serial check compares timing
and checksums at K=2^16, 2^18, 2^19, and 2^20. It is not a shared-CI timing assertion
and cannot compare checksums across different code families. Keep raw logs out
of version control.

## Identity and memory

For every `Code` family, `descriptor()` returns 40 canonical bytes. A consumer can hash them using its
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

The descriptor version and parameter identifier together fix the construction,
map definitions, setup expansion, and coordinate conventions. Backend, tile size,
and packing do not affect the descriptor. A changed realized map requires a new
identifier or version, not a silent implementation update.

`setup_bytes()` and `workspace.bytes()` report owned storage capacities,
excluding inputs, outputs, allocator overhead, and small object headers.
IMT setup prepares both directions and retains routing for bit-packed forward.
Packet setup prepares transpose only.

## Validation and remaining integration

`spin_api_test` checks dense-reference equality, the forward/transpose adjoint
identity, generic byte/64-bit equality, packed-bit equality, and wide lane equality.
It covers natural lengths, partial tiles, exact-size and range-direct paths,
invalid buffers, descriptor fields, workspace ownership, concurrent calls with
separate scratch, and allocation-free successful encoding.
`spin_known_answers` freezes descriptor-v1 forward and transpose outputs for
all three paper-SPIN parameter sets, independent of compiler and selected backend.
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
The packet section above links the RS construction's certificate and measured reference.
