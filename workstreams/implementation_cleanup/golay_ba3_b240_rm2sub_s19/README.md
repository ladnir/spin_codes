# Golay--BA-3/RM2Sub-S19 for libOTe

This directory implements the transposed encoder for the finite `B=240`,
`k=2^20` SPIN design point. It is a separate construction from the frozen
Structured SPIN `B=256` encoder.

The implementation provides two layouts for the BA-3 setup:

- `GolayBa3OuterMode::Reused` uses one pair of BA permutations for all rows.
- `GolayBa3OuterMode::Independent` uses a distinct pair for each row.

Both layouts evaluate the same sequence of linear components. Only the setup
distribution and permutation-table traffic differ.

## Parameters and public map

The parent code has these fixed parameters:

| parameter | value |
|---|---:|
| outer block length `B` | 240 |
| outer block dimension | 120 |
| outer rows `L` | 8,832 |
| parent code coordinates | 2,119,680 |
| parent message coordinates | 1,059,840 |
| returned message coordinates | 1,048,576 |
| fixed shortened coordinates | 11,264 |
| RM2Sub step width | 128 |
| RM2Sub epochs | 16,560 |

`dualEncodeTo` and `dualEncodeUnchecked` compute a transposed map

```text
2,119,680 input blocks -> 1,048,576 output blocks.
```

Each `osuCrypto::block` holds 128 parallel binary instances. The public map
returns the first `2^20` parent-message coordinates. Equivalently, the forward
encoder fixes the last 11,264 parent-message coordinates to zero.

An ordinary forward encoder is not included. The intended libOTe interface is
the transposed block encoder.

## API

Include the public header:

```cpp
#include <libOTe/Tools/RiffleCode/GolayBa3Rm2SubB240K20.h>
```

Initialize one immutable setup object and allocate one workspace per active
call:

```cpp
using osuCrypto::GolayBa3OuterMode;
using osuCrypto::GolayBa3Rm2SubB240K20;
using osuCrypto::block;

GolayBa3Rm2SubB240K20 code;
code.init(
    0x524f5554452d3234ULL, // factored-route seed
    0x42412d332d423234ULL, // two BA permutation layers
    0x524d325355423139ULL, // RM2Sub multiplier schedule
    GolayBa3OuterMode::Reused);

GolayBa3Rm2SubB240K20::Workspace workspace;
std::vector<block> input(GolayBa3Rm2SubB240K20::codeBlocks);
std::vector<block> output(GolayBa3Rm2SubB240K20::messageBlocks);

code.dualEncodeTo(
    input.data(), input.size(),
    output.data(), output.size(),
    workspace);
```

`dualEncodeTo` checks initialization and both buffer lengths.
`dualEncodeUnchecked` removes those checks for a validated hot call. The input
and output buffers must not overlap.

After `init` returns, concurrent calls may share the setup object if no thread
calls `init` again. Each concurrent call requires a separate `Workspace` and
separate input and output buffers. Construction and destruction of a workspace
allocate and release its scratch vectors, so callers should reuse workspaces.

The three seeds have distinct roles:

- `routeSeed` derives the per-row coordinate permutations and per-region row
  permutations.
- `baSeed` derives both BA permutation layers. Independent mode derives two
  layers for every outer row.
- `innerCoefficientSeed` derives the nonzero RM2Sub-S19 field multipliers.

`init` is setup work and is not part of the measured online path.

## Construction evaluated by the hot path

For each outer row, the forward map applies ten extended binary Golay
`[24,12,8]` encoders. It then applies this BA-3 sequence:

```text
Golay^10 -> permutation 1 -> accumulator 1
         -> permutation 2 -> accumulator 2.
```

The global route applies one coordinate permutation per row, followed by one
row permutation in each of the 240 regions. RM2Sub-S19 then processes 16,560
fixed-width epochs.

The transposed encoder evaluates these components in reverse order. It reuses
the existing fixed-width and unrolled RM2Sub-S19 kernel. The outer transpose
uses two fixed 240-block stack arrays, explicit suffix accumulators, byte-packed
BA indices, and an unrolled Golay transpose.

The route is precomposed during setup and split into eight cache buckets. Each
bucket holds 1,104 complete outer rows. RM2Sub emission writes sequentially to
the bucket streams. A bucket-local scatter restores outer-coordinate order in
a 4.04 MiB tile. The bucket layout changes storage order only; it preserves
the sampled permutation exactly.

No hot-path virtual dispatch, `std::function`, type erasure, or per-call heap
allocation is used.

## Memory

The reported sizes include vector capacities after initialization:

| object | reused | independent |
|---|---:|---:|
| setup schedules | 14,308,328 bytes | 18,547,208 bytes |
| workspace scratch | 38,154,240 bytes | 38,154,240 bytes |

Independent mode adds 4,238,880 bytes of byte-packed BA indices. The public
object and workspace handles themselves use fixed inline storage; their owned
vectors hold the large setup and scratch arrays.

## Build and test

The standalone build imports the cleaned Structured SPIN implementation for
the RM2Sub kernel and the performance baseline:

```bash
cmake -S . -B build -G Ninja \
  -DCMAKE_BUILD_TYPE=Release \
  -DSTRUCTURED_SPIN_IMPLEMENTATION_ROOT=/path/to/structured/implementation
cmake --build build -j
ctest --test-dir build --output-on-failure
```

The tests validate all of the following:

- the precomposed route is a permutation of 2,119,680 coordinates;
- each BA table is a permutation of 240 coordinates;
- the Golay generator has spectrum
  `1 + 759 z^8 + 2576 z^12 + 759 z^16 + z^24`;
- the optimized RM2Sub transpose equals an independent scalar recurrence;
- the optimized BA and Golay transpose equals a scalar matrix calculation;
- the shortened output matches its deterministic checksum.

| mode | expected checksum |
|---|---:|
| reused | `0xc4108a1282911d32` |
| independent | `0x2bd00619556f9a86` |

The public-header smoke test compiles without exposing intrinsic-heavy internal
headers to downstream code.

## Upstream libOTe integration

First apply the cleaned Structured SPIN overlay because this implementation
shares its RM2Sub-S19 kernel. Then:

1. Copy this directory's `libOTe/` tree into the libOTe checkout.
2. Apply `integration/libote-after-structured-spin.patch`.
3. Build the normal libOTe target.

The patch adds `GolayBa3Rm2SubB240K20.cpp` to `libOTe/CMakeLists.txt`. It keeps
AVX2, PCLMUL, and VPCLMULQDQ flags local to this source file. These flags do not
become usage requirements on the public `oc::libOTe` CMake target.

The tested implementation requires a processor with AVX2 and VPCLMULQDQ.
Runtime ISA dispatch is not implemented.

## Performance

The final comparison ran on Peach with an AMD Ryzen 9 7950X, GCC 15.2.0,
`-march=znver4`, and CPU 15. Each result is the median of 101 steady-state
calls. Before each process, the benchmark guard scanned `/proc/[0-9]*/exe`
and found no other benchmark executable.

| implementation | code blocks | median | change from B=256 baseline |
|---|---:|---:|---:|
| Structured SPIN B=256 baseline | 2,097,152 | 10.696591 ms | -- |
| Golay--BA-3 B=240, reused | 2,119,680 | 10.263112 ms | -4.052% |
| Golay--BA-3 B=240, independent | 2,119,680 | 10.407942 ms | -2.699% |

Independent mode is 1.411% slower than reused mode. The comparison measures
only transposed online encoding. It excludes setup, allocation, destruction,
conditioning, and schedule authentication.

The complete machine-readable record is in `PERFORMANCE_RECEIPT.json`.

## Setup and proof boundary

`init` deterministically expands three 64-bit seeds. This behavior supplies a
reproducible online implementation, but it does not instantiate the probability
space used by the finite analysis.

The finite analysis samples uniform permutations and independent nonzero field
multipliers. Its conditioned ensembles also require the BA code to satisfy
the event `G_240`: every nonzero BA word has weight from 25 through 215. The
analysis does not provide an efficient algorithm for deciding `G_240`.

Consequently, this implementation makes no conditioned distance claim. A
conditioned sampler or an authenticated fixed BA schedule could reuse the same
online representation and hot path. Such a setup method must receive separate
correctness and provenance evidence.

The original finite-design note has SHA-256
`5861c4878986a483efc7ae0d4287daaf6efe3e8a0fffe6ebe7b0e7ce7f66a876`.
`SOURCE_MANIFEST.json` records this dependency and the exact RM2Sub inputs.

## Directory layout

```text
libOTe/Tools/RiffleCode/
  GolayBa3Rm2SubB240K20.h          public API
  GolayBa3Rm2SubB240K20.cpp        setup and optimized transpose
  GolayBa3Rm2Sub_TestAccess.hpp    test-only access
tests/
  golay_ba3_correctness.cpp        scalar oracle and checksums
  public_header_smoke.cpp          public include boundary
benchmarks/
  golay_ba3_comparison.cpp         guarded serialized comparison
integration/
  libote-after-structured-spin.patch
PERFORMANCE_RECEIPT.json
SOURCE_MANIFEST.json
```

## Recommended next step

Authenticate one reusable BA schedule, or implement an efficient procedure for
the `G_240` conditioning event. Then add a schedule-import API that preserves
the current packed representation and measured online path.
