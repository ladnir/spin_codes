# Using the optimized forward encoders

Use this guide for the selected 128-, 256-, and 512-bit forward paths.
Forward encoding maps K elements to 2K elements. These optional targets are
separate from the default, transpose-only library.

## Build

On Linux x86-64, configure a fresh build directory from the repository root:

```sh
cmake -S workstreams/spin_optimized -B out/spin-forward \
  -C workstreams/spin_optimized/cmake/ForwardRecommended.cmake \
  -DSPIN_BIDIRECTIONAL_SOURCE=/path/to/hypercat/hypercat/native/spin
cmake --build out/spin-forward --target spin_wide -j2
```

The initial-cache file selects four-row BCH, direct routing through K=2^18,
and on-demand BCH evaluation for both single-stream and wide paths.
It does not change the default transpose build or override existing cache values.
Use a fresh directory when switching from an experimental configuration.
`SPIN_TUNE=znver4` reproduces the measured compiler tuning; the recipe defaults
to generic tuning. GCC 15.2 was tested. MSVC and ARM are unsupported here.

Link `spin_half_bidirectional` for 128-bit forward encoding or `spin_wide`
for wide encoding. The latter also links the former. Do not additionally link
`spin_half_transpose`: the libraries define alternative versions of `Spin`.

The upstream source is copied into the build tree, never edited.
Generated import manifests record source hashes and retain the imported license.
Both forward and transpose accept `MessageLength{K}` for an actual length.
The existing integer argument remains an exponent for compatibility.
K must be a positive multiple of 128*t: 16,384 for T128S19 and 8,192 for
T64S12/T64S12R2. No padding, map selection, or certificate search is implicit.
The old 2^20 forward and 2^26 transpose caps are removed. Both use the shared
32-bit routing representation (K < 2^31), with checked arithmetic and ordinary
allocation failures when resources are insufficient.

## Encode native wide input

```cpp
#include "Wide.h"
#include <vector>

using namespace bare_spin;
constexpr unsigned lanes = 2; // 256-bit elements; 4 selects 512-bit elements
// Check wideAvailable(lanes) before selecting a width for the application.
Spin code(Configuration::T128S19, MessageLength{81920}, 1, 2, 256);
code.compact(); // Auto selects Packed24 or Indices32 and retains that layout.
WideWorkspace workspace(code, lanes);

std::vector<block> input(lanes * code.messageBlocks());
std::vector<block> output(lanes * code.codeBlocks());
// Populate input[i*lanes + lane] with lane's 128-bit word at coordinate i.
workspace.forward(input, output);
// output[j*lanes + lane] has the corresponding encoded coordinate.
```

Each 128-bit lane contains 128 independent binary messages. Width lifting
applies the same binary generator independently to every bit; it introduces
no extension-field multiplication or new code setup.

Input and output must be disjoint, with exactly K*lanes and 2K*lanes blocks.
Ordinary 16-byte block alignment suffices. Encoding reuses the workspace and
allocates no per-call scratch. `workspace.bytes()` reports scratch capacity,
excluding the setup, input, and output.

The `Spin` object must outlive its workspace and remain unmoved and unmodified.
Use a separate workspace per concurrent caller. Constructors reject unsupported
widths or CPU features before entering ISA-specific code. The 256-bit path uses
the existing AVX2 baseline; 512-bit elements require AVX-512F/VL/BW/DQ.
Wide workspaces support `T128S19` and `T64S12R2` with either routing layout.
Auto uses Packed24 through K=2^23, then Indices32. Explicit Packed24 requests
beyond its range fail instead of truncating indices. After compaction, Auto
uses the retained layout. Partial tiles have separate kernels; complete-tile
hot loops keep their existing structure, including power-of-two fast paths.

For one 128-bit element stream, use `Spin::Workspace` and
`code.forward(input.data(), K, output.data(), 2*K, workspace)` instead.
No packing or column assembly is part of either encoder API.

## Choose width and tile

Keep separate 128-bit streams separate unless the application benefits from
interleaving them elsewhere. At larger measured sizes, packing only for encoding
costs more than it saves. Native wide input avoids this conversion.

For native wide input, 256 bits is a useful memory-conscious starting point.
The 512-bit workspace is twice as large and does not win at every size.
Width selection therefore remains explicit rather than hidden in dispatch.

| Tested K | Map used in timings | Starting tile, outer rows | Width guidance |
|---|---|---:|---|
| 2^16 | T64S12R2 | 256 | 512-bit won; 256-bit uses less scratch |
| 2^17 | T128S19 | 256 | Prefer 256-bit; a 512-row tile was about 2% faster for native 256-bit input |
| 2^18 | T128S19 | 256 | 256/512-bit totals were close |
| 2^19 | T128S19 | 256 | 256/512-bit totals were close; 512-row tiles were slower |
| 2^20 | T128S19 | 512 wide, 256 single-stream | 512-bit was slightly faster, at twice the wide scratch |

These are performance choices for the measured host, not certificate selection
or universal optima. The caller chooses the map and checks its certificate for
the desired K, distance, and margin. No new certificates are claimed for the
intermediate-size performance samples.

For native 256-bit input, total times for the tested 16-stream batch were
7.39, 16.60, 35.66, 74.06, and 156.37 ms at these starting tiles.
This gives a reasonably smooth measured progression from K=2^16 through 2^20.
It does not establish performance over the full proposed range through 2^26.
See [WIDE_ITERATION.md](WIDE_ITERATION.md) for scope, packing costs, and raw-summary commands.

## Validate and measure

```sh
cmake --build out/spin-forward --target spin_wide_test spin_wide_lengths_test \
  spin_forward_lengths_test spin_bidirectional_test spin_k16_bidirectional_test -j2
ctest --test-dir out/spin-forward \
  -R '^(wide_.*|forward_natural_lengths|bidirectional|k16_bidirectional)$' --output-on-failure -j1
```

The tests check natural lengths, partial/full tiles, generator validation, lane-by-lane equality, compaction,
unaligned buffers, malformed lengths, and full or partial overlap rejection.
The generator also checks every BCH output against the binary matrix on each
on-demand build. Its checks remain active under optimized Python.

To build the comparison harness, configure with `SPIN_BUILD_BENCHMARK=ON` and
build `spin_wide_benchmark`. For example:

```sh
out/spin-forward/spin_wide_benchmark 20 7 12819 512 2
# exponent, trials, configuration, tile rows, lanes (1, 2, 4; 0 tests all)
```

Every row measures 16 independent 128-bit element streams, not one encoder call.
`encode_ms` excludes packing and column assembly; `total_ms` includes them.
The harness requires AVX-512 for column assembly even when testing 256-bit encoding.
It pins CPU 15 and takes shared locks. Never run benchmarks concurrently or build
while a benchmark is timing.

Performance history is preserved in [FORWARD_SCHEDULE.md](FORWARD_SCHEDULE.md)
for single-stream encoding and [WIDE_ITERATION.md](WIDE_ITERATION.md) for wide encoding.
Raw measurements and build outputs remain ignored by Git.
