# Naturally aligned transposed encoding

The `spin_half_transpose` target accepts an actual length K, not only an exponent.
The caller chooses the supplied parameter set and is responsible for its
distance/security certificate. The encoder enforces structural validity and
implementation limits; it neither chooses a different map nor rounds K up.

```cpp
using namespace bare_spin;
Spin code(Configuration::T128S19, MessageLength{81920});
code.compact();
Spin::Workspace work(code);
// 163840 input elements -> 81920 output elements, each currently 128 bits.
code.encode(input.data(), input.size(), output.data(), output.size(), work);
```

The existing `Spin(Configuration::T128S19, 20)` still means K=2^20.
The explicit `MessageLength` wrapper prevents ambiguity with that API.
Callers must rebuild against the generated headers; this is not ABI compatible.

## Parameters and lengths

This target uses BCH [256,128]. There are M=K/128 outer blocks, and the
permutation forms 256 regions of M elements. Each region contains whole inner
steps, so M must be divisible by t. The requirement is therefore K divisible
by 128t, not merely K divisible by 128.

| Supplied configuration | t | s | Rounds | K alignment |
|---|---:|---:|---:|---:|
| `T128S19` | 128 | 19 | 1 | 16384 |
| `T64S12` | 64 | 12 | 1 | 8192 |
| `T64S12R2` | 64 | 12 | 2 | 8192 |

Positive aligned K is accepted within the 32-bit routing representation: K < 2^31.
The former 2^26 policy cap is removed. Allocation can still fail when the requested
setup and buffers exceed available memory. The smallest supplied map permits K=8192;
K=4096 requires another map or a different region-boundary convention and is
not implemented here. Arbitrary user-defined (t,s) pairs are not synthesized.

The existing K16 certificates are documented recommendations. Their bounds
do not automatically apply to other sizes. The old K16-only constructor guard
has been removed from the transpose target, including its exponent API.
No proof search or certification runs during setup or encoding.

## Layout and fast paths

`Layout::Auto` is the default for checked and unchecked encoding and compaction.
It chooses Packed24 when 2K indices fit in 24 bits, otherwise Indices32.
Thus K=2^23 still permits Packed24; the next supported length uses Indices32.
An explicit Packed24 request beyond its range throws in the checked API and
in `compact`, instead of truncating. As usual, `encodeUnchecked` requires a
structurally valid object, buffers, workspace, and a supported retained layout.
After compaction, Auto uses the retained layout. Explicit requests for a
discarded layout remain errors in the checked API.

Internal tiles remain powers of two. Setup accounts for the actual number of
entries in the last tile; K itself need not be a power of two. A separate
partial-tile kernel handles the shorter final tile. The existing full-tile
kernel bodies remain unchanged, without extra per-element bounds checks,
division, or `min` operations. Path selection occurs before the hot loops.

The K16/K18 direct paths remain exact-size specializations. K18's specialized
S19 kernel is explicitly restricted to that map, so selecting a different map
at K18 cannot accidentally execute S19. With AVX-512, other lengths through
K=458752 use a runtime-length direct-routing kernel, with compile-time map
selection. This removes the tiled gather/scatter pass while the working set is
small enough. Larger lengths use the tiled path. The crossover is a performance
choice, not part of the code or certificate; `SPIN_DIRECT_MAX_K` can override it,
or disable the range extension with zero. This does not disable the existing
exact-size specializations.

The direct schedule is validated against the same permutation before
compaction. It uses 32-bit indices even if the logical layout is Packed24,
just as the existing exact-size direct paths do; compaction removes the unused
tiled schedules. There is no runtime type erasure or allocation in encoding.
Non-power-of-two sizes can share the same full-tile path when their tiles divide
evenly; the extra path is needed only for an actual partial tile.

The transposed API still consumes all input into workspace before writing
output. `encodeInplace` preserves the suffix beyond K. Setup is immutable after
compaction; each concurrent caller needs separate workspace.

## Scope and validation

The same `MessageLength{K}` interface and alignment rules now apply to the
optional bidirectional library and wide forward kernels. Both directions share
`LengthGeometry.h`; neither uses a benchmark-size or certificate-size cap.
Wide kernels support both Packed24 and Indices32, including partial tiles.
The older results below document the original transpose generalization.

For the combined natural-length validation, run `natural_lengths_check.sh ROOT UPSTREAM MODE`,
with MODE `release`, `avx2`, `masked`, or `sanitize`. The release build also
provides binaries for `natural_large_check.sh ROOT`, which checks nonzero data
just beyond the packed-index boundary at K=8,404,992. It covers standalone
transpose, bidirectional encoding, and both wide element widths. Allocation-free
geometry tests cover the former 2^26 cap and the 32-bit representation boundary;
they do not claim that allocations near the boundary fit available memory.

After all builds finish, `natural_regression_check.sh ROOT` runs the selected
power-of-two benchmarks serially and audits ISA isolation. Logs and timings
stay ignored under `measurements/natural_lengths`.

The combined change passed all eight selected suites in release, AVX2-only,
and ASan/UBSan builds on Peach. The large nonzero-data checks passed in all
three targets, including both wide widths, and the ISA-isolation audit passed.
The upper representation boundary is checked without allocating enormous buffers.

Power-of-two forward spot checks measured 0.345, 1.604, and 8.966 ms for one
128-bit stream at K=2^16, 2^18, and 2^20. The preceding measurements were 0.348,
1.601, and 9.031 ms. These are consistent with preserving the fast paths;
they are not a simultaneous old/new comparison. The new results use three
processes with 31 trials each, the same maps and compiler tuning, and a
256-row tile. Medians are taken across process medians. Wide measurements
likewise remain in the prior range, with greater small-size timing variation.
The optional `SPIN_FORWARD_FOUR` optimization is disabled by default; see
[FORWARD_PROGRESS.md](FORWARD_PROGRESS.md) for validation and measurements.
The CPU fallback is still the existing AVX2 kernel. A separate, opt-in
[generic element-type fallback](GENERIC_ELEMENTS.md) preserves the same map.

`lengths_test.cpp` checks three maps, aligned lengths, partial/full tiles,
AVX2/automatic backends, both index layouts, compaction, in-place encoding,
and equality with the independent dense oracle. Additional large tests cross
the 24-bit boundary and exercise K=2^26 against AVX2 with nonzero input.
`lengths_check.sh` builds a pre-generalization control, the new implementation,
a feature-masked build, an AVX2-only build, and a sanitizer build. It runs benchmarks serially
after builds/tests finish, with both shared benchmark locks held.

`SPIN_GENERAL_LENGTHS=OFF` exists for regression comparison with the historical
interface. Normal builds default to ON. Results and hashes belong under the
ignored `measurements/lengths` directory; generated benchmark data is not
part of the source change.

## Measured range crossover

On the Ryzen 9 7950X host with GCC 15.2 and `SPIN_TUNE=znver4`, extending direct
routing beyond the exact-size specializations improved these S19 timings.
Each entry is the median of three process medians, with 31 timed in-place
encodes per process after three warmups. Setup is excluded; one logical element
is 128 bits. All before/after output hashes match.

| K | Previous dispatch (ms) | Range direct (ms) |
|---:|---:|---:|
| 16384 | 0.102 | 0.077 |
| 81920 | 0.503 | 0.436 |
| 163840 | 1.031 | 0.884 |
| 245760 | 1.555 | 1.360 |
| 278528 | 1.903 | 1.543 |
| 327680 | 2.402 | 1.864 |
| 393216 | 3.069 | 2.526 |
| 458752 | 3.759 | 3.456 |
| 524288 | 4.274 | 4.681 |

The last row is an experiment beyond the selected crossover: direct routing
loses once the working set grows, so the default keeps tiled routing there.
The crossover was calibrated for S19 on this host, not claimed optimal for
every CPU or map. `lengths_range_check.sh` reproduces this comparison;
`lengths_summary.py` checks hashes and summarizes the measurements.

The final default dispatch was then checked against the pre-generalization
build at the existing tuned points (three process medians, 51 encodes each):

| K / configuration | Previous (ms) | Generalized (ms) |
|---|---:|---:|
| 2^16 / T64S12R2 | 0.343812 | 0.343241 |
| 2^18 / T128S19 | 1.424228 | 1.426953 |
| 2^20 / T128S19 | 9.282557 | 9.188662 |

These differences are within roughly 1%; this run does not indicate a
power-of-two regression. All paired output hashes agree. The new source passed
the ordinary, AVX-512-masked, AVX2-only, and ASan/UBSan test suites. The masked
suite intentionally skips the AVX-512 basis test. ISA inspection confines EVEX
instructions to the AVX-512 objects. The nonzero-input large tests passed at
K=8404992 (past the 24-bit boundary, with a partial tile) and K=67108864.

The subsequent [generic XOR-element fallback](GENERIC_ELEMENTS.md) supports
other element types without extending the forward API.
