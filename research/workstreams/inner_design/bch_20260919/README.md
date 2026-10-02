# Faster BCH-256 transpose for SPIN

2026-09-19. The selected isolated implementation, `packedshare4`, reduces
half-rate SPIN's encoding latency by about 10% at K=2^20 and 22% at the
two smaller measured sizes. It preserves the exact BCH generator, IMT maps,
sampled permutation, and externally visible output order.

The submitted paper, its performance tables, and the pinned implementation
sources are unchanged. This is a new implementation candidate, not a changed
code construction or a new distance certificate.

## Confirmed complete-encoder timings

Times are milliseconds, measured on Peach's Ryzen 9 7950X, one thread pinned
to CPU 15. Each cell is the median of three process medians, each over 101
in-place encodes after three warmups. Inputs are initialized once and never
reset between calls. Candidate order reverses in the middle repetition.

| Message length | Existing implementation | Selected BCH implementation | Latency reduction |
|---|---:|---:|---:|
| K=2^16 | 0.523477 | **0.404305** | **22.77%** |
| K=2^18 | 2.155534 | **1.680456** | **22.04%** |
| K=2^20 | 10.141350 | **9.089597** | **10.37%** |

The table uses route seed 1 and coefficient seed 2. At K=2^20, a second
route seed gives 10.201943 versus 9.129792 ms, a 10.51% latency reduction.
Every paired selected process beats its control. The selected process-median
ranges at K=2^20 are 9.053008--9.149919 ms for seed 1 and
9.122809--9.287176 ms for seed 17. These measurements do not establish a
performance guarantee across machines or all sampled permutations.

The workload is the BCH [256,128] / IMT (t=128,s=19) transpose with feedback
`weight5_seed0`. Each 16-byte block carries 128 parallel binary instances.
No outer, inner, or distance parameter was retuned.

## What changed

The existing BCH implementation shares XOR subexpressions across all 128
outputs, computes all intermediate terms first, and evaluates two code blocks
with 256-bit vectors. Its compiler-generated stack frame is 24,168 bytes.
Assembly shows extensive stores and reloads of vector temporaries.

The selected implementation combines three changes:

1. Compute each output after recursively computing the terms it needs.
   Shared terms remain available for later outputs, but input loads and
   intermediate computations no longer all precede the first output.
2. Process four independent code blocks with 512-bit vectors. The temporary
   tile interleaves their coordinates, so each vector loads contiguously.
   The existing scatter fills this layout directly; there is no extra pass.
3. Share an XOR subexpression only when the synthesis heuristic finds it in
   at least four outputs. This adds arithmetic but reduces retained intermediates.

The selected circuit has 3,758 algebraic XORs per four-block evaluation,
compared with 3,029 per two-block evaluation in the original circuit.
These are symbolic circuit counts, not emitted instruction counts: the compiler
can combine XORs with ternary-logic instructions. The selected stack frame is
19,144 bytes. Four-row batching with the original sharing rule uses 31,688 bytes.

Fewer XORs alone was not the right optimization objective. The effective cost
also includes SIMD width, input packing, live temporaries, and spill traffic.
These experiments show that the combined change helps; they do not attribute
an exact fraction of the gain to individual microarchitectural effects.

## Why the map and proof are unchanged

The generator evaluates every synthesized signal as a 256-bit formal linear
form. It checks all 128 outputs against the committed `BchRows` matrix.
Changing XOR sharing therefore does not change the linear map.

The new tile layout only relabels owned workspace. For an old local tile
offset `local`, the stored position becomes

```text
(local & ~1023) | ((local & 255) << 2) | ((local >> 8) & 3)
```

The high bits select a group of four BCH rows. Within that group, the low two
bits select a SIMD lane and the next eight bits select the coordinate.
The inverse is

```text
(packed & ~1023) | ((packed & 3) << 8) | ((packed >> 2) & 255)
```

Setup stores the transformed offsets in the existing routing tables. No
additional per-coordinate arithmetic is added to the scatter loop. The
four-row BCH kernel reads this layout and stores results in ordinary row order.
The logical route remains identical; it is not a different SPIN permutation.

Retained setup remains 12,713,992 bytes at K=2^20 with packed24 routing.
Workspace remains 40 MiB with the 2048-row tile. No heap allocation or dynamic
dispatch is added to the online path. The polynomial-basis certificate and the
IMT proof need no change when this exact implementation equivalence is retained.

## Other candidates

The short screens use one process of 31 calls, so their small differences
are selection evidence rather than confirmed rankings.

| Candidate | Complete time at K=2^20 (ms) |
|---|---:|
| Original control, first screen | 10.229; ending control 10.123 |
| Demand-driven two-row schedule, unchanged XOR circuit | 10.000 |
| Split original circuit into 32-output routines | 11.155 |
| Resynthesize each 32-output routine | 10.382 |
| Resynthesize each 16-output routine | 10.979 |
| One-row, 128-bit SIMD schedule | 13.182 |
| Four-row SIMD, ordinary tile layout | 9.686 |

Reduced two-row sharing thresholds of 3 and 4 take approximately 10.00 ms.
Threshold 8 takes 10.54 ms. Neither aggressive sharing nor aggressive
recomputation wins universally.

The final K=2^20 comparison, using the same three-process protocol as the
headline, separates the packed four-row finalists:

| Packed four-row candidate | Route seed 1 (ms) | Route seed 17 (ms) |
|---|---:|---:|
| Original sharing threshold 2 | 9.293477 | 9.303576 |
| Sharing threshold 3 | 9.168223 | 9.193551 |
| Selected sharing threshold 4 | **9.089597** | **9.129792** |

## BCH-only measurements

A separate synthetic benchmark processes the equivalent of 8192 rows per call.
The hot mode repeatedly processes one SIMD batch. The streaming mode processes
32 MiB of input and writes 16 MiB of output. Neither mode includes SPIN routing;
both use resident pages without explicit cache flushing.

| Synthetic mode | Original (ms) | Selected (ms) |
|---|---:|---:|
| Repeated hot batch | 3.930055 | 2.483235 |
| Streaming rows | 4.769343 | 3.931539 |

These are medians of three process medians, each over 31 calls. They confirm
a kernel-level gain in both regimes. They are not the BCH contribution inside
SPIN: cache state, interleaving with routing, and the input layout differ.
The full-SPIN table includes every routing and layout cost and is the acceptance
criterion. The larger percentage gain at smaller K is consistent with the hot
kernel result, but this campaign does not isolate its cause with hardware counters.

## Correctness, sanitizers, and integration requirements

All 20 release tests pass. They include full SPIN comparisons with the dense
oracle at K=2^16,2^18,2^20, both index formats, in-place output and suffix
preservation, compaction, alternate tiles, and another sampled setup.
The direct BCH tests exhaust every input basis coordinate separately in each
SIMD lane, then compare random inputs with a dense scalar matrix evaluation.
Output canaries are checked too.

ASan and UBSan pass three selected tests: the full encoder, the BCH basis/random
test, and the minimum-tile/compaction/in-place test. The latter checks tiles of
4, 8, 16, and 512 rows, in both index formats, and rejects a two-row tile.
Full-output hashes agree for all timed implementations with equal message size,
seed, and iteration count.

The selected kernel requires AVX-512 on x86. These measurements use GCC 15.2
with `-O3 -DNDEBUG -march=znver4 -mavx2 -mpclmul -mvpclmulqdq`.
The host had the userspace governor and boost disabled; no machine-wide settings
were changed. Large-page advice remains enabled for owned workspace, as in the
baseline. Setup, allocation, workspace preparation, and output hashing are
excluded from full-SPIN timings. Page-advice success was not separately recorded.

For integration, preserve a CPU-appropriate fallback and the existing public
API. In particular, these isolated four-row builds require at least four tile
rows; an integrated dispatcher should retain the old path for smaller tiles.
Low-level generated kernels require nonoverlapping input and output. The full
encoder continues to support in-place operation through its separate workspace.
Quarter-rate BCH-128 and the ordinary encoder have not been changed or benchmarked
in this campaign. Do not transfer the measured gain to them without testing.

## Reproduction and records

`generate.py` parses the committed BCH circuit and uses the repository's existing
Paar synthesis implementation. It checks exact linear forms and the workspace
relabeling before emitting isolated source files under `generated/`.
`CMakeLists.txt` links the original BCH source for the control and keeps the
remaining SPIN kernel fixed except for the four-row call and internal tile layout.

From a repository root containing the existing source dependencies:

```sh
python3 workstreams/inner_design/bch_20260919/generate.py
# Reuse the scheduling campaign's harness, including its two benchmark locks.
python3 workstreams/inner_design/scheduling_20260919/generate.py
cmake -S workstreams/inner_design/bch_20260919 -B build-bch -DCMAKE_BUILD_TYPE=Release
cmake --build build-bch -j3
bash workstreams/inner_design/bch_20260919/run.sh "$PWD"
bash workstreams/inner_design/bch_20260919/refine.sh "$PWD"
bash workstreams/inner_design/bch_20260919/confirm.sh "$PWD"
bash workstreams/inner_design/bch_20260919/sanitize.sh "$PWD"
python3 workstreams/inner_design/bch_20260919/summary.py
```

Run these commands sequentially. The timing executables acquire the prime-field
BAA lock and SPIN lock in that order and reject other active benchmark executables.
No benchmark ran concurrently with another benchmark in this campaign.

The isolated remote root is `/tmp/spin-scheduling-20260919`; this experiment
has its own `build-bch` and `build-bch-sanitize` directories. Local-only
`measurements/` contains raw samples, correctness and sanitizer logs, environment
information, stack-usage records, and source/binary hashes. It is ignored by Git.
`summary.py` checks 42 measured C++ source/header hashes, sample medians, and
full-output hash agreement. Earlier screens used evolving generators; the final
confirmation has separate source and binary receipts.

Next: integrate `packedshare4` as an explicit AVX-512 implementation option,
preserving the existing fallback, setup identities, and submitted-paper receipts.
