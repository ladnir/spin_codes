# SPIN scheduling exploration, 2026-09-19

## Outcome

Keep the existing half-rate IMT implementation. This bounded experiment found
no improvement from longer prefetching, explicit scatter scheduling, early
index decoding, 32-bit routing indices, or the adjacent tile sizes.
The fresh reference reproduced the paper's approximately 10.11 ms latency.
No production source, certificate, or submitted paper was changed.

The new prime-field BAA result motivated the experiment: could similar memory
scheduling improve SPIN without changing its encoder? Unlike that gather-heavy
implementation, SPIN already uses packed tables, tiled routing, and prefetching.
The negative result applies to the tested kernels, not every possible schedule.

## Workload and method

The encoder is the certified BCH [256,128] / IMT instance with t=128, s=19,
weight-five feedback `weight5_seed0`, and K=2^20. It evaluates the transpose
on 128 parallel binary instances represented by 16-byte blocks.
All candidates preserve the exact route, maps, state recurrence, and output order.

Peach supplies one AMD Ryzen 9 7950X core, pinned to CPU 15. GCC 15.2 uses
Release `-O3 -DNDEBUG -march=znver4 -mavx2 -mpclmul -mvpclmulqdq`.
The host had the userspace governor and boost disabled; no frequency settings
were changed. The environment receipt contains an instantaneous frequency
reading, not a guarantee of constant frequency throughout measurement.

Each process initializes its input once, performs three warmups, then encodes
in place without resetting the input. Setup, allocation, workspace preparation,
and output hashing are outside the interval. Large-page advice is enabled for
owned workspace, as in the reference. There is no explicit cache flush.
The page-advice helper is best-effort; this campaign does not separately
measure successful collapse or workspace-preparation latency.

Benchmarks acquire both `/tmp/prindal-addition-encoder-benchmark.lock` and
`/tmp/bare-spin-benchmark.lock`, in that order. They also reject other active
benchmark executables. All measurements and validation runs were serial.

The first screens use one process of 31 timed calls per cell. Final comparisons
use three processes of 101 calls per cell, reversing candidate order in the
middle repetition. Two route seeds, 1 and 17, test different sampled routes;
the coefficient seed stays 2. Each result below is a median of process medians.

## Confirmation

| Configuration | Route seed 1 (ms) | Route seed 17 (ms) |
|---|---:|---:|
| Reference: packed24, 2048 rows, prefetch 32 | **10.118467** | **10.220417** |
| Prefetch 16, otherwise unchanged | 10.331745 | 10.340992 |
| Tile 1024 rows, otherwise unchanged | 10.305366 | 10.297892 |

Prefetch 16 regresses by 2.11% and 1.18%. The smaller tile regresses by 1.85%
and 0.76%. Every paired candidate process is slower than its control.
The reference process-median ranges are 10.110392--10.181395 ms for seed 1
and 10.188408--10.254410 ms for seed 17. Two seeds do not establish a
distribution-wide performance guarantee.

Retained setup is 12,713,992 bytes with packed24, versus 16,908,288 bytes
with 32-bit indices. Workspace is 40 MiB at 2048 rows, 36 MiB at 1024 rows,
and 48 MiB at 4096 rows. No candidate adds a retained full-size buffer.
Early decoding uses 512 bytes of bounded stack scratch per call.

## Rejected screens

These are single-process selection results, not separately confirmed rankings.
The control before and after the first screen takes 10.102 and 10.116 ms.

| Scheduling change | Complete latency (ms) |
|---|---:|
| No tile-scatter prefetch | 10.591 |
| Tile-scatter read-prefetch 128 | 10.602 |
| Tile-scatter read-prefetch 256 | 10.366 |
| Four-entry scatter scheduling, prefetch 128 | 10.989 |
| Eight-entry scatter scheduling, prefetch 128 | 11.323 |
| Decode an epoch's 128 destinations before IMT emission | 11.178 |
| Early decoding plus future-epoch write-prefetch | 13.491 |
| Future-epoch write-prefetch without early decoding | 12.975 |

The refinement screens try read-prefetch leads 8, 16, and 64, plus four-entry
scatter scheduling at lead 32. They take 10.328, 10.315, 10.516, and 10.791 ms,
against a 10.113 ms control. Prefetch leads count future block accesses,
not bytes. Epoch prefetch leads are 128 or 256 block positions in reverse order.

| Tile rows | Packed24 (ms) | Indices32 (ms) |
|---|---:|---:|
| 1024 | 10.307 | 11.312 |
| 2048 | 10.272 | 11.210 |
| 4096 | 10.860 | 11.772 |

Packed indices remain preferable in this screen. Removing their decoding
does not compensate for the tested implementation's larger index streams.
This comparison does not isolate bandwidth from instruction scheduling.

## Where the time goes

A separate instrumented build times IMT plus initial bucket routing, tile
scatter, and the existing paired BCH transpose. It averages 34 calls,
including warmups. Instrumentation changes generated code and the execution
schedule: its approximately 10.8 ms total must not replace the uninstrumented
approximately 10.1 ms measurement.

| Tile rows | IMT + bucket routing (ms) | Tile scatter (ms) | BCH (ms) |
|---|---:|---:|---:|
| 1024 | 3.438 | 3.259 | 4.102 |
| 2048 | 3.202 | 3.478 | 4.158 |
| 4096 | 3.182 | 4.045 | 4.245 |

The smaller tile improves the local scatter but makes initial routing more
expensive. At 2048 rows, routing feeds four bucket streams and scatters within
an 8 MiB tile. Halving the tile doubles the number of bucket streams; doubling
it makes the local scatter cover 16 MiB. This explains the competing locality
objectives, although these timings do not identify particular cache or TLB stalls.

The important new lead is the BCH stage: approximately 38% of this instrumented
total. The inner is not the only substantial cost. Before changing IMT or the
permutation, investigate register pressure and batching in the existing BCH
transpose. Any improvement there could preserve the exact code and proof too.
No claim is made here that BCH has an available speedup of a particular size.

## Correctness and scope

All 13 non-profile encoder targets pass the existing dense-oracle suite at
K=2^16, 2^18, and 2^20. It covers packed and 32-bit indices, in-place suffix
preservation, compaction, alternate tiles, zero input, and boundary impulses
with another setup. Complete-output hashes also agree across every timed
candidate with the same route seed and iteration count.

ASan and UBSan pass the prefetch-16 and early-decoding/prefetch encoders against
that suite. A separate helper test covers lengths 0 through 513, both index
formats, nonzero base offsets, tail loops, and output canaries under sanitizers.
Future-address prefetches are bounded; no out-of-range metadata loads are used.
The existing adjoint and circuit proofs were not modified or rerun here.

There is no winning half-rate candidate to transfer to quarter rate, so no
quarter-rate benchmark was run. New multilevel staging was also deferred:
the earlier routing experiment rejected several such alternatives, and this
screen did not supply evidence to revisit them yet.

## Files and replay

- `generate.py` creates isolated copies from the pinned weight-five source.
  The `control` target compiles the original source, not a rewritten equivalent.
- `Schedule.h` contains fixed-size, compile-time scheduling helpers.
- `CMakeLists.txt` builds control, candidates, and tests.
- `run.sh`, `refine.sh`, and `confirm.sh` record serial timing campaigns.
- `sanitize.sh` runs selected encoders and helper tails under ASan/UBSan.
- `summary.py` checks final C++ source hashes, sample medians, and output hashes.
  It summarizes results but is not a proof verifier.
- `measurements/` holds local-only samples, validation logs, environment
  information, and source/binary hash receipts. It is ignored by Git.

The isolated remote root was `/tmp/spin-scheduling-20260919`. It contains the
minimal source tree copied from the earlier weight-five experiment, not a new
checkout of the paper. From a repository root with the same dependencies:

```sh
python3 workstreams/inner_design/scheduling_20260919/generate.py
cmake -S workstreams/inner_design/scheduling_20260919 -B build-scheduling -DCMAKE_BUILD_TYPE=Release
cmake --build build-scheduling -j3
bash workstreams/inner_design/scheduling_20260919/run.sh "$PWD"
bash workstreams/inner_design/scheduling_20260919/refine.sh "$PWD"
bash workstreams/inner_design/scheduling_20260919/confirm.sh "$PWD"
bash workstreams/inner_design/scheduling_20260919/sanitize.sh "$PWD"
python3 workstreams/inner_design/scheduling_20260919/summary.py
```

Do not run these scripts concurrently. Preliminary screens used evolving
instrumentation; the final confirmation has separate source and binary receipts.
The later addition of sanitizer test targets does not alter timed source files.
No results have been promoted, committed, or pushed.

Next: inspect the paired BCH-256 transpose's generated assembly and measure
an isolated, cache-policy-matched outer baseline before selecting a new kernel.
