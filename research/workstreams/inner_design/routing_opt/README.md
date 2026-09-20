# Exact-map routing optimization

The experiments below are now followed by an opt-in integration. See
[deployment validation](deployment/README.md) for the supported build option,
size-dependent dispatch, fallback tests, and final integrated measurements.

This workstream optimizes the quarter-rate transposed encoder without changing
its inner map or permutation. The control is `asymmetric_greedy3_2_sparse`:
BCH [128,32,32], t=128, s=19, K=2^20, and 128-way bitslicing.

The best confirmed change combines optional Linux large-page backing for the two
owned workspace allocations with write-prefetch in the tile scatter.
It reduces encoding time by about 5.8% on Peach.
The tested alternative routing algorithms were slower. During this isolated
experiment, supported source files and previous certificate artifacts were unchanged.

## Confirmed measurements

Each entry reports the median of three process medians, each containing 101
in-place calls. Setup, workspace preparation, and input initialization are excluded.
The same input/output buffer is reused, with three warmups per process.
Controls and candidates are interleaved. No benchmarks run concurrently.

| Implementation | Time (ms) | Reduction | Workspace |
|---|---:|---:|---:|
| Control, 4096-row tile | 15.9904 | — | 72 MiB |
| Write-prefetch, same tile | 15.8558 | 0.84% | 72 MiB |
| Large pages, 4096-row tile | 15.4795 | 3.19% | 72 MiB |
| Large pages, 2048-row tile | 15.5002 | 3.07% | 68 MiB |

Retained setup remains 25,427,976 bytes. A row occupies 128 blocks of 16 bytes,
so these tiles occupy 8 MiB and 4 MiB. The small difference between the two
large-page configurations does not establish a robust speed preference.

The host is an AMD Ryzen 9 7950X, pinned to CPU 15, running Linux 6.17.
The build uses GCC 15.2, Release, and the existing znver4/AVX2 configuration.
These results are not a portable promise of the same percentage improvement.

A second three-run comparison tested the combination, using 4096-row tiles:

| Implementation | Process medians (ms) | Median (ms) |
|---|---|---:|
| Control | 16.7164, 15.9720, 15.9249 | 15.9720 |
| Large pages | 15.5710, 15.4589, 15.6210 | 15.5710 |
| Large pages plus write-prefetch | 15.1725, 15.0475, 14.9935 | 15.0475 |

The combination saves 5.79% against the median control and beats every paired
control. One control run was slower than the others; all runs are retained.
The earlier 15.9904 ms control gives a similar comparison. These are three
process repetitions, not a statistical confidence bound across machines.

The same three-repeat check at smaller K used the existing default tiles:

| K | Control (ms) | Combined (ms) | Interpretation |
|---|---:|---:|---|
| 2^16 | 0.65235 | 0.66471 | 1.90% slower; retain the control |
| 2^18 | 3.34262 | 3.27479 | 2.03% faster; smaller benefit |

The optimization should therefore be size-dependent, not enabled universally.

## What changed

`generate_pages.py` produces isolated implementations that request
`MADV_HUGEPAGE` and, where available, `MADV_COLLAPSE` during workspace creation.
The requests cover only whole pages within owned allocations. They do not
change system settings, caller-owned buffers, the route, or the inner circuit.
An unsuccessful request leaves ordinary memory available for encoding.

Both requests succeeded for both allocations in every final large-page run.
Synchronous collapse is preparation work, not encoding work. The existing
benchmark's `setup_ms` field measures the code constructor only, so it does not
measure that workspace preparation cost either.

The generated implementations emit advice diagnostics during setup. They are
experimental opt-in variants, not a new production default or a revised
distance certificate. Exact-map routing changes introduce no new mathematical
distance loss; they also do not complete the asymmetric inner's pending
outward certificate.

The combined variant also changes the existing 32-position read-prefetch to a
write-prefetch hint. The scatter loop and its destination indices are unchanged.
It adds no hot-path allocation or runtime dispatch.

## Why the other changes lost

The ordinary-page sweep covered tile sizes from 64 through 8192 rows and both
packed-24 and 32-bit schedules. The existing 4096-row, packed-24 choice won.
Small tiles improved local processing but made initial bucket routing expensive.

An instrumented control separates the existing stages approximately:

| Tile rows | Inner plus bucket routing | Tile scatter | BCH |
|---|---:|---:|---:|
| 512 | 10.27 ms | 5.54 ms | 2.52 ms |
| 2048 | 8.16 ms | 6.28 ms | 2.79 ms |
| 4096 | 6.79 ms | 6.63 ms | 2.93 ms |

Instrumentation perturbs execution, and these entries average warmups and
measured calls. They explain the tradeoff; they are not the final benchmark.

Single-process, 31-call screens rejected the following variants:

- Eight-way unrolled scattering: approximately 18.8–21.3 ms.
- Gathering within coarse buckets and immediately feeding BCH: 25.6–26.0 ms.
- Two-level buckets with 64 KiB, 256 KiB, or 1 MiB final subtiles: 22.6–24.8 ms.
- Longer-distance scalar prefetch and first-stage write-prefetch: 17.2–18.5 ms.

The two-level construction preserves each coordinate's final destination.
It first partitions a coarse bucket into contiguous sub-bucket streams, then
scatters each sub-bucket locally and evaluates BCH. It adds an intermediate
pass and 8 MiB of retained packed schedule. That cost outweighed its locality
benefit in these implementations.

Hardware-counter runs compared the control against large pages with a
2048-row tile, using 501 calls. Reported dTLB load misses fell from 34.65 million
to 13.94 million, while instruction counts remained about 83.3 billion.
This supports investigating address translation. The counters include process
setup, and both page backing and tile size changed; they do not isolate one cause.

## Correctness and reproduction

Every screened implementation passes the existing independent dense oracle at
K=2^16, 2^18, and 2^20. Tests cover both index layouts, in-place output and suffix
preservation, alternate tile sizes, compaction, linearity, and a second setup.
All 31-call runs have output hash `b180fa625fba455f`. All K=2^20 final 101-call
runs have output hash `6f9101ba915d544a`, matching the control. Smaller-size
runs match their respective controls as well.

Generate candidates from the repository root:

```text
python workstreams/inner_design/routing_opt/generate.py
python workstreams/inner_design/routing_opt/generate_hierarchical.py
python workstreams/inner_design/routing_opt/generate_pages.py
python workstreams/inner_design/routing_opt/generate_profile.py
```

Use the existing `workstreams/inner_design/CMakeLists.txt`. Generated targets
have names `basis_routeopt_*_bench` and `basis_routeopt_*_test`.
`run_screen.sh` covers the initial alternatives; `run_final.sh` repeats the
confirmed finalists and collects hardware counters. Run scripts sequentially
on the otherwise idle benchmark host. Each benchmark acquires the shared lock.
`run_combined.sh` repeats the combined optimization; `run_sizes.sh` checks the
two smaller supported message sizes with their existing default tiles.

`measurements/` retains timing receipts, page-advice results, correctness logs,
hardware counters, binary hashes, and hashes of compiled source files.
`summarize.py` validates source bindings, correctness receipts, and output hashes,
then writes `RESULT.json`. It is a performance-artifact checker, not a proof verifier.

## Recommendation

Retain the existing bucket algorithm. The combined variant is worth making an
explicit, optional deployment optimization, with ordinary pages as the fallback.
Keep the ordinary implementation at K=2^16, where the combination regressed.
For large pages alone, the 2048-row choice saves 4 MiB for nearly the same time;
the combined winner was measured with the existing 4096-row choice.
Before promotion, check allocation behavior and another machine.
Do not change the permutation ensemble to obtain this modest gain.
