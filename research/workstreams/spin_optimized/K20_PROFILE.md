# K20 profile: three comparable costs

2026-09-19. The current bidirectional implementation spends approximately
one third of its transposed-encoding time in each phase: IMT plus initial
routing, tile scatter, and BCH. The next bounded optimization target is
BCH register pressure, not another tile-size or prefetch sweep.

This is a diagnostic result, not a new speedup. The unchanged implementation
still takes approximately 9.1 ms. It retains BCH [256,128], IMT
(t,s,r)=(128,19,1), 10% relative distance, and its 50.062-bit certified margin.
Neither the submitted artifact nor the paper was modified.

## Measurements

The workload is the in-place transpose on 128 parallel binary instances at
K=2^20. Peach supplies one Ryzen 7950X core, pinned to CPU 15, with GCC 15.2,
Release, znver4 tuning, and the runtime-selected AVX-512 BCH backend.
Setup, allocation, input initialization, and output hashing are outside
the measured phases. All benchmark processes run serially under both locks.

The diagnostic build places clocks and hardware-counter reads at phase
boundaries. It does not instrument individual inner steps or BCH calls.
Each process has three excluded warmups and 101 measured in-place calls,
without input reset. The table averages phase means from four processes:
two repetitions for each route seed, 1 and 17. The coefficient seed is 2.

| Tile rows | Inner + bucket routing | Tile scatter | BCH | Sum |
|---|---:|---:|---:|---:|
| 512 | 3.685 ms | 2.537 ms | 2.830 ms | 9.052 ms |
| 1024 | 3.538 ms | 2.700 ms | 2.885 ms | 9.123 ms |
| 2048, default | 3.081 ms | 2.946 ms | 2.920 ms | 8.947 ms |

Smaller tiles improve scatter but increase the number of initial bucket
streams. Going from 2048 to 512 rows saves about 0.41 ms in scatter while
adding about 0.60 ms to the initial phase. BCH changes comparatively little.
This explains why shrinking the tile did not improve complete encoding in
the earlier [routing screen](K20_RESULTS.md).

Instrumentation can change compiler scheduling as well as add measurement
overhead. At the default geometry, the median of process medians is
8.943 ms instrumented versus 9.127 ms uninstrumented for the bidirectional
library. The transpose-only library moves in the other direction:
9.532 ms instrumented versus 9.044 ms uninstrumented. Therefore the phase
numbers diagnose scale and tradeoffs; they are not replacement performance
numbers or evidence of a new implementation gain.

## Hardware counters and unchanged-binary sampling

Each phase reads one group of four per-thread, user-space hardware counters.
All measured groups ran without multiplexing. Default-geometry phase means
for the bidirectional library are:

| Phase | Instructions, millions | Cycles, millions | Instructions/cycle | Generic cache misses |
|---|---:|---:|---:|---:|
| Inner + bucket routing | 36.320 | 13.723 | 2.647 | 26,052 |
| Tile scatter | 33.555 | 13.083 | 2.565 | 1,989,019 |
| BCH | 13.026 | 13.035 | 0.999 | 286,314 |

The cache counts are Linux generic `PERF_COUNT_HW_CACHE_MISSES` events.
They are not labeled as DRAM transactions or misses at a particular cache
level. Read and clock overhead contributes a small amount of work to the
counter intervals. These data do not isolate memory latency from bandwidth,
instruction dependencies, or execution-port pressure.

A separate 999 Hz cycle-sampling run uses the unchanged bidirectional binary,
three warmups, and 501 measured calls. BCH receives 32.84% of sampled periods;
the fused inner/routing/scatter function receives 66.12%. Setup receives
0.67%. This supports the approximately one-third BCH share independently of
the instrumented phase boundaries. The sampling profile includes setup;
the phase counters and encoder timing do not.

The subsequent [locality screen](BCH_LOCAL_RESULTS.md) tested this lead.
It found that reducing spills through local output groups or forced input
reloading did not improve the complete encoder. The original kernel remains
selected; the optimization hypotheses below should be read with that result.

## Why investigate BCH register pressure next?

The optimized four-row BCH circuit still has substantial stack traffic:

- Its annotated assembly contains 6,335 instruction lines, of which 4,540
  reference stack addresses.
- Instructions referencing the stack receive 73.02% of its sampled periods.
- Hot locations include 512-bit spill stores, reloads, and XOR operations
  that read stack operands.

These counts include memory-operand arithmetic, not only explicit moves.
Cycle samples can land after the event that caused a stall. Consequently,
73.02% is neither a causal stall attribution nor a recoverable speedup.
Together with BCH's low instruction throughput, it nevertheless identifies
a concrete issue worth testing: too many live vector temporaries.

The existing circuit was already improved by reducing XOR sharing and
processing four rows together. Its generator still retains shared signals
across many outputs. A bounded next experiment should preserve four-row
SIMD and the exact matrix while limiting the lifetime of those signals.
Candidates include locally shared output groups and explicit recomputation
of cheap values instead of keeping them live across distant outputs.

Earlier two-row chunking and global sharing-threshold screens are documented
in [the BCH report](../inner_design/bch_20260919/README.md). They must not be
presented as new candidates. The new question is whether a register-budgeted
four-row schedule improves the complete encoder, not merely its XOR count
or stack-frame size. Any candidate must pass symbolic matrix equivalence,
full-encoder tests, and uninstrumented matched timings before adoption.

For scale, reducing the approximately 2.9 ms BCH phase by 25% would save
about 0.73 ms, or 8% of complete encoding, if other costs stayed unchanged.
That is an arithmetic scenario, not a forecast. This profile provides no
reason to change the mathematical inner or weaken its margin.

## Reproduction and checks

Build the retained uninstrumented control with `k20_check.sh` first. Then run:

```sh
bash workstreams/spin_optimized/k20_profile.sh "$PWD"
python3 workstreams/spin_optimized/profile_summary.py \
  workstreams/spin_optimized/measurements/k20-profile
python3 workstreams/spin_optimized/annotate_summary.py \
  workstreams/spin_optimized/measurements/k20-profile/annotate.txt
```

The diagnostic build requires Linux perf-event access. It is selected only
with `SPIN_STAGE_PROFILE=ON`; ordinary builds contain no profiling code.
Set `SPIN_PROFILE=1` when running its benchmark to collect phase statistics.
The default CMake option is off. The scripts copy the existing upstream
snapshot without changing it and acquire the same benchmark locks as the
other experiments.

All six diagnostic-build correctness tests passed. Complete output hashes
match across tile sizes, libraries, instrumentation, and uninstrumented
controls for matching seeds and iteration counts. The summary checks that
the reported call count equals the measured iteration count and rejects
multiplexed counter groups. The independent sampler reported no lost samples.
Raw measurements, annotated assembly, and source/binary hashes remain ignored
under `measurements/k20-profile/`. No data files were committed.
