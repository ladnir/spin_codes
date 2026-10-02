# Packed GL32: Optimized-Driver Transplant

2026-09-30. The GL32 encoder now takes **6.203 ms** in a minimal driver with
explicitly inlined routing. The matched optimized legacy encoder takes
**5.805 ms**; the retained old binary remeasures at **5.829 ms**. Thus GL32
costs about **0.40 ms, or 6.86%**, over the optimized same-build legacy code.
The stronger construction is close to the 6 ms target, but is not a speedup
over the optimized legacy construction.

This investigation also identified a concrete cause of the earlier context
differences: GCC outlined routing callbacks in the large comparison harness.
Those calls occur inside the statically expanded inner loop. The new driver
forces the fixed-width emitter to inline without changing the code or its
sampling distribution. Earlier measurements remain in `packed_mixer_REPORT.md`.

## What Changed

`packed_driver.cpp` includes the existing `joint.cpp` helpers with its main
disabled. It contains only three selected full-encoder modes:

| Mode | Inner and routing implementation | Outer work |
|---|---|---|
| 0 | Original shared-GF16 R2 `encode` source | Unchanged GFNI BCH |
| 1 | Same operations through a force-inlined emitter | Unchanged GFNI BCH |
| 2 | Force-inlined emitter, without redundant GF16 maps | Full packed GL32 + BCH |

Mode 2 uses the same independent GL32 matrices, canonical eight-coordinate
blocks, coefficient layout, and packed BCH kernel as the prior experiment.
This task does not introduce another permutation or mixer distribution.
The new emitter retains four explicit row lanes, fixed-width SIMD assembly,
and streaming 64-byte stores. It adds no dynamic dispatch or timed allocation.

The driver restores the retained benchmark's vector allocation order and
input alignment. Only routing scratch receives the workspace memory advice.
It allocates and samples the GL32 setup outside timing for all modes, including
legacy controls. Its 16 MiB expanded coefficient table is read only by mode 2.
The original source, production kernels, existing proof inputs, and defaults
are unchanged.

## Static Evidence

Disassembly was inspected before running the new timings. The retained old
shared-GF16 R2 encoder has no calls to its route callback. Its function text
is 15,428 bytes. In the previous large packed executable, the same source
produced a 12,738-byte body with **32 static route-callback call sites**.
Its GL32 encoder had three such call sites. These are static counts, not
counts of calls for one complete encoding.

The outlined callback stores its captured lane state through pointers, and
the calling function materializes a capture object. This creates actual
call/capture overhead rather than only changing instruction alignment.
No pointer aliasing rule or numerical operation had to change to eliminate it.

The minimal driver uses `SPIN_FORCEINLINE` on a concrete emitter's
`operator()`. Disassembly contains no emitter calls. The original-source
control also inlines fully in this smaller translation unit.

| Function in minimal driver | Text bytes | GCC stack-usage report |
|---|---:|---:|
| Original-source legacy control | 15,396 | 3,584 B |
| Force-inlined legacy | 14,852 | 3,584 B |
| Force-inlined GL32 | 12,150 | 3,456 B |

These frames include working arrays and are not spill-byte counts. The
15,396-byte legacy body is close to the retained binary's 15,428-byte body.
The change in callbacks explains a substantial performance artifact. It does
not establish a complete causal decomposition of every earlier timing change:
allocation, function layout, and register allocation also changed between
executables. We therefore use the new matched measurements directly.

## Matched Measurements

Workload: `K=2^20`, rate one-half, 128-bit XOR elements, IMT `(128,19,2)`,
precomputed setup. The shared route's four-row tile has 1,028-element stride.
Each process warms up three times, then measures 101 complete in-place calls.
The table gives the median of four process medians: seeds 1 and 17, each in
forward and reverse process order. All three new modes share one executable.
The retained binary is interleaved as a fourth control.

| Variant | Median (ms) | Process-median range (ms) |
|---|---:|---:|
| Retained optimized legacy binary | 5.829211 | 5.819539--5.837030 |
| Original-source legacy, minimal driver | 5.904256 | 5.898786--5.931898 |
| Force-inlined legacy, minimal driver | 5.805222 | 5.801745--5.814298 |
| Force-inlined GL32, minimal driver | 6.203224 | 6.176253--6.229052 |

No samples or processes were discarded. Matrix/permutation generation,
allocation, validation, input filling, and final checksums are outside the
timed interval. As in the retained benchmark, each measured call consumes
the preceding call's in-place buffer; there is no timed input reset.

A preceding two-seed, 31-call screen instrumented the two phases. These
are separate per-process means including warmups, not a decomposition of
the 101-call medians:

| Variant | Inner + route range (ms) | BCH or packed GL32+BCH range (ms) |
|---|---:|---:|
| Retained optimized legacy | 3.163950--3.167013 | 2.670424--2.681072 |
| Original-source minimal legacy | 3.146890--3.158044 | 2.763925--2.769438 |
| Force-inlined minimal legacy | 3.043682--3.052099 | 2.755113--2.765636 |
| Force-inlined GL32 | 2.677570--2.681257 | 3.513034--3.547716 |

The large packed harness had measured approximately 2.99 ms for GL32's
inner/route phase. The minimal force-inlined driver restores approximately
2.68 ms, close to the initial packed build's 2.67 ms. The current packed BCH
phase is approximately 3.53 ms instead of that initial build's 3.42 ms.
No historical offset is subtracted from the measured total.

Removing GF16 saves about 0.37 ms from the optimized same-build route here.
Packed GL32 adds about 0.77 ms to its BCH phase. These separate phase
measurements explain why the complete stronger encoder is modestly slower
than legacy after fixing the comparison harness's inlining problem.

## Correctness and Conditions

All three modes passed independent full-encoder checks at `K=2^14` and
`K=2^20`, seeds 1 and 17, with dense, sparse, and alternate dense inputs.
The reference independently computes the dense inner transpose, materializes
the route, applies the ordinary-layout scalar mixer when present, and calls
the production XOR BCH reference. It compares the entire in-place buffer,
including the untouched suffix. Inner forward/transpose adjoint checks pass.
The force-inlined legacy mode also matches the original encoder directly.

The GL32 timed final checksums equal the earlier packed implementation's
checksums for each seed and call count. Legacy checksums equal the retained
binary's checksums. Existing exhaustive packed-kernel physical-basis and
local-adjoint checks are documented in `packed_mixer_REPORT.md`; this task
reuses that unchanged kernel object. No new sanitizer campaign was run.

Machine: Peach Ryzen 7950X, GCC 15.2, pinned CPU 15. Driver flags match the
retained driver: `-O3 -DNDEBUG -std=gnu++20 -mavx2 -mavx512f -mavx512vl
-mavx512bw -mtune=znver4 -fno-strict-aliasing`, with the existing BCH and
workspace-routing macros. GFNI is enabled in the separately compiled kernel
object. Setup sampling remains reproducible benchmark PRNG sampling, not
literal ideal uniform setup randomness. All benchmark batches held the three
existing locks, with no concurrent compiler or benchmark job.

## Reproduction

New research sources: `packed_driver.cpp` and `packed_driver_run.sh`, alongside
this report. Remote build: `/tmp/spin-packed-driver-G573ww`.

```sh
bash /tmp/spin-packed-driver-G573ww/packed_driver_run.sh /tmp/spin-packed-driver-G573ww build
bash /tmp/spin-packed-driver-G573ww/packed_driver_run.sh /tmp/spin-packed-driver-G573ww check
bash /tmp/spin-packed-driver-G573ww/packed_driver_run.sh /tmp/spin-packed-driver-G573ww screen
bash /tmp/spin-packed-driver-G573ww/packed_driver_run.sh /tmp/spin-packed-driver-G573ww confirm
```

The script defaults to headers in `/tmp/spin-joint-9n57TT` and the unchanged
kernel object in `/tmp/spin-packed-mixer-NnskDD`; explicit alternatives can be
passed as the third and fourth arguments. It links the production BCH
reference from `/tmp/spin-bch-compare-2e2UPS/BchAvx512.o`.
Those prerequisites' build instructions are in the prior comparison reports.

Ignored local raw records: `tmp/packed-driver-results/`, including full
disassembly, symbols, stack reports, checks, and measurements. The measured
executable SHA-256 is
`6ee10f75d6c52cfb58f4be65b11fe11bead769cc29e951ccfcb0f07a277bcf95`.
Its source SHA-256 is
`5a1cb2fce7a129d897340026696db2af6c59e7abbe5429ca36dd80f5779f6d20`.
The unchanged packed kernel object SHA-256 is
`ebe4ef88a815da7ceb8abfb486fc199dbd6d7a6b27d35d0d227a98a694b3a3b5`.

The subsequent [whole-code closure](packed_mixing/FIRST_CLOSURE.md) certifies
>9.5% relative distance with >36.68 bits of setup-failure margin for this
construction, without changing the implementation or these measurements.
Keep the force-inlined emitter as the stable performance baseline. Compact
coefficient storage or packed-kernel tuning can target the remaining roughly
0.2 ms above 6 ms; that tuning is not assumed in the current timing.
