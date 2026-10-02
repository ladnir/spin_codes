# Packed GL32: Two Versus Four Inner Updates

2026-09-30. Four inner updates cost **6.784301 ms**, versus **6.170696 ms**
for two updates in the same executable. The difference is 0.613605 ms,
or 9.94%, at K=2^20 with precomputed setup and 128-bit XOR elements.
The retained earlier R2 measurement, 6.203224 ms, is unchanged.

A later [exact-map fusion](packed_fused_r4_REPORT.md) reduces R4 to
6.212096 ms against a same-build sequential R4 control of 6.753111 ms.
It preserves the sampled encoder. The measurements below remain the
original update-count comparison, not the current optimized R4 timing.

This is a construction comparison. R4 applies four independent transvections
per IMT step, not two. Its distance proof is tracked separately in the
[hill-climb ledger](packed_mixing/HILL_CLIMB.md). Favorable timings do not
establish a distance certificate, and an R2 certificate does not transfer
automatically to R4.

## Implementation and Correctness

The existing `packed_driver.cpp` now instantiates update counts 2, 3, and 4
at compile time. The command-line choice dispatches before setup and timing.
The default remains R2. All masks, reference calls, adjoint checks, and full
encoders use the selected count. The fixed-width emitter, explicit SIMD
lanes, streaming stores, allocation order, and unchanged packed BCH kernel
are retained. No production default changes.

The forward recurrence is

    Y_i = X_i + A Q_i,
    Q_(i+1) = T_(i,R) ... T_(i,1) Q_i + C X_i.

Thus output uses the entering state; raw-input feedback is added after the
transvections. This order agrees with the proof's birth-class operators.

R2, R3, and R4 each pass independent dense-inner, adjoint, scalar routing/mixer,
production BCH, full-buffer, and preserved-suffix checks at K=2^14 and
K=2^20, seeds 1 and 17. Each check-only process tests dense, sparse, and
alternate dense inputs. Local GCC syntax and mask-index checks also pass.
No new sanitizer campaign was run.

Disassembly contains no calls to the routing emitter. The R2 GL32 encoder
body is 12,150 bytes, identical in size to the retained optimized driver.
R3 is 12,278 bytes and R4 is 12,534 bytes. These sizes are static evidence
that the added template instantiations did not restore callback outlining;
they do not alone establish equal instruction schedules.

## Matched Measurements

Peach Ryzen 7950X, GCC 15.2, pinned CPU 15. Compiler flags and kernel objects
match the [retained driver](packed_driver_REPORT.md). Each process warms up
three times, then measures 101 complete in-place calls. Setup, allocation,
input filling, validation, and checksums are outside timing. No timed input
reset is added between calls. Every benchmark is serial and holds all three
existing benchmark locks; no compiler runs concurrently on the host.

| Seed | Process order | R2 median (ms) | R4 median (ms) |
|---|---|---:|---:|
| 1 | R2, R4 | 6.168716 | 6.759239 |
| 1 | R4, R2 | 6.170621 | 6.785739 |
| 17 | R2, R4 | 6.170771 | 6.796699 |
| 17 | R4, R2 | 6.202619 | 6.782863 |
| | Median of process medians | **6.170696** | **6.784301** |

No samples or processes were discarded. R2's final checksums agree with
the retained measurements for the same seeds and call count.

A separate 31-call instrumentation batch gives these phase means, including
warmups. They are not a decomposition of the 101-call medians.

| Variant | Inner + route range (ms) | Packed GL32 + BCH range (ms) |
|---|---:|---:|
| R2 | 2.649043--2.673122 | 3.535878--3.557683 |
| R4 | 3.247947--3.269453 | 3.522340--3.530846 |

The extra cost is in the inner/route phase. The outer implementation and
route distribution are unchanged. R4 remains near the original performance
target, but this report does not round 6.78 ms down to a 6 ms result.

## Reproduction and Provenance

`packed_variant_run.sh` supplies locked serial build, correctness, matched
measurement, and phase-profile batches. `packed_driver_run.sh` retains its
original default-R2 behavior. Explicit R4 invocation:

```sh
taskset -c 15 <build>/packed_driver 20 2 1 101 0 4
```

Build directory: `/tmp/spin-packed-r4-01sr3y`.
Matched logs: `measurements/updates-confirm-I8Cr5D`.
Correctness logs: `measurements/updates-check-jE2qRV`.
Phase logs: `measurements/updates-profile-1d8cy0`.
Local raw copies are under ignored `tmp/packed-hill/r4-performance/`.

Source SHA256:
`fda1849261c6777d42d13c753452790a3401c33b783174bcc4e07fb29efb3f8b`.
Executable SHA256:
`0eed0757beb8d3f471e4684148d74b603c38ba3a13cd37ac72ebed87c1866adf`.
Unchanged packed-kernel object SHA256:
`ebe4ef88a815da7ceb8abfb486fc199dbd6d7a6b27d35d0d227a98a694b3a3b5`.

Follow-on completed: [fresh whole-code replay](packed_mixing/R4_CLOSURE.md)
certifies >10% distance with >52.05 bits of margin. The
[exact-map fused implementation](packed_fused_r4_REPORT.md) removes most of
the extra inner cost. The measurements above remain the earlier matched
update-count comparison, not the current fastest R4 implementation.
