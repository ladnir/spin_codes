# Bidirectional port and K16 tuning

The optimized bidirectional target now uses four-row AVX-512 BCH in its
transposed path. Its forward implementation and exported routing tables are
unchanged. The K=2^16 specialization reduces transposed latency by another
14--16% relative to the tiled four-row implementation.

## Final comparison

Ryzen 7950X, CPU 15, GCC 15.2, Release, `SPIN_TUNE=znver4`. Each entry is the
median of three process medians. Each process performs three warmups and
101 timed in-place encodings, without copying or resetting the input.
Backend order is reversed in the middle repetition. All benchmarks run serially.

| Target | Route seed | Tiled four-row control | Direct K16 default | Reduction |
|---|---:|---:|---:|---:|
| Transpose-only | 1 | 0.413182 ms | 0.346677 ms | 16.10% |
| Transpose-only | 17 | 0.411007 ms | 0.346517 ms | 15.69% |
| Bidirectional, transposed call | 1 | 0.403634 ms | 0.344063 ms | 14.76% |
| Bidirectional, transposed call | 17 | 0.407120 ms | 0.348650 ms | 14.36% |

These timings process 128 parallel binary instances. Setup and workspace
allocation are excluded. This is not a forward-encoding speedup claim.

The larger sizes still use the tiled path:

| Target | K | Control | Build with K16 specialization |
|---|---:|---:|---:|
| Transpose-only | 2^18 | 1.688682 ms | 1.685347 ms |
| Transpose-only | 2^20 | 9.189947 ms | 9.178315 ms |
| Bidirectional, transposed call | 2^18 | 1.674456 ms | 1.682912 ms |
| Bidirectional, transposed call | 2^20 | 9.067478 ms | 9.095010 ms |

All differences at larger sizes are within 0.51% in this comparison. We treat
those small differences as comparable performance, not demonstrated speedups.
The direct routine has a separate function, avoiding changes to the larger-size
loop's register allocation. It does not allocate direct tables at larger sizes.

## Why the small path helps

The tiled path first stores the inner's output in buckets, then scatters those
buckets into a tile before applying BCH. At K=2^16, all routed values occupy
2 MiB. The direct path scatters once into that full buffer and then applies BCH.
It removes 2 MiB of reads and 2 MiB of writes for the intermediate copy.
One 32-bit destination stream also replaces two packed24 index streams.

For inner coordinate i, setup stores `packFour(mRoute[i])`. This sends each
emitted value to the same outer coordinate, represented as a four-row SIMD lane.
The four-row BCH kernel reads that representation and produces the original
row-major output. Neither the mathematical permutation nor the inner changes.

After Packed24 compaction, retained setup changes as follows:

| Target | Tiled four-row control | Direct K16 |
|---|---:|---:|
| Transpose-only | 794,632 bytes | 532,480 bytes |
| Bidirectional | 1,720,332 bytes | 1,851,400 bytes |

The bidirectional target preserves its original tables for forward encoding
and `WideView`. The direct table replaces its extra packed BCH table but is
slightly wider. Both targets keep the existing 3 MiB workspace geometry;
the direct transposed path does not use its 1 MiB tile allocation.

## Alternatives examined

A preliminary tile sweep tested 32, 64, 128, 256, and 512 rows with packed24
and 32-bit indices. Its best single-process result was approximately 0.39 ms,
slower than direct routing. That sweep is screening evidence, not a repeated
ranking of all tile choices.

Repeated direct-route screens used both route seeds and three processes per
case. Packed24 direct indices took approximately 0.35--0.36 ms. Plain 32-bit
indices were faster; adding store-prefetch leads of 16 or 64 did not help.
We therefore selected the simpler, no-prefetch 32-bit loop. The final comparison
above measures that selection after separating it from the larger-size routine.

## Checks and scope

The final control, default, AVX-512-disabled, and feature-masked builds each
pass three correctness suites. All three also pass under ASan and UBSan in
the enabled default build. Tests cover:

- dense transpose references at K=2^16, 2^18, and 2^20;
- both backends and schedule formats, alternate tiles, and in-place suffixes;
- bidirectional forward references and the adjoint identity at K=2^14, 2^16,
  2^18, and 2^20;
- packed-bit forward encoding after compaction;
- byte-for-byte equality of exported `WideView` tables;
- legacy configurations and rejection of unsupported tile sizes.

Equivalent benchmark inputs produce identical complete-output hashes.
Regenerated local sources match the sources used for the final measurements.
Disassembly checks place EVEX instructions only in the dispatched AVX-512
objects. Feature masking tests fallback selection; no physical AVX2-only host
was used in this campaign.

`final_check.sh` replays the final build/test/benchmark matrix.
`production_summary.py` checks the results; `small_summary.py` checks the
earlier candidate screen. Raw measurements and build files are ignored by Git.
The submitted paper, supplemental archive, and upstream application checkout
have not been modified by this port.
