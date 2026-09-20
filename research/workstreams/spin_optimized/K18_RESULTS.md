# K18 direct-routing optimization

The selected K18 specialization uses direct 32-bit routing followed by the
existing four-row BCH kernel. It preserves the half-rate BCH [256,128] outer,
the `(t,s,r)=(128,19,1)` IMT inner, setup sampling, and output ordering.
The full margin remains **50.1890763816 bits** at 10% relative distance and
`K=262144`; no new distance claim is needed.

Matched transpose latency is **1.410–1.424 ms**, down from 1.678–1.694 ms:
**15.9–16.5% less time** than the previous optimized implementation.
The paper reports 2.154 ms on the same host; the new timing is approximately
1.51–1.53 times faster than that separately measured result.

## Final matched confirmation

Ryzen 7950X, CPU 15, GCC 15.2, Release, `znver4` tuning, AVX-512 BCH.
Each entry is the median of three process medians. Each process has three
warmups and 101 timed in-place calls, without input copying or resetting.
Setup and workspace allocation are excluded. All benchmarks run serially;
candidate order reverses in the middle repetition.

| Library | Route seed | Tiled control | Direct K18 | Latency reduction |
|---|---:|---:|---:|---:|
| Transpose-only | 1 | 1.693672 ms | 1.424078 ms | 15.92% |
| Transpose-only | 17 | 1.693061 ms | 1.423146 ms | 15.94% |
| Bidirectional, transpose call | 1 | 1.678093 ms | 1.410313 ms | 15.96% |
| Bidirectional, transpose call | 17 | 1.692119 ms | 1.413729 ms | 16.45% |

These are 128-way bitsliced transpose timings, not forward timings.
Output hashes match across routing variants, libraries, and repetitions.
Single-process K16/K20 control measurements differed by at most 2.23%;
their kernels and retained setup capacities are unchanged. Those controls
do not establish new performance improvements at either length.

All six tests passed in each of five confirmation builds: control, selected,
AVX2-only, AVX-512 detection disabled, and ASan/UBSan. The K18 tests include
small and alternative tile sizes, both layouts, compaction, in-place suffix
preservation, forward/transpose dense oracles, adjoint equality, bit-packed
forward encoding, and unchanged `wideView()` schedules. K16 and K20 also pass.
The object-code audit confines EVEX instructions to the dispatched objects.
The fixed map and inner-kernel headers match byte-for-byte between the
selected and control builds.

## Why this helps

The previous path writes inner outputs into buckets, reads those buckets,
and scatters into a tile before applying BCH. The new path scatters each
inner output directly into the four-row BCH representation of the full
8 MiB routed buffer. BCH then reads this representation sequentially.
This removes an intermediate 8 MiB read and 8 MiB write per encoding.

For inner position i, setup records `packFour(mRoute[i])`. This is the
same outer coordinate as before, with its row and coordinate bits arranged
for the four-row SIMD kernel. It changes the storage representation, not
the permutation. All inner inputs are consumed before BCH writes the output,
so the in-place transpose interface remains valid.

The hot function is separate from K16 and K20. It retains the existing
fixed-map circuits and unrolled inner computation. There are no hot-path
allocations, extra permutations, or runtime polymorphism.

## Candidate screen

The initial screen used route seed 1, two process medians per variant,
and 101 calls per process. These are candidate-selection measurements,
not the final matched confirmation.

| Route | Transpose-only, ms | Bidirectional transpose, ms |
|---|---:|---:|
| Tiled control | 1.6822 | 1.6665 |
| Direct 32-bit | **1.4292** | **1.4121** |
| Direct packed 24-bit | 1.4486 | 1.4584 |
| Direct 32-bit, prefetch 16 | 1.4420 | 1.4367 |
| Direct 32-bit, prefetch 64 | 1.5163 | 1.5163 |

The tiled screen also checked 32, 64, 128, 512, 1024, and 2048 rows
with both index formats. Its best additional single-process measurement
was 1.6652 ms, at 1024 rows with 32-bit indices. No tested tile geometry
approached direct routing. Plain 32-bit routing therefore avoids the
extra unpacking or prefetch instructions used by the slower candidates.

## Memory and interface

At the default tile geometry, workspace remains 9 MiB. The direct path
uses the existing 8 MiB bucket allocation and leaves the tile unused.
The current workspace API retains its original geometry.

| Library | Tiled setup, bytes | Direct setup, bytes |
|---|---:|---:|
| Transpose-only | 3,178,504 | 2,129,920 |
| Bidirectional | 6,881,292 | 7,405,576 |

These are retained capacities after Packed24 compaction. The transpose-only
library discards the old slot and offset schedules. The bidirectional
library must retain its forward schedules, so its net setup grows by
approximately 512 KiB. Forward encoding and `wideView()` are unchanged.

## Reproduction

`SPIN_K18_DIRECT=ON` is the default and applies only at K18 with the
AVX-512 BCH backend. Set it to `OFF` for the tiled control. AVX2 and
two-row tiles retain their existing paths. K16 has its own specialization;
K20 remains tiled.

`k18_screen.sh ROOT` builds and tests the routing variants before timing
them serially. `k18_confirm.sh ROOT` checks the control, selected, AVX2-only,
feature-disabled, and sanitizer builds, then runs matched serial timings.
The scripts expect the native source snapshot at `ROOT/hypercat/native/spin`.
That snapshot is copied, not edited. `k18_summary.py DIRECTORY` validates
matching output hashes and summarizes the results.

Raw samples, test logs, and source/binary hashes stay ignored in
`measurements/k18-screen/` and `measurements/k18-confirm/`.
The submitted paper and artifact are unchanged; no data files were committed.

## Next step

The K18 data-movement improvement is ready to use. The next engineering
target is the larger K20 working set; do not assume direct routing retains
this benefit once the routed buffer grows to 32 MiB. Revisit inner parameters
only as a separate experiment with its own certificate checks.
