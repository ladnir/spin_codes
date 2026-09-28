# Direct cache-line routing: about 5.10 ms

Follow-up: TWO_COLUMNS.md uses the same kernel with smaller bundles.
It retains approximately the same speed and verifies all single-group ranks
for that new distribution. Multiple active groups remain open.

2026-09-24. The uniform-shared-coordinate candidate now takes approximately
5.10 ms at K=2^20, versus 9.62--9.64 ms for matched production controls.
This is a 1.89x speedup. Another 5.5--5.7% latency reduction would reach
2x against these controls. The full distance certificate remains open.

## Remove the local copy without changing the code

The preceding implementation stored four consecutive shuffled columns together.
BCH then copied each four-row group into a 16 KiB local buffer before
reading its columns through a permutation table. Profiling still attributed
about 19% of total cycles to the function containing that copy.

The new implementation writes each complete 64-byte column directly into
the canonical BCH layout. Each column contains one 128-bit element from
each of four outer rows. The inner produces those elements in registers;
a vector permutation restores row order before one non-temporal store.
The BCH kernel can then read its inputs directly, without a local copy
or indirect column addressing. Its existing 8 KiB preparation buffer remains.

The tradeoff is four independent cache-line destinations instead of one
contiguous four-line destination. Measurements show essentially unchanged
inner/routing time and lower BCH time. This does not establish that arbitrary
scattered writes are cheap: every store here fills an aligned cache line.

The route still uses the distribution from RANDOM_BLOCKS_GFNI.md:
four-row groups share a uniform permutation of all 256 coordinates,
and each four-column macroregion independently permutes groups and column lanes.
The new layout changes only where intermediate values reside in memory.
It preserves the sampled route, BCH map, IMT map, and complete encoded output.

The hot route uses one 32-bit destination and one eight-byte shuffle control
per column: 6 MiB at K=2^20. It needs no per-element lane array. The experiment
also retains comparison layouts and reference buffers; 6 MiB is not a claim
about its total retained allocation. Encoding adds no allocations or dynamic
dispatch. Production sources and the paper remain unchanged.

## Confirmed performance

Peach Ryzen 7950X, GCC 15.2, CPU 15, Release/znver4. The workload is precomputed
half-rate BCH[256,128], IMT(128,19), K=2^20, with 128-bit elements and mask seed 2.
Each process performs three warmups and 101 timed in-place encodes.
Input is not reset between calls. Setup, allocation, initialization, validation,
and checksums are outside timing. All benchmarks run serially under the three
shared benchmark locks.

Each table entry is the median of three process medians. Mode order reverses
in the middle repetition.

| Full encoder | Route seed 1 (ms) | Route seed 17 (ms) |
|---|---:|---:|
| Production control | 9.637502 | 9.614819 |
| Previous register route with local copy | 5.413987 | 5.423335 |
| **Direct canonical cache-line route** | **5.108006** | **5.089621** |
| Same canonical route, retain local copy | 5.247957 | 5.232719 |

The selected variant's process medians range from 5.077629 to 5.111012 ms
for seed 1, and from 5.082729 to 5.106683 ms for seed 17.
Relative to the preceding implementation, its latency drops 5.65--6.15%.
The matched speedups are 1.8867x and 1.8891x. Do not use an older, slower
paper timing to claim that the 2x target has been met.

A separate instrumented run reports 2.434508 ms for inner/routing and
2.649673 ms for BCH. These are phase means including warmups, not sums
of independently measured kernels. The preceding mapped-layout profile gave
2.443258 ms and 3.013705 ms respectively.

## Rejected variants

Prefetching did not improve the tested implementations. With the old mapped
layout, prefetching one, two, or four groups ahead into L2 increased full
time from 5.441839 ms to 5.715400--5.739575 ms in the screen.
Reducing GFNI output batches from two tiles to one gave 5.402636 ms,
too close to the control to replace the simpler established choice.

With direct canonical routing, prefetching the entire current group into L1
took 5.726620 ms, versus 5.085514 ms without prefetching in that screen.
Keeping the mapped layout, removing the copy, and prefetching the current
group was worse still: 7.581453 ms. These are bounded negative results for
specific schedules, not a claim that prefetching can never help.

## Validation and reproduction

Both route seeds pass all 24 column permutations on all 512 input basis vectors.
All ten generated GFNI variants match the original BCH routine on all
131072 input basis vectors. Every timed encoder compares its entire in-place
output against the materialized-route reference, including the untouched suffix.
The new route also passes full-output checks at K=2^14 and K=2^18,
with both two- and four-column macroregions, for seeds 1 and 17.
Those single-call checks are not additional performance claims.

ASan/UBSan passes for all newly added variants at K=2^14, seeds 1 and 17.
The two-column direct route also passes for both seeds. The retained dense
transpose, adjoint, route, and suffix checks pass, along with all four library
tests. This is a small-size sanitizer campaign, not a full-size one.
The remote log is `/tmp/spin-locality-JjY7yR/sanitize-physical-columns.log`.

Run `run_copy_overlap.sh` for the initial prefetch/tile screen,
`run_physical_columns.sh` for the layout screen, and
`confirm_physical_columns.sh` for the matched confirmation.
`profile_vector.sh` profiles the preceding 5.40 ms implementation.
`run_sanitize.sh` includes all newly tested variants and the two-column
direct route at K=2^14, seeds 1 and 17.

Raw measurements remain on Peach under
`/tmp/spin-locality-JjY7yR/measurements`, outside version control.
The confirmed release executable SHA-256 is
`78f60158b994f8601451fb0b10e6e0e59c1f33ba5e1263e92b0da656c932e6fd`.
The unchanged production library SHA-256 is
`ecbebbfa16d15b4b5460c8811053720d1dd534f2aa81ea5b82757bebd153f0d7`.

## What remains

RANK_THREE.md gives a verified 40.584-bit margin for the union of the
one-group rank-one, rank-two, and rank-three classes. The layout changes
preserve that partial result. Rank four and multiple active groups remain
uncovered; they must jointly fit within the remaining failure budget.

The implementation is now close enough to the speed target that the missing
proof classes deserve priority. For further speed work, target the BCH
bit-transpose and matrix evaluation, now the larger measured phase.
Retain `random-gfni-physical` as the working implementation, not the rejected
copy or prefetch variants. The goal requires both the full certificate and 2x.
