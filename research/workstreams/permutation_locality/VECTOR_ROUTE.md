# Register-based column routing: about 5.40 ms

Follow-up: PHYSICAL_COLUMNS.md removes the local copy through direct
cache-line routing and confirms about 5.10 ms, or 1.89x, for the same distribution.

2026-09-24. The uniform-shared-coordinate g=4,c=4 candidate now encodes
in about 5.40 ms, versus 9.50--9.60 ms for matched production controls.
The speedup is 1.76--1.78x, not yet 2x. The change preserves the candidate's
linear map and permutation distribution; it does not close the missing
rank-four or multiple-group proof cases. Production and the paper are unchanged.

## Profile and implementation

The starting candidate used a 16 KiB local copy before its mapped GFNI
BCH kernel. A phase run measured 3.45 ms for inner/routing and 3.24 ms
for copy/BCH. Sampling 1001 encodes at 499 Hz attributed about 48% of
cycles to inner/routing, 19% to the library copy, and most remaining cycles
to GFNI preparation and output tiles. The sample includes setup and checks,
although the 1001 encodes dominate it. Kernel symbols were restricted;
user-space symbols resolved and no samples were lost.

Two changes help:

1. Replace the 16 KiB library copy with an explicit loop over four 64-byte
   loads and stores per iteration. The same aligned local buffer is reused.
2. Keep each four-element column in registers during reverse inner traversal.
   It produces lanes 3,2,1,0. Assemble them into one 512-bit vector, apply
   the precomputed lane permutation, and issue one non-temporal cache-line
   store. This replaces dynamically indexed 128-bit stores into a stack
   staging buffer followed by 512-bit reloads.

Each column's permutation is represented by eight bytes, expanded to eight
64-bit lane indices for one vector permutation. The new controls occupy
4 MiB at K=2^20. The isolated experiment retains the older lane table for
comparisons; a production integration could choose one representation.
The experiment adds no encode-time allocations, indirect callbacks, or
virtual dispatch. The aligned copy buffer remains in a non-coroutine helper.

We also tested packing BCH outputs into full-vector stores with an in-register
transpose. This matches the original BCH map, but adds no consistent speed
advantage after the routing change. Next-group prefetching was slower in
the initial copy-path screen. Retain the simpler mapped GFNI output kernel
as the working candidate: `random-gfni-vector-mapped`.

## Matched performance

Peach Ryzen 7950X, GCC 15.2, CPU 15, Release/znver4, K=2^20, half-rate
BCH[256,128], IMT(128,19), 128-bit elements, mask seed 2. Each process uses
three warmups and 101 in-place encodes, without resetting the input.
Setup, allocation, initialization, validation, and checksums are outside
timing. Three processes per mode and route seed; reverse mode order in the
middle repetition. All three shared benchmark locks are held, and no
benchmarks run concurrently.

The table reports medians of the three process medians:

| Full encoder | Route seed 1 (ms) | Route seed 17 (ms) |
|---|---:|---:|
| Production control | 9.500325 | 9.599120 |
| Previous uniform-shuffle implementation | 6.511175 | 6.491819 |
| Explicit copy loop only | 6.253284 | 6.281345 |
| Explicit copy and register-based routing | **5.402445** | **5.403888** |
| Same, with packed BCH output stores | 5.394460 | 5.426610 |

The selected variant's process medians range from 5.381817 to 5.438843 ms
for seed 1 and from 5.399179 to 5.412183 ms for seed 17. It reduces latency
by approximately 17% relative to the preceding uniform-shuffle path.
The matched production controls imply a remaining reduction of roughly
11--12% to reach 2x; do not substitute an older, slower paper baseline.

A separate phase run gives inner/routing 2.419453 ms and copy/BCH
3.031179 ms. These are instrumented means including warmups, not additive
predictions from independently timed kernels. The shift makes copy/BCH
the larger remaining phase.

## Validation and reproduction

Both route seeds pass all 24 column permutations on all 512 input basis
vectors. All nine GFNI variants, including packed outputs, match the original
BCH routine on all 131072 input basis vectors. Every full timed encoder
checks its entire in-place output against the materialized-route reference,
including the unchanged suffix.

ASan/UBSan passes for all new variants at K=2^14 and seeds 1 and 17.
The two-column register routes also pass for both seeds. The existing dense
transpose, adjoint, direct-route, and suffix checks pass, as do all four
library tests. These sanitizer runs are correctness checks, not performance
measurements or a full-size sanitizer campaign. Their remote log is
`/tmp/spin-locality-JjY7yR/sanitize-vector-route.log`.

Use `profile_current.sh` to profile the preceding implementation,
`run_copy_packed.sh` for the screen, and `confirm_vector_route.sh` for the
matched comparison and exhaustive basis checks. `run_sanitize.sh` includes
the new variants at K=2^14, both seeds, and both two-/four-column geometries
for the register-based routes. Raw measurements remain under
`/tmp/spin-locality-JjY7yR/measurements` on Peach, outside version control.

The release executable SHA-256 for the confirmed comparison is
`cd38b14e40970ba90a685b3b0e713cf59d7fd0d9131f16fcfd3f20f753bd2761`.
The unchanged production library SHA-256 is
`ecbebbfa16d15b4b5460c8811053720d1dd534f2aa81ea5b82757bebd153f0d7`.

## Next step

Focus engineering on the copy/BCH phase, now about 3 ms. Profile the winning
path and test whether preparing its dense inputs can overlap more memory
traffic or avoid part of the local copy, while retaining the fast systematic
input access. The earlier direct scattered path was much slower, so removing
the copy alone is not a supported optimization.

For the proof, RANK_THREE.md records the verified 40.584-bit union over
one-group ranks one through three. The code-equivalent changes here preserve
that coverage, but rank four and multiple active groups remain open. Both
the full-distance certificate and the 2x target are still required.
