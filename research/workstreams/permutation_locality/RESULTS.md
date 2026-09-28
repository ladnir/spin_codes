# First iteration: row grouping is not enough

2026-09-23. See [SECOND_ITERATION.md](SECOND_ITERATION.md) for subsequent
distributions and a corrected ISA comparison. This first harness compiled its
new inner kernels for AVX2, whereas the library's existing kernels use AVX-512;
custom-kernel timing losses here must not be treated as optimized ceilings.
The 2x goal remains active. No new full distance certificate or
production change has been made. See [ANALYSIS.md](ANALYSIS.md) for proof progress.

## Measurements

Ryzen 7950X, Linux, CPU 15, GCC 15.2, Release, znver4 tuning, AVX-512 BCH.
Half-rate BCH [256,128], IMT (128,19), K=2^20, 128-bit elements, in place.
Setup, allocation, page preparation, input initialization, validation, and
checksums are excluded. Three warmups precede each timed sequence; input is
not reset between calls. All benchmarks hold the three shared locks and run
serially. No paper numbers serve as matched controls.

Initial 31-call screens at route seed 1:

| Group | Coordinate shuffle | Tiled schedule | Direct scalar scatter |
|---:|---|---:|---:|
| 1 | independent, original distribution | 9.647 ms | 17.204 ms |
| 2 | independent | 9.501 ms | 16.839 ms |
| 2 | shared in group | 8.741 ms | 15.294 ms |
| 4 | independent | 9.324 ms | 16.359 ms |
| 4 | shared in group | 8.661 ms | 14.707 ms |
| 8 | independent | 9.163 ms | 16.043 ms |
| 8 | shared in group | 8.542 ms | 14.643 ms |
| 16 | independent | 9.257 ms | 16.124 ms |
| 16 | shared in group | 8.851 ms | 14.179 ms |

These single-process screening medians do not establish rankings between nearby
group sizes. At g=1, the route matches the public seed-to-route hash exactly.

The second experiment specializes shared g=4. Its four values belong to one
64-byte destination line. The new kernels assemble that line in a small aligned
array and issue two 256-bit ordinary or non-temporal stores. The latter includes
a store fence before BCH reads the buffer.

Three 101-call processes per mode, reversing order in the second pass, seed 17:

| Implementation | Median of process medians |
|---|---:|
| Original distribution, existing tiled implementation | 9.651189 ms |
| Shared g=4, existing tiled implementation | 8.716685 ms |
| Shared g=4, direct scalar scatter | 13.805927 ms |
| Shared g=4, assembled cache-line ordinary stores | 15.749466 ms |
| Shared g=4, assembled cache-line streaming stores | 9.874887 ms |

The tiled gain is 9.68% less time, about 1.11x throughput. The 2x target requires
at most 4.8256 ms against this control. Streaming stores recover much of the
direct-scatter regression but do not beat the tiled candidate.

All modes allocate the same reference and candidate buffers before timing.
Direct scratch receives the library workspace's large-page advice policy and
is aligned to 64 bytes in the second experiment. Retaining diagnostic buffers
makes this a comparison harness, not a memory-minimal deployment implementation.
The original library and its tiled schedules remain unchanged.

## Validation

All four existing library test targets pass. At K=2^14, seed 17, each tested
group/shuffle combination passes the dense transpose oracle, forward/transpose
adjoint check, direct-route equivalence, and in-place suffix check. The two
cache-line kernels also pass these checks for shared g=4. Each timing process
compares the candidate to its materialized reference at the measured size.
Checksums match across implementations for each fixed route and call count.

This is not yet a sanitizer campaign or exhaustive testing across lengths.
Experimental length/ISA restrictions do not change public library support.

Build this directory's CMake project with SPIN_BUILD_TESTS=ON. Run the tests in
BUILD/spin, then `run_screen.sh ROOT` or `run_lines.sh ROOT`. The scripts expect
ROOT/build/locality and retain ignored measurements. SPIN_LOCK_WAIT sets a
bounded lock wait. Raw CSVs are currently under `/tmp/spin-locality-JjY7yR` on
Peach. No raw data have been committed.

Final executable SHA256:
`e3c70295e76ffe62c641aaeea034ba8c85580696f98bc845553303914065d5d2`.
Linked library SHA256:
`ecbebbfa16d15b4b5460c8811053720d1dd534f2aa81ea5b82757bebd153f0d7`.

## Next iteration

Investigate two-dimensional routing tiles and the outer circuit's input layout
together. Aim to remove a routing pass or its random accesses. Keep long-range
spreading explicit: confining whole rows to small chunks or simply increasing
row-group size encounters known cancellation obstructions. Preserve this
iteration as a control; it has not exhausted the new-distribution design space.
