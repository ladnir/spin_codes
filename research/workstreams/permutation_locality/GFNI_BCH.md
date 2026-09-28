# GFNI tiles reduce the full encoder to about 6.04 ms

Follow-up: RANDOM_BLOCKS_GFNI.md implements the indirect-load proposal below.
A sequential local copy recovers a roughly 6.49 ms uniform-shuffle candidate
with a stronger one-group activation bound. Neither has a full certificate.

2026-09-24. The canonical four-row/four-column route now measures about 1.6x
faster than the matched production encoder at K=2^20. The permutation is
unchanged from CANONICAL_BLOCKS.md; this iteration evaluates exactly the same
BCH map using a different implementation. The full distance certificate for
that permutation and the 2x performance objective remain open. Production
sources and the paper are unchanged.

## Why change the BCH evaluation?

After the routing improvements, BCH consumed more time than inner/routing.
Changing common-subexpression thresholds reduced stack size but not runtime.
The new kernel evaluates small binary submatrices with GFNI rather than
retaining hundreds of intermediate XOR values.

The fixed BchRows matrix has 128 output rows and 256 input coordinates.
Its first 127 columns are systematic. Column 127 is a separate parity
contribution; input columns 128 through 255 form the dense part. The last
output row has no systematic contribution. The generator checks this exact
structure instead of assuming a 128-column identity matrix.

As before, one 512-bit vector holds the same coordinate from four codewords
on 128-bit elements. First transpose bits across groups of eight dense input
coordinates. Each byte then holds the corresponding bit from those eight
coordinates. A GFNI affine instruction with zero additive constant applies
one fixed 8-by-8 binary matrix to every byte. XOR the contributions from
the sixteen input groups, undo the bit transpose, and add the systematic
and parity contributions. The GFNI use is a binary linear operation; it
does not change SPIN to an extension-field code.

The best version computes two groups of eight output coordinates together.
It has sixteen vector accumulators and an 8 KiB prepared-input array.
The output indices are compile-time constants. Two input-group contributions
are combined with each accumulator using a ternary XOR instruction. The
bit-transpose network uses shifts and ternary bit selection. Tile helpers
are non-inlined to retain bounded working storage; the bit-transpose helper
is forcibly inlined. Encoding performs no new heap allocation.

Full unrolling of the input-tile loop did not help in the screen. Keep the
fixed eight-iteration paired-input loop, with explicitly emitted SIMD lanes.
The implementation requires GFNI and AVX512BW in addition to the experiment's
existing AVX512 requirements. The harness checks the added features before
calling these variants; production dispatch and ISA requirements are unchanged.

## Measurements

Peach, Ryzen 7950X, GCC 15.2, CPU 15, Release with znver4 tuning. The workload
is precomputed, half-rate BCH[256,128] with IMT(128,19), K=2^20, on 128-bit
elements. Both reference and candidates encode in place. Setup, allocations,
input initialization, correctness checks, and checksums are outside timing.
Each process performs three warmups and 101 timed calls without resetting
input. Scripts hold all three shared benchmark locks and run serially.

A one-process screen at seed 17 compares the implementation choices:

| BCH implementation, same canonical route | Median (ms) |
|---|---:|
| Existing nonaliasing XOR circuit | 7.258292 |
| GFNI, one output tile | 6.457316 |
| GFNI, two output tiles | 6.382055 |
| Compile-time output tiles | 6.286165 |
| Paired contributions with ternary XOR | 6.123041 |
| Paired loop fully unrolled | 6.165109 |
| Ternary XOR plus cheaper bit transpose | 6.039555 |

The final comparison uses three independent processes per mode and route
seed, reversing mode order in the middle repetition. Mask seed remains 2.
The following entries are medians of process medians:

| Full encoder | Route seed 1 (ms) | Route seed 17 (ms) |
|---|---:|---:|
| Production route and implementation | 9.655060 | 9.609194 |
| Canonical route, nonaliasing XOR BCH | 7.288458 | 7.274041 |
| Canonical route, GFNI with ternary XOR | 6.138580 | 6.122489 |
| Canonical route, GFNI with cheaper transpose | **6.040556** | **6.042750** |

The last row gives speedups of 1.5984x and 1.5902x over the matched
production controls. Its process medians range from 6.018064 to 6.048421 ms
for seed 1, and 6.030307 to 6.058490 ms for seed 17. Relative to the same
canonical route with XOR BCH, the reduction in full-encoder time is about
17%. These are complete encoder measurements, not isolated BCH timings.

The matched 2x target is about 4.8 ms. Reaching it requires roughly another
20% reduction from this candidate, as well as the missing distance proof.

## Correctness and reproduction

`gfni_bch.py` emits build-local C++ and checks all 65536 actions of the
256 submatrices on eight-bit inputs. It also checks the original transpose
network on all 512 basis vectors of an eight-by-64-bit block, in both
directions. The generated maps have zero affine constants.

The compiled `block-gfni-check` mode checks all six GFNI variants against
the original BCH routine on all 131072 input basis vectors of a four-row
batch. All pass, including the revised transpose network. Both computations
are binary linear maps with fixed addressing, so agreement on this complete
basis establishes equality of their linear maps. Each timed full encoder
also compares its whole in-place buffer with the materialized-route
reference, including the suffix that must remain unchanged.

All six new variants pass ASan/UBSan at K=2^14, seed 17. The retained route
checks and all four library tests also pass. These are small-size sanitizer
checks, not a full-size sanitizer campaign. The remote log is
`/tmp/spin-locality-JjY7yR/sanitize-gfni.log`.

Run `run_gfni.sh ROOT` for the basis checks and screen, and
`confirm_gfni.sh ROOT` for the repeated two-seed comparison. The benchmark
driver is the existing C++ experiment; Python only generates the BCH code.
The CMake project now tracks its generator scripts as configure dependencies.
Raw measurements remain under
`/tmp/spin-locality-JjY7yR/measurements/gfni-confirm` on Peach, outside the
repository. No raw measurements were committed.

The confirmed release executable SHA-256 is
`66b8db3303e7930eab0c3b80dc8f2068bbe0af0d61be1f4d0d01146a5c55d47b`.
The unchanged production library SHA-256 is
`ecbebbfa16d15b4b5460c8811053720d1dd534f2aa81ea5b82757bebd153f0d7`.

## Next step

Measure the new phase balance before choosing the next optimization.
Candidate directions include smaller BCH batches that fit better in L1,
and a route layout that reduces the conversion work. Keep the exact map
checks as a gate.

There is also a specific opportunity to improve the proof/implementation
tradeoff. The GFNI kernel loads each dense input coordinate only once into
its prepared array. It might tolerate indirect coordinate loads better than
the earlier XOR circuit did. Test the shared uniform 256-coordinate shuffle
within each four-row group, using contiguous shuffled-order route writes and
an indirect GFNI preparation pass. This is the earlier rectangular g=4,c=4
distribution, not the canonical-block distribution. Its occupied blocks are
sampled from each tuple's bit support, avoiding a fixed canonical-block
enumerator. Measure the cost before investing in that proof. This proposed
variant is not implemented or measured in this report.

The canonical permutation's one-group activation bound
is still the one documented in CANONICAL_BLOCKS.md; no full-distance
probability bound follows from this implementation improvement.
