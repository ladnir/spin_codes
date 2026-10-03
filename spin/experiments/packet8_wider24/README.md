# Byte-packet encoder with a wider outer and 24-bit state

This isolated research encoder evaluates a rate-one-half construction with
precomputed setup and 128-bit XOR elements. It is not a production backend
and does not replace either public code family.

At K=65536, the complete transpose measured 98.1625 microseconds in the
fresh-seed holdout. The matching [complete outward certificate](../../../research/workstreams/packet8_codesign/iteration5/README.md)
gives **68.893748662558 bits** of whole-code margin at 10% distance.
An independent global replay also passes the exact 40-bit threshold.
See the [performance record](../../../research/workstreams/packet8_codesign/iteration5/PERFORMANCE.md)
for every run median and the measurement conditions.

## Construction

Each 256-bit input group is encoded by eight parallel GF(16) RS[16,8] rows.
Their outputs form sixteen 32-bit symbols. Each symbol receives an independent
nonzero-transitive binary randomizer: the adjoint of multiplication by a
nonzero GF(2^32) scalar. The transpose therefore uses ordinary tower-field
multiplication on native byte planes.

Each group supplies 64 byte packets. A group permutation assigns one packet
to each of 64 regions; a separate permutation orders the groups in each
region. These are complete setup-time permutations, not a small heuristic
permutation bank.

The inner consumes eight bytes per step. Over the AES polynomial basis of
GF(256), its 24-bit state is s=(a,b,c), and its maps are

\[
 (As)_h=a+hb+h^2c,
 \qquad Cx=\left(\sum_hx_h,\sum_hhx_h,\sum_hh^2x_h\right),
 \qquad h=0,\ldots,7.
\]

The forward step emits y=x+As, then updates s'=L_r^T s+Cx. Here L_r is
multiplication by a nonzero scalar in GF(256)[z]/(z^3+z+1).
The transpose uses a six-byte-product circuit for L_r. Both A/C adjoints
use the literal binary coordinate dot product; ordinary field multiplication
is not substituted for a binary adjoint.

State starts at zero, persists between regions, and is not flushed.
Setup draws a separate update for every physical step and fixes the resulting
code for reuse. The fixed-state action law needed by the proof is unchanged
when GL(3,GF(256)) updates are replaced by these adjoint scalar updates.
The [implementation audit](../../../research/workstreams/packet8_codesign/iteration5/IMPLEMENTATION_AUDIT.md)
gives the coordinate conventions and distribution-transfer argument.

## Interface and layout

`Plan(K, seed)` accepts positive multiples of 2048 within its checked storage
and routing limits. The benchmark target is K=65536; acceptance of another
size does not assert a distance certificate for that size.

`transposeFast(input, output, scratch, plan)` maps 2K input elements to K
output elements. Input and output need 16-byte alignment; output may equal
input. Scratch must be disjoint, 64-byte aligned, and contain at least
`plan.scratchBlocks()` elements. The hot operation allocates nothing.

`reverseRouteFast`, `outerFast`, and `routeOnlyFast` expose the measured
stages. `Scalar.cpp` supplies independent literal forward and transpose
oracles, including explicit packing helpers. Scratch contains native byte
planes throughout the inner-to-outer boundary; there is no intermediate
unpack/repack pass.

The SIMD experiment requires AVX-512, VBMI, and GFNI. It is a standalone
CMake target and uses an existing `spin::spin` package:

```sh
cmake -S spin/experiments/packet8_wider24 -B build/packet8_wider24 \
  -DCMAKE_BUILD_TYPE=Release -Dspin_DIR=/path/to/spin/cmake
cmake --build build/packet8_wider24 -j
ctest --test-dir build/packet8_wider24 --output-on-failure -j1
```

The executable arguments are

```text
spin_packet8_wider24 [K=65536] [seed=1] [calls=0] [mode=0] [phases=0]
```

Zero calls runs correctness checks only. Modes 0, 1, 2, and 3 select the
complete transpose, route-only, outer-only, and inner-plus-route operation.
Use phases=0 for headline timings; phases=1 inserts a diagnostic clock
between the inner and outer. Calls are capped at one million.

For a matched three-way campaign, `run_ab.sh` runs the certified nibble
control, the earlier byte control, and this candidate sequentially. It
requires a fresh output path and takes the shared encoder benchmark locks:

```sh
bash spin/experiments/packet8_wider24/run_ab.sh \
  /path/to/nibble-control /path/to/byte-control \
  /path/to/spin_packet8_wider24 /path/to/fresh-results.txt \
  2001 15 41 113 257 997
```

Do not run a separate benchmark concurrently. Setup, allocation, and
correctness checks are outside the timed region; each run uses five warmups
and normal pages. Preserve this matched proof/performance checkpoint before
further tuning. It is not a production promotion or a replacement of the
previous certified baselines.
