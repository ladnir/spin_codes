# Independent paired-map certificate review

The completed receipt `../proof/disjoint-pairs-whole-v2-p256.json` passes this
review. Its whole-code first-moment upper bound is below
`2^(-62.84015235334637)`, hence below the requested `2^-40`. This is a bound
over the stated ideal independent setup ensemble for `K=65536`, `N=131072`,
and output weight at most `13107`. It is not a certificate for a particular
deterministic seed or a timing claim.

The review found no blocking mathematical or source-geometry mismatch. The
proof agent ran the complete 256-bit replay. This review independently
checked its assembly and recomputed q=1; it did not repeat all q=2..512
interval calculations.

## Map and implementation interface

The fixed map has the constant row, six linear rows, and nine quadratic
rows: three single monomials and six sums of disjoint monomial pairs. The
fifteen quadratic monomials occur exactly once. Independent truth-table
evaluation agrees with both the proof's feedback columns and the literal
C++ columns and feedback expressions in `T64PairedMap.h`.

The complete 65536-state census gives the declared expansion spectrum.
The rows satisfy `CA=0` with `C=A^T`. The proof regenerates packet restrictions
and birth data from these actual columns; equality with the old selected
map's weight spectrum is not its only premise.

The proof macro comprises two chronological t64 steps with separate fresh
updates and retained state, not one t128 step. Direct placement of 32
physical steps per region agrees with placement of 16 two-step macros for
regional occupancies 0, 1, and 2. State is retained across region boundaries;
all terminal coordinates are summed, with no final flush.

The paired C++ reverse recurrence, update indexing, and first/last-step
handling agree with the scalar recurrence on inspection. The review tests
model packing and maps; they do not execute compiled SIMD instructions.
Compiled scalar/SIMD and adjoint checks remain a separate implementation
validation obligation.

The native outer field family remains admissible: if the reverse symbol
map is `M_c Q`, its forward map is `Q^T M_c^T`. For uniform nonzero `c`, each
fixed nonzero input has an exactly uniform nonzero image. The fixed
invertible coordinate permutation does not change that property. With
independent choices at every symbol and state update, the fixed-message
process required by the first-moment proof is unchanged. No independence
between different messages is needed.

## Occupancy bounds and final sum

The q=1 calculation uses exact expected outer shell counts and a separate
minimum over valid tilts for each support size. This review derived the
RS[16,8] shell counts independently by shortened-support inclusion-exclusion
over the 65536-symbol alphabet, then convolved the four-packet symbol law.
Their total is exactly `2^128-1`.

A fresh q=1 calculation using 32 physical t64 steps per region agrees with
the receipt's macro-based value to relative error below `2^-200`. The exact
rational fold of saved support endpoints does not exceed the saved q=1
endpoint.

On source inspection, q=2 retains both supports in the bivariate regional
recurrence and folds unordered group pairs with the correct off-diagonal
factor. The q=3..512 calculation uses the pointwise uniform outer-measure
majorant, not an independence claim about actual RS packet occupancies.
Its artificial zero output retains the active-group label, as required.
Every integer occupancy from 1 through 512 has a saved positive endpoint.

Independent integer arithmetic reproduces the exact sum of all 512
endpoints, including the saved hexadecimal mantissa and exponent. One
upward rounding to a 256-bit significand reproduces the saved union:

```
32339837625175123994719208883880423090321278423104167124509294407714646079197 * 2^-317
```

An exact rational comparison confirms that this value is below `2^-40`.
The displayed decimal margin is not used for acceptance. All saved local
source hashes match the reviewed files.

## Reproduction and identity

Run the review checks from the repository root:

```
C:/Python314/python.exe -B -m unittest discover -s research/workstreams/k16_codesign_100us/review -p 'test_*.py' -v
```

All 14 tests pass: five native-field mapping tests, four paired-map/interface
tests, and five completed-receipt tests. No benchmark or remote workload was
run in this review.

The s16 receipt audit accepts the `SPIN_PAIRED16_RECEIPT` environment variable
for a fresh replay saved at another path. Its default remains the original
receipt path shown above.

Reviewed SHA256 values:

| File | SHA256 |
| --- | --- |
| `disjoint-pairs-whole-v2-p256.json` | `9610cf7d7ec3225a7908f3ea08a8239882304c42615ae4ce7bed19647b9ea9c2` |
| `disjoint_pair_maps.py` | `ce7c90cab39796fe03e6edc76fddd9ace863f7ae9071e91780f570125c3307c0` |
| `reproduce_disjoint_pairs.py` | `7e62799774b7f1180573fb6291dbc24b6e0f6e573062268bf1851f26a27c6f59` |
| `K16Paired.cpp` | `a4e553919f2ce8142859d7b28d03eae019605ec0d5d7d85393359c9d808965dc` |
| `T64PairedMap.h` | `0ea07a15ca52917aaa286d4410a81a73d09c64037d6cafc91500ccf5cdc25ac7` |

Next, keep this fixed map and randomizer law while checking the final
compiled kernel, then run the serialized timing gate. If a later
optimization changes the map, update independence, or routing distribution,
re-establish the corresponding proof premise before reusing this receipt.
