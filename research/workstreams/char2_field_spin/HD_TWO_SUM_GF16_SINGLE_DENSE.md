# G16 single-dense-group controls

The broad small-support screen left one cohort unresolved: a single
active source group whose fixed row set has arity 3–16. This note retains
the failed common-potential control and describes a full-matrix control.
Both concern the frozen Group16 construction and use existing marks only.
Neither changes native code or claims a whole-code certificate.

The final deterministic arity schedule gives a **−53.22305320340571-bit
FLOAT union** over all 8,192 source positions and 65,399 dense row sets.
The common-potential control, the two-old-mark control, and the four-mark
control remain separate raw artifacts. Directed replay is still required.

There are 8,192 possible active source groups and 65,399 dense row sets.
For this source size, every regional output subset is empty or a
singleton. Its full-pair count is zero under both matchings, so no source
geometry maximum is needed. The fixed-arity Hölder weights \(W_a(w)\) are
those defined in [HD_TWO_SUM_GF16_ARITY_DENSE.md](HD_TWO_SUM_GF16_ARITY_DENSE.md).
They include the row-set counts and dense RS normalizers once.

## Failed inherited common-potential control

The first control keeps the 14 arities separate, uses `weight_vectors`
with zero sparse groups, and retains at least 9 nonempty columns among
16. Its endpoint is bounded using the same common potential as the larger
support screens. No arity receives an additional row-set factor.

| Frozen witness | Union log2 FLOAT screen |
|---|---:|
| `sparse-inherited`: \(z=996248217/10^9\), \(\log_2p=-127/2\) | +2686.0474126448876 |
| `dense-inherited`: \(z=8970988817957597/9007199254740992\), \(\log_2p=-2266954815/33554432\) | +2954.108400323131 |
| `h2-inherited`: \(z=1116892707587883/1125899906842624\), \(\log_2p=-40\) | +7203.669447595483 |

Every arity selects the first witness and still fails. Arity 3 is worst,
at +2686.05 bits; arity 14 is best, at +1985.46 bits. Every peak is at nine
nonempty columns. The run took 0.653 seconds. Separating the arities alone
does not close this cohort, and these positive upper bounds are not
counterexamples to the construction.

The helper is
[hd_two_sum_gf16_single_dense.py](theory/hd_two_sum_gf16_single_dense.py).
Its six tests cover the scalar positive tail, the singleton geometry,
count normalization, witness selection, and resource guards. The raw
artifact `hd_two_sum_gf16_single_dense_three_frozen_float.json` is preserved.

## Full matrix polynomial from the initial zero state

The next control retains both inner states through all 256 regions.
Let \(R_i(b)\) be the existing marked regional upper matrix for occupancy
\(b\in\{0,1\}\). For each dense weight \(w\), compute

\[
Q_w=\sum_{|T|=w}\prod_{i=0}^{15}R_i(\mathbf1_{i\in T}),
\]

where the product follows regional order. The coefficient DP starts at
the identity matrix and appends \(R_i(0)+U R_i(1)\). It counts each labelled
mask once without commuting the matrices.

For each fixed arity \(a\), define

\[
B_{a,0}=W_a(0)Q_0,\qquad
B_{a,1}=\sum_{w=1}^{16}W_a(w)Q_w.
\]

The positive matrix polynomial \((B_{a,0}+YB_{a,1})^{16}\) retains the
chronology of the outer columns. The comparison starts in state zero and
sums both final states. Its arity bound is

\[
\frac{8192\,z^{-734003}}{255}\;
e_0\left(\sum_{j=9}^{16}[Y^j]
(B_{a,0}+YB_{a,1})^{16}\right)\mathbf1.
\]

There is no common potential and no \(1/p\) endpoint factor. The fixed-row
Hölder comparison applies to each mask sequence before this positive
matrix sum. It does not sample new row sets between columns.

Only the two old r=1 marks are allowed in this control:
`499451537/500000000` and `1199/1200`. No mark is fitted. The new helper is
[hd_two_sum_gf16_single_dense_matrix.py](theory/hd_two_sum_gf16_single_dense_matrix.py).
Seven tests pass, including exact enumeration of noncommuting matrix
words, the initial-state restriction, and domination by the earlier
common-potential bound on tiny fixtures. Independent source review found
no issue. At the two old marks the matrix screen takes 0.667 seconds and
has union +104.16509799860285 bits. Arity 3 gives +104.17, arity 4 gives
+51.31, and arity 5 gives −5.68; arities 6–16 pass the −40-bit threshold.
Both the mark and the operator treatment differ from the first control,
so their individual effects are not separated by these results.

## Four fixed rational marks and the selected schedule

The isolated wrapper
[hd_two_sum_gf16_single_dense_marks.py](theory/hd_two_sum_gf16_single_dense_marks.py)
tests exactly `1999/2000`, `2499/2500`, `9997/10000`, and `4999/5000`.
It reuses the unchanged full-matrix primitives, with no optimizer. The
four-bank screen takes 0.730 seconds. Its best standalone union is still
positive, at +355.39978581747823 bits: marks suited to small arities worsen
large arities.

The first new mark closes arities 3–5. Selecting a fixed mark separately
for each arity from the two preserved matrix artifacts gives:

| Arity | Selected exact mark | Failure log2 FLOAT screen |
|---:|---|---:|
| 3 | `1999/2000` | −63.170602248 |
| 4 | `1999/2000` | −57.986580563 |
| 5 | `1999/2000` | −53.279273631 |
| 6 | `1199/1200` | −64.544006130 |
| 7 | `1199/1200` | −116.273994211 |
| 8 | `1199/1200` | −119.577264131 |
| 9 | `1199/1200` | −119.450251030 |
| 10 | `1199/1200` | −119.820424546 |
| 11 | `499451537/500000000` | −164.565400823 |
| 12 | `499451537/500000000` | −175.855142536 |
| 13 | `499451537/500000000` | −177.587314138 |
| 14 | `499451537/500000000` | −179.753178792 |
| 15 | `499451537/500000000` | −172.153131267 |
| 16 | `499451537/500000000` | −165.646193146 |
| Union | Fixed three-mark schedule | **−53.22305320340571** |

Every selected term peaks at nine nonempty columns. The arity classes
partition the fixed row sets, so mark selection within each class occurs
before their union and needs no mark-bank factor. The schedule does not
depend on sampled construction keys. No exact-A enumeration was launched.

`hd_two_sum_gf16_single_dense_selected_float.json` freezes this schedule.
It was assembled from the two pinned matrix artifacts after checking
their scope and hashes. It retains all 14 rows, their exact marks and
moment pins, and both source-artifact hashes. This selection performed no
new numerical evaluation and explicitly requires fresh directed replay.

Reproduce each control into a fresh output path:

```text
python -B -m unittest hd_two_sum_gf16_single_dense_tests hd_two_sum_gf16_single_dense_matrix_tests hd_two_sum_gf16_single_dense_marks_tests
python -B hd_two_sum_gf16_single_dense.py --output NEW_POTENTIAL_OUTPUT.json --max-seconds 20
python -B hd_two_sum_gf16_single_dense_matrix.py --output NEW_MATRIX_OUTPUT.json --max-seconds 20
python -B hd_two_sum_gf16_single_dense_marks.py --output NEW_FOUR_MARK_OUTPUT.json --max-seconds 20
```

All CLIs refuse to overwrite an artifact. The three test modules contain
16 tests. The last wrapper's tests pin the four rational marks and check
that per-arity selection introduces no extra bank factor.

Frozen SHA-256 pins:

| File | SHA-256 |
|---|---|
| `hd_two_sum_gf16_single_dense.py` | `1da7ce59e7b8241b49d7410ac44a9c62d3207e24267c4727318bb3a32c36e0ed` |
| `hd_two_sum_gf16_single_dense_matrix.py` | `5b395eef84e0760c56fd22bc3d3a7ae9bb2e9d060ed9d767517b2e4c3f13b9a0` |
| `hd_two_sum_gf16_single_dense_marks.py` | `1c23700767186d014bcf9c8a85fbb0e56bd6666195ceaa5f18ef722d111a9c8b` |
| `hd_two_sum_gf16_single_dense_three_frozen_float.json` | `d144bd41bcbbf13531a4dbb78476bf24ccf8fc4aa57246959f73a0c7b90c35c5` |
| `hd_two_sum_gf16_single_dense_matrix_two_old_marks_float.json` | `ec16d6fe548e8bca9e5486663acbb60013bac5e8edf66888eb24dd09e4e823f8` |
| `hd_two_sum_gf16_single_dense_four_rational_marks_float.json` | `c3f33706ddd4896c68676e370ecef9cfb6f5a85d3a71759caddb5753abcbe9ad` |
| `hd_two_sum_gf16_single_dense_selected_float.json` | `f0bd63b95bed054500796427cc1a7f9f07c96afce911b137ad85f39a0ee44a15` |

The next step is fresh directed replay of the selected matrix schedule.
These FLOAT values must not enter a certified coverage ledger before that
replay succeeds. No production or native work occurred in these controls.
