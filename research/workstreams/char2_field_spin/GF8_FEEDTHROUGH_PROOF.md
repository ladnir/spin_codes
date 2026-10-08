# Nonidentity feedthrough: local proofs and rate-three-quarters screens

This checkpoint tests whether stronger local output mixing closes the gap at nominal rate three quarters.
The target is 10% relative symbol distance and a 40-bit whole-code margin at requested dimension `K=2^20`.
The route remains one independent regional permutation.
No result in this note is a whole-code certificate or a measured runtime.

## Construction and probability space

The outer applies a fixed GF(256) MDS `[128,96]` code to each disjoint message row.
Natural padding gives 10,944 rows, padded dimension 1,050,624, and output length 1,400,832.
The requested rate is `128/171`, approximately 0.74854, rather than exactly three quarters.
The bad-output threshold is 140,083 symbols.

After routing, each packet contains eight field elements.
Fix the source quotient and the incoming state before sampling the current packet's diagonals.
Its input coordinates and state coordinates may be dependent linear forms on that quotient.
Each diagonal entry is sampled independently from GF(256) minus zero.
Write `X=D_input x` and `R=D_state r` for the diagonal-scaled packet and state.
The local update has the form

\[
Y=FX+2AR,\qquad r'=R+CX.
\]

The transfer counts witnesses with nonzero input coordinates on an exact support.
It does not treat those coordinates as independent source variables.
Two-state transfers distinguish actual zero state from actual nonzero state.
Rank-cap transfers additionally bound the dimension of the incoming state-form space.

## Cubic maps with nonidentity feedthrough

The first candidate retains the repaired cubic width-eight, state-four maps, which satisfy `CA=0`.
It changes the feedthrough to `F=I+4AC`.
For an input of weight `w=1,...,8`, exhaustive support checks give the following minimum output weights:

| Branch | Minimum output weights |
|---|---|
| Zero-state birth | 8, 6, 5, 5, 3, 2, 2, 1 |
| Nonzero-state return to zero | 8, 6, 6, 4, 4, 2, 2, 1 |

The return map is `I+6AC`: setting `r'=0` gives `R=CX` in characteristic two.
These minima improve output support but do not by themselves provide a rank charge.

### Elimination without identity pivots

Fix an exact input support `I` and a selected set `Z` of zero output coordinates.
Let `B=F[Z,I]` and let `E=(2A)[Z,:]`.
Exact row reduction of `[B|E]` leaves `rank(B)` input-pivot equations and state-only equations `GR=0`.

First impose the state-only equations using the fresh state diagonal.
For a fixed nonzero state vector whose support admits a kernel vector of `G`, the supported columns contain a circuit.
Consequently their rank is at least `spark(G)-1`.
Independent nonzero diagonal pivots bound the probability by `255^(-(spark(G)-1))`.
If `G` has full column rank, no nonzero state survives.

Suppose the incoming state-form space has rank `r` and the surviving space has rank `H`.
Counting its nonzero vectors gives

\[
\mathbb E[256^H-1]
\le (256^r-1)255^{-(\operatorname{spark}(G)-1)},
\qquad H\le\min(r,4-\operatorname{rank}(G)).
\]

The scalar moment problem bounds the rank potential `tau^(H-r)` for `1<=tau<=255`.
The rank-cap variant solves the same problem for each upper tail of `H`.
Differences of those tail bounds form a nonnegative measure on surviving ranks.
This measure dominates nonnegative continuations that are nondecreasing in the rank cap; it is not an asserted transition law.

Next expose the nonpivot input diagonals and charge the input-pivot equations.
Each pivot uses a distinct fresh scalar and a coordinate required to be nonzero.
This contributes `tau^(-rank(B))` in the rank cover.
The outgoing state lies in the span of the surviving state forms and the unpivoted inputs.
Thus its rank cap is at most `H+|I|-rank(B)`.
For a zero-state birth, the cap is at most `|I|-rank(B)`.

The script enumerates numerical profiles for every input support and selected output-zero set.
Translation orbits only identify profiles after checking the simultaneous coordinate symmetry.
The proof does not substitute identity-feedthrough pivots for `rank(B)`.

### Current cubic screening result

`gf8_feedthrough_screen.py` uses the preceding positive transfer and the existing regional source cover.
At `[128,96]`, its rank-cap version remains negative:

| Class | Screened margin in bits |
|---|---:|
| 64 active outer rows | -1,884.83 |
| 128 active outer rows | -3,775.15 |
| Input density 0.03 | -19,341.09 |
| Input density 0.10 | -23,596.49 |
| Input density 0.14 | -18,003.01 |

All five optimizations reported convergence.
These failed upper bounds do not show that the construction has low distance.
They show that this local change does not close the current proof.

## Full-Cauchy packets

The second candidate starts with a Cauchy matrix of size `8+s`, partitioned as

\[
T=\begin{pmatrix}B&A_0\\C_0&E\end{pmatrix},
\qquad
P=\begin{pmatrix}I&0\\0&E^{-1}\end{pmatrix}T
=\begin{pmatrix}B&A_0\\C&I\end{pmatrix},
\quad C=E^{-1}C_0.
\]

It sets `F=B` and `2A=A_0`.
Nonzero diagonal scalings of the output rows and input columns are also allowed.
The concrete helper validates state sizes 4, 8, and 16.

The normalized packet `P` is **not** fully superregular: its state-to-state block is identity.
The needed rank statements concern selected output rows, or selected output rows together with **all** state rows.
For `z` selected output rows, `w` input columns, and `h` state columns, their ranks are

\[
\min(z,w+h),\qquad \min(s+z,w+h),
\]

respectively.
The first follows from Cauchy minors.
The second survives the invertible change of basis on all bottom rows.
The birth map `B` and return map `B-A_0C` are both superregular.
They therefore give output weight at least `9-w` for nonzero input weight `w`.

Eliminating `min(z,w)` input pivots leaves a superregular state Schur complement with `max(0,z-w)` rows.
When this number is at least `s`, no nonzero state survives.
Otherwise the existing diagonal-state moment bound applies.
This replaces support enumeration by a profile depending only on `w,z,s`.

### Return bound without incoming zero-flat enumeration

The current return refinement treats state return and selected output zeros jointly.
Fix the joint form space `U` before exposing either current diagonal.
Its dimension `d` is at most the incoming state-rank cap plus `w`.
Only the `w` input coordinates are required to be nonzero.
Incoming state coordinates may vanish.

For a fixed vector in `U`, let `h` be its state support size.
If `w+h<=s+z`, the selected equations cannot hold for that nonzero vector.
Otherwise the Cauchy rank formula gives probability at most `255^(-(s+z))`.
The surviving kernel has dimension at most `w-z`.

If a kernel of dimension `j` contains a vector with all designated input coordinates nonzero, it contains at least

\[
L_j=256^j-1-w(256^{j-1}-1)
\]

such vectors.
This is the union bound on the nonzero parts of the `w` coordinate hyperplanes.
The expected number of surviving designated-live vectors is at most `(256^d-1)/255^(s+z)`.
The helper `partial_live_rank.py` solves the two-point moment linear program on

\[
(0,0),\qquad (L_j,\tau^{j-d}),
\quad 1\le j\le\min(d,w-z).
\]

It then maximizes over all allowed `d`.
No incoming zero-flat factor is multiplied into this bound.
No second charge using the same diagonal samples is multiplied into it either.
The positive selected-zero expansion and the pure minimum-output bound give two return bounds; the transfer takes their minimum.
The earlier flat-based return covers remain available as controls.

### Full-Cauchy screening result

The following rows use the two-state cover and the designated-live return refinement.
All sparse optimizations converged.
The density-0.14 runs reached their 65-iteration limit; their displayed values are candidate witnesses, not certified optima.

| State size | 64 active rows | 128 active rows | Density 0.10 | Density 0.14 |
|---:|---:|---:|---:|---:|
| 8 | -875.98 | -2,312.06 | +359.88 | +38,198.35 |
| 16 | +355.66 | -1,041.33 | +386.94 | +38,202.83 |

At the state-16, 128-row witness, retaining all 17 rank-cap states improves the margin to -1,007.68 bits.
That single evaluation is not a separate optimization.
The remaining gap is much larger than this observed cap refinement.

The state-four rank-cap control, before the designated-live refinement, gave -1,882.93 and -3,762.13 bits for the two sparse classes.
Its density-0.10 margin was -2,197.75 bits; density 0.14 was positive.
Adding a mark for the outer minimum input weight did not improve its sparse bounds: the optimizer selected the neutral mark.

Increasing state size improves isolated sparse classes, but these gates do not cover the full sparse-to-middle transition.
An independent width-eight, state-32 screen retains the older conservative return cover.
Its margins at input densities 0.03, 0.05, and 0.07 are -14,880.97, -17,188.51, and -14,180.13 bits; all three optimizations converged.
That geometry screen assumes a full-Cauchy matrix but does not validate a specific packed implementation.
At the density-0.05 witness, its marked packet matrix is approximately

\[
\begin{pmatrix}1&0.814803886\\5.23016\cdot10^{-25}&0.999957828\end{pmatrix}.
\]

The return term is already negligible there.
The logarithmic bound has source contribution -101,295.21 bits and output-threshold contribution +118,462.25 bits.
Its inner-matrix contribution is only +14.24 bits; modal conditioning contributes +900.44 bits.
Thus suppressing returns alone cannot materially improve this particular witness.
The next step is to test alternative packet widths and diagnose the intermediate-density loss before attempting a directed certificate.

## Reproduction and validation

The new scripts are `theory/gf8_feedthrough_screen.py` and `theory/gf8_feedthrough_mds_screen.py`.
The latter accepts `--state 4`, `8`, or `16`; `--two-state` disables rank-cap tracking.
`--no-partial-live` restores the older return interface.
Together with `--no-live-kernel`, it restores the conservative rank-only return charge.
For example:

```text
python research/workstreams/char2_field_spin/theory/gf8_feedthrough_mds_screen.py --state 16 --two-state --density .10 .14 --iterations 65
python research/workstreams/char2_field_spin/theory/gf8_feedthrough_tests.py
python research/workstreams/char2_field_spin/theory/gf8_feedthrough_mds_tests.py
```

The cubic suite has seven tests.
It includes exact elimination ranks, full-support support enumeration, finite-field moment checks, successor-rank inequalities, and continuation-cone checks.
The full-Cauchy suite has five tests.
It checks restricted numerical ranks for states 4, 8, and 16, state Schur minors, branch exclusions, and exact-rational continuation cones.
The separate map suite validates the concrete matrices and their factored forward/transpose formulas.
The numerical gates use floating arithmetic and cover only the listed classes.
All frozen certificates remain unchanged.
