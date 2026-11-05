# Covering profiles with normalized regional bounds

The first adaptive-cell probe loses substantial slack by bounding a marked
coefficient and its inverse normalizer separately. Those two factors can
vary greatly while their ratio changes much less. At a fixed source weight,
we can compute the ratio first and then bound a set of profiles.

This note retains the construction, setup distribution, and all-pair bound
from [the adaptive-cell lemma](HD_FIELD128_ADAPTIVE_CELL_COVER.md). It changes
the numerical enclosure, not the code. The method below covers one fixed
weight; it does not interpolate between weights.

## Normalize before taking maxima

Fix an active-group weight `r`, exact witnesses `z, alpha`, and one directed
all-pair bank. For each feasible full-pair count `f`, compute a nonnegative
matrix `U(r,f)` that bounds the normalized regional matrix `R(r,f)`
entrywise. Use the existing directed extraction, exact matching normalizer,
and entrywise minimum over independently valid positive-marker bounds.
Every feasible `f` must be present in this table.

Let `B` be a native profile cell contained in this selected weight. Let
`E_B` and `O_B` contain all integer even and odd full-pair counts in `B`.
Both sets must be subsets of the feasible table index set
`F_r = {max(0,r-4096), ..., floor(r/2)}`. An integer-empty cell represents
no supports and receives the zero bound without a maximum query.
For a nonempty cell and each matrix entry, define

\[
 V_0(B)_{ij}:=\max_{f\in E_B}U(r,f)_{ij},\qquad
 V_1(B)_{ij}:=\max_{f\in O_B}U(r,f)_{ij}.
\]

The maxima may include full-pair counts not attained by the cell. This
only enlarges the bound. A range-maximum table may answer these queries;
its entries remain directed upper bounds, not interpolated values.

Let `C_B` be the exact projected support-count upper bound from the
adaptive-cell lemma. For the ideal sampled setup, the probability that
a bad message has its active support in `B` is at most

\[
 C_B\left(\frac{\beta(q)^r}{(q-1)z^L}\right)^\alpha
 e_0\bigl(V_0(B)V_1(B)\bigr)^{128}\mathbf1.
\]

Indeed, each profile in the cell has both regional matrices bounded by
the corresponding `V_i`. Products of nonnegative matrices preserve this
entrywise order. The outer factor is constant because the weight is fixed.
Finally, sum over the represented supports, whose number is at most `C_B`.
No independence between the two matching profiles is needed.

At an integer singleton cell, the method uses the same regional matrices
as the selected-weight profile verifier. It therefore avoids the restriction
to one common marker for all four entries of a regional matrix. Directed
rounding may make the displayed rational ceilings differ slightly.

## Coverage obligations

Keep the two singleton-free faces separate at even weights. Use the same
exact support counts and disjoint integer split tree as the adaptive-cell
lemma. A selected-weight claim requires every leaf to pass and a directed
sum of their bounds. An incomplete table, failed cell, or timeout leaves
that weight unresolved. Combining different weights requires separate
valid tables or a further proved enclosure; sampled tables are insufficient.

The implementation and its independent checks are research proof tools.
Their unit tests do not themselves certify any new weight or distance.
