# Exact-region composition boxes

This method covers a fixed occupation Q by partitioning the possible
outer-weight-band counts. It retains the exact region support coefficients
R_j and relaxes only the type choices that remain unresolved inside a box.
The same outer constituent is reused across all rows and regions.

Let band g have the fixed counting-measure bound

    nu_g <= Gamma_g Bernoulli(p_g)^B,
    rho_g = Gamma_g^(1/B).

A composition box specifies integer lower and upper counts l_g,u_g,
with sum_g c_g=Q. Put r=Q-sum_g l_g. Let pi_l be the distribution of the
sum of l_g independent Bernoulli(p_g) variables for each band g. Start with

    S_j = (prod_g rho_g^l_g) sum_i pi_l(i) R_(i+j),  0 <= j <= r.

Apply the adaptive recurrence r times:

    S'_j = max_g rho_g ((1-p_g) S_j + p_g S_(j+1)).

The maximum is entrywise. Each application dominates every fixed choice
of its next band. Ignoring the upper count constraints in this recurrence
is a relaxation; the bounds u_g remain part of the counting factor and the
partition. The final matrix S_0 dominates the density-weighted mixture for
every composition in the box. If r=0, the matrix is the exact reference
mixture for that composition, including its density cost.

For each composition c, the number of ordered band assignments is
multinomial(Q;c). Let C_box upper-bound the sum of these multiplicities in
the box. The occupation contribution of that box is at most

    choose(L,Q) C_box exp(lambda H) e_Z S_0^B 1.

The implementation bounds C_box by the largest multinomial coefficient in
the box times an upper bound on the number of feasible count vectors. The
maximum multinomial coefficient has the most balanced feasible integer
counts: exchanging one unit from a larger free count to a smaller free
count decreases the factorial denominator. Integer water filling finds
that vector. Any choice of all but one coordinate determines the remaining
coordinate, so the product of all box widths except the largest is a
valid upper bound on the number of count vectors.

Each split assigns every feasible count vector to exactly one child. A
parent may retain its own bound or use the sum of its children's bounds.
Different boxes may choose different witnesses; a witness stays fixed
inside its box. Tests enumerate small feasible compositions, verify the
counting bound, compare every fixed-type matrix with the relaxed matrix,
and check that the selected boxes cover every composition exactly once.

All calculations are binary64 diagnostics. This method requires the R_j
coefficients through Q. It complements the coefficient-free typed method
for large L; it does not make a large full coefficient table inexpensive.

## Batched whole dense ranges

`batched_typed_ranges.py` evaluates the typed Cauchy inequality from
`NEXT_GRID_REFINEMENTS.md` over a deterministic partition of counts that
includes zero rows. The root has c_0 <= L-Q_min and sum_g c_g=L, so its
children jointly cover every occupation Q_min..L. Every box has an
independent tilt, band-probability scale, and positive coefficient proposal.
The proposal starts from the mean of the box's feasible vertices; a fixed
log-odds adjustment varies its total nonzero-row probability.

The default outer bands have weights below 3B/16, between 3B/16 and 13B/16,
and above 13B/16. Empty bands are omitted. An all-one singleton keeps its
exact probability-one law. The batched implementation evaluates the same
vertex maxima as the scalar formula and stores every selected witness.

`run_typed_range_grid.py` records both the raw typed margin and the selected
upper bound. If counting all messages is smaller, that fallback is explicit.
A completed cover can therefore remain too loose to support a distance
claim. Such a row establishes evaluation coverage and remains excluded from
positive full-bound parameter choices.
