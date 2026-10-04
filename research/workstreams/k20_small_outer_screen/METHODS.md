# Assembling the K20 proof

The target event is that some nonzero message has output weight at most
D=209715, for the fixed small-outer construction described in the README.
All probabilities are over its ideal independent setup ensemble. The code
has G=8192 outer groups and R=64 routing regions. Let E_q be the event
restricted to messages with exactly q active outer groups.

Every component bounds Pr[E_q], including the union over the choices of
those groups. We can therefore take the smaller valid endpoint for each q,
then sum those endpoints over q=1,...,8192. A complete certificate requires
every occupancy; successful sampled occupancies alone are insufficient.

## Sparse components

The existing q=1 and q=2 engines use exact expected outer support counts.
They average the inner comparison operators over all ordered support
placements, preserving state between regions. Their first moments bound
Pr[E_1] and Pr[E_2]. The ordinary tail engine uses the pointwise outer
comparison and exact placement for each requested occupancy.

## Routing-conditioned components

Fix an active-group set S of size q. The expected outer-word measure is
pointwise bounded by beta^q times independent uniform 256-bit group words,
where beta=2^256/(65535)^8. This domination holds for each fixed route.
It does not assert independence of the actual outer packets.

Condition on the number of potential slots from S in each physical step.
The within-step slot sets are independent uniform subsets given those counts;
the group identities need not be independent. The uniform-word comparison
removes those label dependencies before the local operators are applied.
For a conditional bad-message first moment A, use
Pr[bad | counts] <= min(1,A) <= A^alpha, with 0<alpha<=1.

The [K18 derivation](../k18_margin60/SLACK_REVIEW.md) expands the nonnegative
matrix paths to bound this fractional moment. Uniform four-bit packet
values, including zero, are averaged before powering. Two chronological
physical matrices are multiplied before taking their entrywise power.
The probabilities of the potential slot counts remain outside that power.

Write T_j for the resulting two-step operator when j of its 32 potential
slots are assigned to S. Its regional operator is

\[
 R_q=\frac{[x^q](\sum_{j=0}^{32}\binom{32}{j}T_jx^j)^{256}}
             {\binom{8192}{q}}.
\]

If e_0 selects the initial zero state and 1 sums all terminal coordinates,
the component bound is

\[
 \Pr[E_q]\leq\binom{8192}{q}\,
   \beta^{\alpha q}e^{\alpha\lambda D}\,
   e_0R_q^{64}\mathbf1.
\]

Here lambda is the positive weight tilt used to construct the local operators.
Both lambda and alpha are proof witnesses, not construction parameters.
Matrix powering preserves the state across region boundaries. The outer
group-set union remains outside the fractional power.

## Dense components

Exact coefficient computation becomes costly for large q. Fix a marker
probability p in (0,1) and define

\[
 M_p=\sum_{j=0}^{32}\binom{32}{j}p^j(1-p)^{32-j}T_j,
 \qquad c_q=\binom{8192}{q}p^q(1-p)^{8192-q}.
\]

All matrix-polynomial coefficients are nonnegative. Consequently
R_q <= M_p^256/c_q entrywise. Nonnegative matrix multiplication is monotone,
so the preceding bound remains valid after replacing R_q^64 by
M_p^16384/c_q^64. This argument does not require commuting matrices.
The marker is a coefficient-bound witness, not a change to routing.
At p=1, only q=8192 is covered; then R_8192=T_32^256 exactly.

Each dense receipt evaluates every integer occupancy in its declared range.
It selects among finitely many local tilts, fractional exponents, and marker
probabilities, using outward scalar bounds. There is no interpolation across q.

## Arithmetic and verification

Local physical operators are authenticated dyadic upper endpoints from the
exact s22 census. Global replay uses 256-bit Arb arithmetic, including
fractional powers, marker corrections, and chronological matrix products.
Each final component endpoint is a positive exact dyadic.

The assembler authenticates the source and receipt hashes, checks the
construction and coverage, and selects per-occupancy minima by integer
comparison. It sums all selected dyadics exactly and compares that sum
with 2^-50 or 2^-40. Displayed decimal margins are not used for acceptance.
An independent verifier checks provenance, selections, coverage, and the
exact union; it does not independently regenerate the local census.
