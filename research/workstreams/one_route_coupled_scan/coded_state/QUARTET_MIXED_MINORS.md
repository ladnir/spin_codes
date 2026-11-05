# Mixed minors under the existing quartet setup

2026-10-27. This is an actual-setup lemma, not a new acceptance rule or a
whole-code distance certificate. It is preserved because it permits stronger
local constraints without changing the construction. The current numerical
test finds too little benefit to make this the next proof direction.

## The missing joint constraint

Use J=GF(2^32), the retained nonzero delta, and the accepted multiplication
matrix A=L_c. Setup already requires A and A+delta I to be superregular.
For one quartet, let U and V be the incoming and next-state supports.
Let T mark input coordinates belonging to inactive outer rows. Let Z be a
selected set of zero-output coordinates, contained in U intersect V.

On T outside V, the input and next state are zero, so (A Q_prev)[i]=0.
On T intersect Z, the input and output are zero, so
(A Q_prev)[i]=delta Q_prev[i]. These are homogeneous equations on Q_prev,
but they use different diagonal shifts in different rows.

Suppose every square minor of A+delta D is nonzero for every binary
diagonal matrix D. If Q_prev is nonzero, the joint necessary condition is

    |T intersect ((complement V) union Z)| < |U|.

Otherwise choose |U| of the homogeneous rows and the columns U. Their
mixed-shift minor is invertible, contradicting the support U. This condition
implies the two separate minor constraints used previously.

## How often can the stronger condition fail?

Fix row and column index sets I,J of size k. The shifted minor depends only
on the diagonal choices in I intersect J. The empty and full choices are
already covered by acceptance of A and A+delta I. Thus this pair has

    max(0, 2^|I intersect J| - 2)

additional minors. Over all I,J and k=1,...,4, there are 74 additional
polynomials. The sum of their degree bounds is 224.

Each polynomial is nonzero. Its highest-degree homogeneous part is the
corresponding minor of L_c, known to be a nonzero polynomial from the
retained actual-field witness. The new helper also checks all 74 additional
minors directly at that same witness key.

For raw c sampled from GF(q^4), q=2^32, the polynomial root bound gives
failure probability at most 224/q. Original acceptance has probability at
least 1-280/q. Therefore, for one actually accepted key,

    Pr[mixed-minor failure] <= 224/(q-280).

Call such a key exceptional. The N/4=2^19 accepted keys are independent.
If B is their total exceptional count, the seventh binomial moment gives

    Pr_setup[B>=7]
      <= binom(2^19,7) (224/(2^32-280))^7
      <= 2^(-53488800536375/1099511627776)
       < 2^-48.64778.

This event is defined by setup alone and applies simultaneously to every
message. No additional rejection is performed.

## Accounting for the remaining six exceptional positions

Fix proof parameters theta in (0,1] and ell>0. Write P_old(ell) and
P_mix(ell) for the permitted selected-zero polynomials under the old and
mixed constraints. Their feasibility at Z=empty is identical. Every feasible
polynomial has constant coefficient one, and there are at most four
candidate zero coordinates. Hence

    P_old(ell) <= (1+ell)^4 P_mix(ell).

For a fixed set S of at most six exceptional quartet positions, use the
old constraints at S and the mixed constraints elsewhere. The raw
fixed-input point bound q^-4 still applies at every live quartet: intersecting
with a specified minor event can only decrease its probability. Do not
condition and renormalize the distribution at exceptional positions.
The original global acceptance-density allowance remains unchanged.

Consequently the theta-moment count for S is at most
(1+ell)^(4 theta |S|) times the corresponding all-mixed counting expression.
Summing possible S gives the factor

    H = sum_(j=0)^6 binom(2^19,j) (1+ell)^(4 theta j).

If U_mix denotes that counting upper bound for a specified message family,
including the original acceptance allowance, then

    Pr_setup[bad message in this family]
      <= Pr_setup[B>=7] + H U_mix.

This is a transfer between counting expressions, not a claim that setup
was sampled from the stronger acceptance distribution. For families using
different proof markers, retain the corresponding factor H for each family.

At theta=7/40 and ell=13592/664659, the directed calculation gives

    log2 H <= 28760690499381/274877906944 < 104.631.

Thus a bound U_mix<=2^-146 would suffice for a 40-bit bound on that family.
The setup-exception lemma alone is not a distance guarantee.

## Checks and measured value for the current proof

`quartet_mixed_minors.py` checks all 4096 local support/target triples.
The mixed rule removes 290 selected-zero choices while preserving constant
coefficients. It also checks 5120 exact polynomial-ratio inequalities and
all 74 additional minors at the original GF(2^32) witness. The two directed
exponents above are frozen assertions. No frozen helper was edited.

In the joint-regional FLOAT screen at tau=0.1 and the original markers,
the change improves the r=2100 exponent only from about +23068.31 to
+23067.98, before paying H. That is not useful at this point. Preserve the
lemma, but do not replace the construction or pursue this refinement as
though it had resolved the long-quiet gap.
