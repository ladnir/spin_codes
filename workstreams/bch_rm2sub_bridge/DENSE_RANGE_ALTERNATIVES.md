# Dense-range work: stronger shells and fixed row weights

The construction and full-range goal are unchanged. The fixed BCH [256,128]
outer uses RM2Sub t128_s15, fresh independent nonzero field multipliers,
N=2^21 output bits, and bad-weight cutoff H=209716. The verified range
remains Q=1,...,1024, with partial failure contribution below 2^-49.
`verify_coverage_1024.py` rechecks that ledger.

This continuation strengthens the available BCH constraints and develops
an alternative counting bound. Neither change yet certifies Q>1024.

## Exact shell objectives

The earlier strength-29 shell caps are useful but not optimal for the full
retained BCH constraint model. That model includes the audited Q/P sandwich,
common coset spectrum, hull constraints, and independently proved moments.
For a weight w, its physical C-shell objective is q_w+31 h_w.

`exact_joint_shell_caps.py` reconstructs the audited model and changes only
the objective. Each successful solve is checked with rational arithmetic:
all primal constraints, all dual signs, every dual objective inequality,
and equality of the primal and dual objectives. Scaling back and rounding
down gives an integer upper bound on A_w(C).

The new caps currently include:

| Weight | Previous cap, log2 | New exact cap, log2 |
|---:|---:|---:|
| 50 | 69.129 | 60.775 |
| 52 | 75.855 | 60.982 |
| 54 | 77.122 | 62.370 |
| 56 | 78.438 | 64.660 |
| 58 | 79.808 | 67.676 |
| 60 | 81.237 | 70.833 |
| 62 | 82.733 | 73.664 |
| 128 | 126.335 | 124.673 |

The logarithms are diagnostics. Each `joint_shell_exact_w*/cap.json` stores
the integer cap, exact rational bound, and source hashes. The original
caps and occupancy certificates are unchanged. Weights 50,52,54,128 were
also replayed with the read-only command:

```powershell
python -B workstreams/bch_rm2sub_bridge/exact_joint_shell_caps.py --weights 50 52 54 128 --verify
```

The initial binary64 LP screen reported infeasibility or numerical failure.
It is not a mathematical infeasibility certificate; the exact feasible
primal and dual witnesses supersede that diagnostic.

The first sweep stopped at weight 64 after its 120-second solver limit.
That attempt is preserved without a cap. `continue_exact_shell_sweep.py`
uses new folders, a nearby successful basis, and a longer limit. It keeps
failed attempts and continues other objectives instead of deleting results.

## What the convex experiment established

An entrywise convex sequence can replace a Poisson-binomial region count
by a binomial count with the same mean. To prove the comparison, fix two
Bernoulli probabilities with sum s. Conditional on all other trials, the
expectation of a convex sequence has the form

\[
f(k)+s\bigl(f(k+1)-f(k)\bigr)
 +p_1p_2\bigl(f(k+2)-2f(k+1)+f(k)\bigr).
\]

Replacing p_1,p_2 by their average increases their product. The final term
is nonnegative. Repeated pair averaging therefore increases the expectation
until all probabilities equal their common mean.

`convex_region_majorant.py` constructs a valid entrywise convex upper
sequence by backward inequalities. Exact tests check domination, convexity,
and all mixtures from a small probability grid.

At Q2048, convexification adds little loss for several useful tilts. But
some pure-group bounds remain vacuous even after selecting different tilts.
At Q8192, this particular convex upper sequence is far too loose: a tested
loss exceeds 880,000 bits. Thus this implementation is not a uniform
solution for the dense range. The original adaptive screen at Q8192 also
remains vacuous at the tested tilt. These diagnostics do not show code failure.

A separate precision check at Q2048 compared region coefficients with
192-, 512-, and 1024-bit Arb arithmetic. Their reported logarithms agree.
The tested dense gap is not explained by insufficient Arb precision.

## An alternative that retains fixed row weights

The Bernoulli comparison dominates a fixed-weight row by a product measure.
An alternative is to dominate the region transfer by finitely many positive
exponentials, then average the actual fixed-weight row supports directly.

Fix a Chernoff parameter z in (0,1). Let R_j(z) be the previously proved
3-by-3 region envelope for a region containing j ones. Fix an integer d>=0,
numbers a>=0 and Delta>=0, and nonnegative 3-by-3 matrices P_0,...,P_d.
Write r_i=a+i Delta. Suppose a numerical certificate establishes

\[
R_j(z)\le\sum_{i=0}^{d}P_i r_i^j
\quad\text{entrywise for every }0\le j\le L,
\qquad L=8192.
\tag{1}
\]

Use the convention 0^0=1 if a=0. The matrices P_i are bounding coefficients,
not transition probabilities. Rates above one are allowed. Condition (1)
is the outstanding numerical obligation, not an established property of
the present fixed construction.

Let n=256. Define the matrix polynomial coefficients

\[
C_b:=[u^b]\left(\sum_{i=0}^{d}P_i u^i\right)^n,
\qquad 0\le b\le nd.
\]

Matrix multiplication retains region order; the matrices need not commute.
Let S_C(x)=sum_{w>0} A_w(C)x^w be the nonzero outer weight enumerator.
For any nonnegative upper bound \widehat S(x)>=S_C(x), condition (1) implies

\[
\boxed{\mathbb E[Z_Q]\le
 \binom LQ z^{-H}
 \sum_{b=0}^{nd} e_Z C_b\mathbf1\,
 \widehat S\left(a+\frac{\Delta b}{n}\right)^Q,}
\qquad 1\le Q\le L.
\tag{2}
\]

Here Z_Q and the setup probability space are those of `README.md`.
The polynomial \widehat S(x)=sum_{w>0} U_w x^w from certified shell caps is
one admissible choice. A certified joint linear objective can instead bound
S_C at each evaluation point without maximizing all shells separately.

### Derivation

Fix Q row positions and an ordered sequence of region mode labels
(i_1,...,i_n). Expanding the positive products in (1) gives the matrix
factor P_{i_1}...P_{i_n}. A selected outer row of fixed weight w contributes

\[
\frac1{\binom nw}\sum_{T\subseteq[n],\,|T|=w}
       \prod_{j\in T}r_{i_j}.
\tag{3}
\]

Expression (3) averages its original uniform w-subset support. There is
no Bernoulli substitution. Summing over nonzero row messages replaces it by
the same expression weighted by the exact counts A_w(C). Independent row
permutations make the contribution of Q selected rows the Qth power of
that sum.

For nonnegative numbers r_{i_j}, their order-w elementary symmetric mean
is at most the wth power of their arithmetic mean. One proof averages two
coordinates at a time: with their sum fixed, each elementary symmetric
polynomial has a nonnegative coefficient on their product. Averaging those
two coordinates increases that product, and repeated averaging gives equal
coordinates. Thus (3) is at most

\[
\left(\frac1n\sum_{j=1}^n r_{i_j}\right)^w
=\left(a+\frac{\Delta\sum_j i_j}{n}\right)^w.
\]

The row-counting sum is consequently bounded by S_C at this arithmetic
mean. It depends on the mode sequence only through b=sum_j i_j.
For fixed b, the sum of all ordered matrix products is exactly C_b.
Summing these nonnegative bounds, choosing the Q row positions, and applying
the Chernoff factor proves (2).

This replacement is an inequality on fixed-weight support averages. It
does not assume independence between the coordinates of an encoded row.

### Verification and remaining computation

`exponential_mode_bound.py` implements the coefficient identity for exact
small-instance tests. Tests exhaust noncommuting mode sequences and compare
the bound with direct enumeration of all permuted two-row messages of a
[4,2] toy outer code. Before the symmetric-mean inequality, the expansion
equals that direct enumeration exactly. A further test checks every weight
for 256 nonnegative rate vectors, including rates above one.

To turn (2) into a full-range certificate, the next step is to find useful
mode matrices and rates and verify (1) at every required integer j. That
step cannot be replaced by testing a few region weights or interpolating
between them. Afterward, the coefficient calculation and all occupancy
sums require outward arithmetic. Neither numerical obligation is complete.

All new code, proofs, and generated outputs stay in our integration
directory. No file in the other worktree was changed.
