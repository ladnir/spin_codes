# Combining complete bounds for each matching profile

One choice of proof parameters need not work well for every support class.
We can use several choices without changing the encoder or its setup
distribution. The choice occurs in the analysis, not in the sampled code.

Fix the original active-group count `r`. Let `(a,b)` record the numbers of
full pairs in the ordered even and odd matchings. Retain every feasible
ordered pair. Let `C_r(a,b)` be the common integer upper bound on the number
of supports in that class, as in the
[midpoint-count argument](HD_FIELD128_MIDPOINT_COUNTS.md).

Fix a finite, nonempty set of proof witnesses. Witness `j` consists of
`0<z_j<1`, `0<alpha_j<=1`, and positive marker pairs.
For this witness, the existing regional argument supplies nonnegative
matrices `R_(j,r,a)` and an outer prefactor `P_j(r)=2^(c_j+d_j r)`.
Here `c_j` and `d_j` are the constant and slope of the
[lift prefactor](HD_FIELD128_LIFT_PROOF.md), including its projective divisor.
Write `e_0=(1,0)` and `1=(1,1)^T`. The
[selected-weight transfer](HD_FIELD128_MIDPOINT_BOXES.md) bounds the
contribution of class `(a,b)` by

    B_j(r,a,b) = C_r(a,b) P_j(r)
                e_0 (R_(j,r,a) R_(j,r,b))^128 1.

Let `E_(r,a,b)` be the event that some nonzero message with outer support
in this class has an output below the fixed weight threshold.
Probability is over the construction's fixed setup distribution.
Each witness proves `Pr[E_(r,a,b)] <= B_j(r,a,b)`, so this probability
is also at most `min_j B_j(r,a,b)`.
A union bound over the full partition of supports gives the selected-weight bound

    sum_(a,b) min_j B_j(r,a,b).

There is no additional union factor for the number of witnesses:
we are choosing among valid deterministic upper bounds on one contribution,
not taking a union of distinct failure events. Independence between these
bounds is neither assumed nor needed.
The support classes partition supports; their failure events need not be disjoint.

The minimum must follow each witness's complete chronological transfer and
its own outer prefactor. Taking matrix entries from witnesses with different
`z` or `alpha` would not implement this inequality. The common support count
remains outside the fractional power. Every regional matrix and prefactor
uses the original `r`; complementation is permitted only inside the support
count calculation.

The directed producer builds each coefficient bank afresh, retaining only
its regional upper endpoints and prefactor. It releases that bank before
building the next. It then evaluates complete class bounds with directed
arithmetic, takes minima of their exact upper endpoints, and sums them.
Saved row ceilings already include their prefactors. A separate 384-bit
sum of those ceilings must not add a second prefactor.

A complete result covers one weight only. An audit of saved rows checks
their partition and sum; it does not itself replay coefficient generation
or regional transfers. Numerical screens and results at nearby weights
are not interval certificates. The construction's distance claim still
requires a directed union over all active-group counts.
