# M22 closure reduced to one BCH shell

Updated: 2026-09-04.

Under the original RandomStepConv-M22 model, it now suffices to prove

\[
A_{38}(C)\le10^{13}.
\]

With this one assumption, the certified full first-moment bound has about
40.14372 bits. The bounds at weights 40 and 42 are deterministic; the earlier
statistical assumptions at those weights are unnecessary for this result.
The inequality above remains unproved. The full-closure goal is still active.

Without that assumption, the current deterministic BCH envelope gives a
certified full bound of approximately 2^-36.11873. These exponents describe
upper bounds, not estimates of the true bad-setup probability.

## Model and cellwise comparison

Keep C, the route, and the setup defined in RANDOM_INNER_M25_CLOSURE.md, but
use memory 22: the independent matrices are 23 by 23. There are 8192 outer
rows and 256 regions of length 8192. The message length is 2^20 bits, the
output length is 2^21 bits, and the bad-distance cutoff is D=209716.

Fix one nonzero outer row of weight w. Its occupied regions form a uniform
w-subset of the 256 regions. Each occupied region has one input bit at an
independent uniform position. These facts follow from the specified row and
region permutations. Let p_w be the probability that this routed input has
output weight at most D, over those permutations and the inner setup.

For this fixed input, the marginal state process has a simple representation.
Write a_i for the indicator of a nonzero state before position i, and x_i for
the input bit. Independently sample fair bits B_i and bits R_i with
Pr[R_i=1]=1-2^-22. Initialize a_1=0 and set

\[
u_i=a_i\lor x_i,\qquad y_i=u_iB_i,\qquad a_{i+1}=u_iR_i.
\]

A fresh uniform linear map sends each fixed nonzero input vector to a uniform
output vector. Its next-state coordinates and output coordinate are
independent. This proves the representation above for one fixed message.
It does not assert independence between different messages under shared
setup matrices.

Let H=sum_i u_i and Y=sum_i y_i. Conditional on the route and the R_i,
Y has distribution Bin(H,1/2). Partition each region into r=16 cells of
length ell=512. Define two counts on the same route and state coins:

- C_plus charges a cell if its entry state is active or it contains the input one.
- C_minus charges a cell if its entry state is active and every R_i in the cell is one.

Every actually active position belongs to a charged C_plus cell. Every
charged C_minus cell is active throughout. Therefore

\[
\ell C_-\le H\le\ell C_+.
\]

Binomial random variables are stochastically increasing in their number of
trials: couple them by appending independent fair bits. Consequently

\[
\mathbb E[\Pr(\operatorname{Bin}(\ell C_+,1/2)\le D\mid C_+)]
\le p_w\le
\mathbb E[\Pr(\operatorname{Bin}(\ell C_-,1/2)\le D\mid C_-)].
\]

## Exact transition formulas and numerical certification

Put rho=1-2^-22, s=rho^ell, and

\[
\eta=\frac1\ell\sum_{j=1}^{\ell}\rho^j
     =\frac{\rho(1-\rho^\ell)}{\ell(1-\rho)}.
\]

For a cell containing one input at a uniform position, eta is the probability
of an active exit state, regardless of the entry state. Let z mark one charged
cell. Matrix rows and columns index the entry and exit activities, in order
0,1. The zero-input and one-input transition polynomials are

\[
Z_+(z)=\begin{pmatrix}1&0\\(1-s)z&sz\end{pmatrix},\quad
A_+(z)=z\begin{pmatrix}1-\eta&\eta\\1-\eta&\eta\end{pmatrix},
\]

\[
Z_-(z)=\begin{pmatrix}1&0\\1-s&sz\end{pmatrix},\quad
A_-(z)=\begin{pmatrix}1-\eta&\eta\\1-\eta&\eta-s+sz\end{pmatrix}.
\]

In A_minus, an active entry and no reset throughout the cell has probability
s. That event charges one full cell and exits active. The remaining active
exit probability is eta-s, which is nonnegative. An inactive entry never
charges C_minus, even if some positions later become active.

For either sign, the region matrices are

\[
\mathcal Z=Z^r,\qquad
\mathcal A=\frac1r\sum_{j=0}^{r-1}Z^jAZ^{r-1-j}.
\]

Starting with the row vector (1,0), the support-count recurrence appends either
an unoccupied region through \(\mathcal Z\) or an occupied region through \(\mathcal A\).
It records the number of occupied regions and charged cells. Division by
\(\binom{256}{w}\) at the end yields the conditional count distribution for weight w.
Truncating the charge degree does not affect retained coefficients, because
charges never decrease.

The computation uses 192-bit Arb enclosures for these polynomials. For C_plus,
all probability entries and every subsequent product and sum are rounded
downward. For C_minus, they are rounded upward. The binary64 recurrence is
unnormalized and uses only nonnegative operations. No cancellation occurs.
Positive underflow is harmless downward and explicitly rounded up when an
upper bound is required. Final normalization and aggregation use exact
fractions.

No floating binomial CDF is trusted by the certificate. For h trials,
Pr[Bin(h,1/2)<=D] is bounded below by

\[
\begin{cases}
1,&h\le D,\\
1-\exp(-2(D+1-h/2)^2/h),&D<h<2(D+1),\\
0,&h\ge2(D+1),
\end{cases}
\]

and above by 1 for h<=2D, and exp(-2(h/2-D)^2/h) otherwise.
These Hoeffding inequalities follow from
E exp(t(B-1/2))=cosh(t/2)<=exp(t^2/8), followed by exponential Markov and
optimization in t. The bound on cosh follows by integrating
(log cosh x)''=sech^2 x<=1 from zero twice.

The charge cutoff is K=floor((2D+25000)/ell). For the upper bound, all omitted
mass contributes at most the upper binomial bound at ell(K+1). For the lower
bound, omitted terms are simply discarded. Arb encloses every exponential.

The resulting brackets include the 8192 choices of the active outer row:

| Weight | log2 lower bound on 8192 p_w | log2 upper bound on 8192 p_w |
| --- | ---: | ---: |
| 38 | -86.85971 | -85.45846 |
| 40 | -93.97186 | -92.37104 |
| 42 | -101.35840 | -99.52719 |

These logarithms are diagnostics. Exact rational endpoints are retained in
`generated/bch256_q1_activity_cells_outward.json`.

The toy audit independently enumerates 6,400 complete bit-level route/coin
paths, including inactive-state coins that must be marginalized. It checks
both count distributions and the binomial comparison against exact fractions.
A complete directed-arithmetic replay reproduces every saved array and the
exact aggregate.

## The remaining deterministic inequality

Use the OA21 caps from BCH_OA21_REFINEMENT.md. Retain the other deterministic
shell bounds and all higher-occupation certificates. Replace only the
weight-38, weight-40, and weight-42 transfer coefficients by the improved
upper endpoints above. Their complementary shells retain their old certified
coefficients.

After removing the paired weights 38 and 218, the full bound, including all
higher occupations, consumes approximately 0.69834034 of the 2^-40 budget.
Thus the paired weight-38 coefficient leaves an exact one-variable test.
The largest integer cap accepted by the saved bound is

\[
A_{38}\le14,584,006,710,777.
\]

The simpler sufficient cap A38<=10^13 gives the stated 40.14372-bit bound.
The current deterministic cap is 678,661,807,513,927; it does not meet either
threshold. Existing statistical evidence is not substituted for a proof.
This reduction removes two shell assumptions but does not prove the remaining one.

## What exact inner tails alone cannot fix

The saved exact OA21 primal solution maximizing h_38 is feasible for every
current spectrum constraint. Multiplying its A38=31h38 by the certified lower
endpoint for 8192 p38 gives a contribution above 2^-40, with an exponent of
approximately -37.59003. Therefore even the exact inner-tail probabilities
cannot make the first-moment objective of this LP relaxation uniformly below
2^-40. Stronger BCH information or a different proof argument is necessary.

An independent check replaces the binary64 recurrence with integer arithmetic.
It uses four cells per region, rounds region coefficients downward to dyadics
with 64 fractional bits, and propagates masses with 192 fractional bits.
Each multiplication truncates downward by an integer shift; additions are
exact. The event of at most 204 charged cells, followed by the same analytic
binomial lower bound, still gives a feasible-primal contribution above 2^-40,
with exponent approximately -39.19516. This establishes the same obstruction
without trusting the binary64 recurrence.

The feasible primal spectrum need not belong to any code. A lower bound on
its hypothetical first moment is neither a lower bound on our code's bad-setup
probability nor a disproof of the SPIN conjecture.

Verification commands, run sequentially:

    python -B code/certify_bch_oa21_probe.py --verify
    python -B code/certify_bch_q1_activity_cells.py --verify
    python -B code/audit_bch_q1_cell_budget.py --verify

Next: concentrate deterministic BCH work on A38<=10^13. The looser exact
threshold above is available if useful. Further cell refinement may relax
that threshold, but cannot by itself make the current OA21 relaxation suffice.
Do not restart the weight-40/42 statistical experiments for this reduction.
