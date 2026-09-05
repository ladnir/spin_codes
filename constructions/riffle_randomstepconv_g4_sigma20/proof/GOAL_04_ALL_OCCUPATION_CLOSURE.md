# Goal 04: all-occupation closure at the 5% threshold

## Objective

Goal 03 bounds selected outer occupation shells.  Goal 04 must include every
number of active data blocks without evaluating 16384 separate shells.

The target is the expected number of nonzero outer words whose RandomStepConv
output has binary weight at most 5% of the output length.  The expectation is
over the local block permutations, the global packet permutation, and the
independent RandomStepConv setup matrices.

## Summing the occupation variable

Fix a packet-support moment parameter $u>0$.  Use the local BCH moments from
Goal 03.  For an outer support of size $r\ge2$, restricted Cauchy gives

\[
Z_r(u)
\le
M_1(u)^{r-2}M_{2,*}(u)M_{2,\mathrm{all}}(u).
\]

There are $\binom Kr$ data supports, where $K=16384$.  Therefore,

\[
\sum_{r=2}^K\binom Kr Z_r(u)
\le
\frac{M_{2,*}(u)M_{2,\mathrm{all}}(u)}{M_1(u)^2}
\left((1+M_1(u))^K-1-KM_1(u)\right).
\]

The $r=1$ contribution is at most $K M_3(u)$ by Goal 03.  Their sum defines
the all-occupation moment $Z_{\mathrm{all}}(u)$.

## Support localization

On one Goal 02 interval $I=[a,b]$, the analytic inner certificate has the
form

\[
B_I(h)=\frac{C_I c_I^h}{\binom Nh}.
\]

Define

\[
f(h):=-\log\binom Nh.
\]

The function $f$ is convex.  Let $A_I+B_Ih$ be the secant through
$(a,f(a))$ and $(b,f(b))$.  Then

\[
\binom Nh^{-1}\le \exp(A_I+B_Ih)
\]

for every integer $h\in I$.  This inequality replaces the reciprocal
binomial by an exponential support factor.

For $0<s\le1$, the lower-tail inequality is

\[
\sum_{H\in I}c_I^H
\le s^{-b}Z_{\mathrm{all}}(c_Is).
\]

For $s\ge1$, the upper-tail inequality is

\[
\sum_{H\in I}c_I^H
\le s^{-a}Z_{\mathrm{all}}(c_Is).
\]

The sums include the expected multiplicity of all nonzero outer words.  Apply
them with the adjusted support base $c_I\exp(B_I)$.  The smaller optimized
tail bound gives the interval contribution

\[
C_I\exp(A_I)
\min\left\{
\inf_{0<s\le1}s^{-b}Z_{\mathrm{all}}(c_Ie^{B_I}s),
\inf_{s\ge1}s^{-a}Z_{\mathrm{all}}(c_Ie^{B_I}s)
\right\}.
\]

The implementation subdivides wide Goal 02 intervals.  It retains the same
valid inner tilts on every subinterval.  Reflection around $N/2$ prevents the
last interval from replacing all binomial denominators by $\binom NN=1$.

## Acceptance criteria

Goal 04 is complete as a diagnostic when the evaluator:

1. sums all occupations $1\le r\le K$ through the closed form above;
2. evaluates both support tails on every refined interval;
3. validates the analytic inner form against Goal 02;
4. identifies the dominant interval and occupation saddle;
5. reports whether the complete first moment is below one at 5%.

The algebraic inequalities are exact.  Floating-point optimization and
evaluation remain diagnostic until a separate outward-rounded pass succeeds.
