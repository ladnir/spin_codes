# Uniform certificate for every fixed RM2Sub occupation

## Result and scope

The optimized fixed-\(Q\) thresholds increase with \(Q\). That increase does
not require a sublinear outer block. A single nonoptimized tilt proves a
uniform lower bound on the exponent for every fixed \(Q\ge3\).

Consider the random-outer RM2Sub-S19 ensemble from
RANDOM_OUTER_RM2SUB_RAMP.md. Thus \(N=LB\), \(B\) is even, and \(128\mid L\).
Setup samples independent uniform outer injections, coordinate permutations,
region permutations, and nonzero epoch multipliers. Let \(Z_{D,Q}\) count
messages supported on exactly \(Q\) outer blocks whose encoded output has
weight at most \(D\).

**Theorem 1 (uniform certificate over fixed occupations).** Let
\(\delta=11/100\). For every fixed integer \(Q\ge3\), take any admissible
sequence with \(L\to\infty\) and

\[
  B=9\log_2N+O(1).
  \tag{1}
\]

Then

\[
  \Pr[Z_{\lfloor\delta N\rfloor,Q}>0]
  \le
  N^{-\varepsilon Q+o_Q(1)},
  \qquad
  \varepsilon>0.04439431461547.
  \tag{2}
\]

The existing exact certificates cover \(Q=1,2\), even at the larger target
\(\delta=0.11002\). Consequently, \(B=9\log_2N+O(1)\) closes every fixed
occupation at relative distance \(0.11\).

The quantifier order is important. The theorem first fixes \(Q\), then sends
\(N\) to infinity. It does not control an occupation \(Q=Q_N\) that grows
with \(N\).

## Coefficient bound

This section bounds the exact continuum matrix from
RM2SUB_FIXED_OCCUPATION_CONTINUUM.md. Define

\[
  M:=2^{19}-1,
  \qquad
  r:=\frac1M,
  \qquad
  p:=\frac{2^{18}}M,
  \tag{3}
\]

and

\[
  H:=
  \begin{pmatrix}
    1&1\\
    r&2-r
  \end{pmatrix},
  \qquad
  D_y:=\operatorname{diag}(1,y).
  \tag{4}
\]

For \(Q\ge1\), define nonnegative coefficient matrices \(C_{Q,k}\) by

\[
  D_y(H D_y)^Q
  =\sum_{k=0}^{Q+1}C_{Q,k}y^k.
  \tag{5}
\]

If \(S_{Q,k}\sim\operatorname{Beta}(k,Q+1-k)\), with the usual degenerate
endpoint conventions, the exact region matrix is

\[
  T_Q(\theta)
  =\sum_{k=0}^{Q+1}
  C_{Q,k}\,\mathbb E[e^{-p\theta S_{Q,k}}].
  \tag{6}
\]

Fix \(y\in(0,1]\) and \(\tau>0\). Set

\[
  G_Q(y,\tau):=
  \sup_{0\le x\le1}
  \left\{
    -\tau x+
    \left(1+\frac1Q\right)
    \ln(1-x+x/y)
  \right\}.
  \tag{7}
\]

For every \(x\in[0,1]\), the definition of \(G_Q\) gives

\[
  e^{-Q\tau x}
  \le
  e^{QG_Q(y,\tau)}
  (1-x+x/y)^{-(Q+1)}.
  \tag{8}
\]

Put \(n:=Q+1\). For \(X\sim\operatorname{Beta}(k,n-k)\), the beta integral
and the hypergeometric identity give

\[
\begin{split}
  \mathbb E[(1-X+X/y)^{-n}]
  &={}_2F_1(n,k;n;1-1/y)\\
  &=y^k.
  \tag{9}
\end{split}
\]

The same equality holds at \(k=0,n\) by the endpoint conventions. Applying
(8) inside (6), with \(\theta=Q\tau/p\), proves the entrywise bound

\[
  T_Q(Q\tau/p)
  \le
  e^{QG_Q(y,\tau)}D_y(H D_y)^Q.
  \tag{10}
\]

Because \(D_y\le I\) entrywise,

\[
  D_y(H D_y)^Q\le(H D_y)^Q.
  \tag{11}
\]

Perron monotonicity applied to (10) and (11) yields

\[
  \lambda_Q(Q\tau/p)
  \le
  \left(e^{G_Q(y,\tau)}\rho(H D_y)\right)^Q,
  \tag{12}
\]

where \(\lambda_Q\) is the Perron root of \(T_Q\).

## One tilt for all \(Q\ge3\)

Choose

\[
  y:=\frac12,
  \qquad
  \tau:=\frac{11}{10}.
  \tag{13}
\]

For \(Q\ge3\), the coefficient \(1+1/Q\) in (7) decreases with \(Q\).
The logarithm in (7) is nonnegative. Hence

\[
  G_Q(1/2,11/10)\le G_3(1/2,11/10).
  \tag{14}
\]

The maximizing point for \(Q=3\) is \(x=7/33\). Therefore

\[
  G_3(1/2,11/10)
  =-\frac7{30}+\frac43\ln\frac{40}{33}.
  \tag{15}
\]

For the matrix in (12), direct evaluation gives

\[
  \rho_*:=\rho(H D_{1/2})
  =1-\frac r4+
  \frac12\sqrt{2r+\frac{r^2}{4}}.
  \tag{16}
\]

Recall the fixed-occupation exponent

\[
  \eta_Q(\delta)
  :=\frac Q2-
  \inf_{\theta>0}
  \left(
    \log_2\lambda_Q(\theta)
    +\frac{\delta\theta}{\ln2}
  \right).
  \tag{17}
\]

Equations (12)--(17), evaluated at \(\theta=Q\tau/p\), imply

\[
  \frac{\eta_Q(\delta)}Q
  \ge
  \frac12-
  \frac{
    \ln\rho_*+G_3(1/2,11/10)+\delta\tau/p
  }{\ln2}.
  \tag{18}
\]

The exact certificate described below proves that the right-hand side of
(18) exceeds

\[
  0.11604381273505275.
  \tag{19}
\]

Thus the strict logarithmic-block threshold is less than

\[
  \frac1{0.11604381273505275}
  <8.617434884556627.
  \tag{20}
\]

The fixed-occupation theorem gives

\[
  \Pr[Z_{\lfloor\delta N\rfloor,Q}>0]
  \le
  N^{Q-9\eta_Q(\delta)+o_Q(1)}.
  \tag{21}
\]

Substitution of (19) into (21) proves (2).

## Exact arithmetic receipt

The script certify_rm2sub_uniform_fixed_occupation.py checks (19) without
floating-point arithmetic. It uses the positive series

\[
  \ln\frac{1+x}{1-x}
  =2\sum_{j=0}^{m-1}\frac{x^{2j+1}}{2j+1}+R_m,
  \qquad
  0<R_m<
  \frac{2x^{2m+1}}{(2m+1)(1-x^2)}.
  \tag{22}
\]

The checker uses \(x=1/3,m=14\) for \(\ln2\), and \(x=7/73,m=10\) for
\(\ln(40/33)\). It rounds the square root in (16) upward on a \(10^{-40}\)
grid and uses \(\ln\rho_*\le\rho_*-1\). Every subsequent operation uses
exact rational arithmetic.

The receipt rm2sub_uniform_fixed_occupation_d11.json records

\[
\begin{split}
  \text{objective per active row}
    &<0.38395618726494724,\\
  \eta_Q(0.11)/Q
    &>0.11604381273505275,\\
  9\eta_Q(0.11)/Q-1
    &>0.04439431461547471.
  \tag{23}
\end{split}
\]

The decimals in (23) display exact rational inequalities from the receipt.
They are not binary64 evidence.

## Scope and completed finite bridge

The new bound is uniform in the numerical value of a fixed \(Q\). The
continuum reduction is not uniform in \(Q\). In particular, its collision
bound contains

\[
  \binom Q2\frac{127}{L-1}.
  \tag{24}
\]

The proof also contains boundary errors and a fixed-\(Q\) Perron asymptotic.
Their dependence on \(Q\) is not controlled in this continuum proof. Thus
Theorem 1, by itself, cannot be summed over all \(1\le Q\le L\).

`RM2SUB_DENSE_OCCUPATION.md` supplies the completed finite bridge. It uses a
Bernoulli conditioning law rather than a uniform continuum approximation.
An exact four-state transfer covers every growing occupation with
\(Q/L\le10^{-4}\), and outward three-state boxes cover the remaining bulk.
Together the two notes give the full random-outer minimum-distance theorem.

Transfer to a structured outer still requires a uniform spectrum or
transfer-weighted comparison for that outer family.
