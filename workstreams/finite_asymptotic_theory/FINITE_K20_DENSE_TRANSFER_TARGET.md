# Finite dense-transfer target at (k=2^{20})

## Purpose

This note records the next proof object for occupations

\[
  65\le Q\le8832.
\]

It is a target, not a proved certificate. The single Bernoulli envelope used
for occupations through 64 maximizes over the BA weight before applying the
inner transfer. At dense occupation that order can charge one unlikely
endpoint weight to every active row. The corrected Holder calculation has the
same endpoint sensitivity at large conjugate powers.

## Exact outer counting measure

Let \(\overline A(w)\) be the exact expected spectrum of one B=240
Golay--BA-3 draw, and let \(g\) be the proved lower bound on
\(\Pr[\mathcal G_{240}^{23}]\). For \(23\le w\le217\), define

\[
  a_w:=\frac{\overline A(w)}{g}.
\]

After the independent route permutation, the expected counting mass of one
active row on a support \(S\subseteq[240]\) satisfies

\[
  \nu(S)\le \frac{a_{|S|}}{\binom{240}{|S|}}.
  \tag{1}
\]

This is a pointwise conditional-expectation bound. Products of (1) are valid
because the finite theorem uses independent conditioned BA draws across rows.

## Weight-band mixture

Partition \([23,217]\) into bands \(I_1,\ldots,I_m\). Choose a binary-rational
probability \(p_b\in(0,1)\) for every band and define

\[
  c_b:=\max_{w\in I_b}
  \frac{a_w}
  {\binom{240}{w}p_b^w(1-p_b)^{240-w}}.
  \tag{2}
\]

Let \(\mu_p\) be the Bernoulli-\(p\) product measure on
\(\{0,1\}^{240}\). Equations (1) and (2) give the pointwise domination

\[
  \nu\le\sum_{b=1}^m c_b\mu_{p_b}.
  \tag{3}
\]

Unlike a single worst-shell envelope, (3) retains which weight band paid each
likelihood factor.

For exactly \(Q\) active rows, expand the product of (3). By symmetry, label
a term by its composition
\(\boldsymbol n=(n_1,\ldots,n_m)\), where
\(\sum_b n_b=Q\). If
\(P_{Q,\boldsymbol n}(d)\) is a valid inner bad-weight upper bound when
exactly \(n_b\) rows use Bernoulli probability \(p_b\), then

\[
\begin{split}
  \mathbb E[Z_{d,Q}]
  \le{}&\binom{8832}{Q}
  \sum_{\boldsymbol n:\,\sum n_b=Q}
  \binom{Q}{n_1,\ldots,n_m}
  \left(\prod_b c_b^{n_b}\right)
  P_{Q,\boldsymbol n}(d).
  \tag{4}
\end{split}
\]

Equation (4) is an exact finite reduction. It includes mixed BA row weights;
it does not replace them by one endpoint likelihood.

## Inner interface

Under the product reference measure in one term of (4), bits are independent
across the 240 regions. Within each region, there are \(n_b\) candidate rows
of type \(p_b\), placed by a uniform permutation of the 8,832 positions. The
finite RM2Sub transfer can therefore be computed from a multitype
hypergeometric recurrence:

1. select how many candidates of every type enter each 128-position epoch;
2. average the exact impulse matrix over the corresponding convolution of
   binomial active-bit counts;
3. multiply the 69 epoch matrices while retaining the zero/live boundary
   state; and
4. raise the resulting region matrix to the 240th power and apply a rational
   Chernoff tilt.

This recurrence defines \(P_{Q,\boldsymbol n}(d)\) without an asymptotic
approximation. Direct enumeration of every composition is too large. The
remaining lemma must compress the positive composition sum in (4).

## Required compression lemma

A successful certificate may use either of the following forms.

**Composition boxes.** Partition the simplex
\(\{\boldsymbol n/Q\}\) into rational boxes. For one fixed set of positive
Collatz and Chernoff witnesses per box, prove that the logarithm of the term
in (4) is convex in the independent composition coordinates. Checking box
vertices then covers every integer composition.

**Positive operator generating function.** Introduce fugacities
\(u_1,\ldots,u_m\) for the type counts. Construct a nonnegative
matrix-valued generating transfer whose coefficient of
\(u_1^{n_1}\cdots u_m^{n_m}\) is the multitype region transfer. Apply a
positive weighted norm before coefficient extraction. If the resulting norm
factorizes by type, the multinomial sum in (4) becomes one scalar power.

The proof must retain the exact finite factors
\(\binom{8832}{Q}\), \(\binom{Q}{\boldsymbol n}\), 69 epochs per region,
240 regions, and threshold \(d=233{,}164\). An exponent with unspecified
\(o(N)\) terms is insufficient at this length.

## Acceptance criterion

The dense verifier must cover every integer \(Q\in[65,8832]\) and sum their
outward upper bounds. Together with the existing receipts, it must prove

\[
  \mathtt{0x1.6f3c66666f370p-47}
  +\mathtt{0x1.4dc4b8b530941p-1}\,2^{-84}
  +\sum_{Q=65}^{8832}\mathbb E[Z_{233164,Q}]
  \le2^{-40}.
\]

The occupation-one receipt has 6.479 bits beyond the target. The dense sum
must fit in the remaining probability budget. The dense verifier must also
use the probability space in
`FINITE_K20_INDEPENDENT_SETUP.md`; a reused outer requires a different
higher-moment argument.
