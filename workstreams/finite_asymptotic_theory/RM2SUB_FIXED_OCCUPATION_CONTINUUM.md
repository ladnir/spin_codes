# Fixed-occupation continuum for the RM2Sub-S19 inner

## Purpose

The occupation-one and occupation-two calculations use different region
matrices. Repeating that derivation for every fixed occupation would not
produce a uniform sparse theorem. This document gives one matrix formula for
all fixed occupations.

The formula also tests the block schedule. At relative distance \(0.11\), the
common-tilt threshold increases with occupation. The binary64 diagnostic first
exceeds \(c=18/5\) at sampled occupation \(Q=16\). This is a limitation of the
current sufficient condition, not a counterexample to minimum distance.

## Continuum reset process

Fix an outer occupation \(Q\ge1\). The symbol \(Q\) counts active outer
blocks. It is unrelated to the RM2Sub field state used in the construction
documents.

Set

\[
  M:=2^{19}-1,
  \qquad
  r:=\frac1M,
  \qquad
  s_0:=1-r,
  \qquad
  p:=\frac{2^{18}}M.
  \tag{1}
\]

Here \(r\) is the termination probability at an impulse entering a live
state. Define the impulse transition

\[
  P:=
  \begin{pmatrix}
    0&1\\
    r&s_0
  \end{pmatrix}.
  \tag{2}
\]

Rows and columns index zero and live state, in that order. A zero state
becomes live at an impulse. A live state terminates with probability \(r\).

Fix \(\theta>0\), set \(\gamma:=p\theta\), and define the interval transfer

\[
  D_\gamma(x):=
  \begin{pmatrix}
    1&0\\
    0&e^{-\gamma x}
  \end{pmatrix},
  \qquad 0\le x\le1.
  \tag{3}
\]

The live entry in (3) is the limiting output-weight transform for a live
interval occupying fraction \(x\) of one region.

## Region matrix for \(a\) impulses

Condition on a region that receives \(a\) impulses. For fixed \(a\), the
region permutation makes their locations a uniform \(a\)-element subset of
\([L]\). As \(L\to\infty\), the normalized locations converge to the order
statistics

\[
  0<U_1<\cdots<U_a<1
  \tag{4}
\]

of \(a\) independent uniform points. Define the spacings

\[
  X_0:=U_1,\quad
  X_i:=U_{i+1}-U_i\ (1\le i<a),\quad
  X_a:=1-U_a.
  \tag{5}
\]

The vector \((X_0,\ldots,X_a)\) is uniform on the \(a\)-simplex.

Define

\[
\begin{split}
  K_a(\theta)
  :=a!\int_{\substack{x_i\ge0\\\sum_{i=0}^a x_i=1}}
  D_\gamma(x_0)P D_\gamma(x_1)P\cdots
  P D_\gamma(x_a)\,d\boldsymbol x.
  \tag{6}
\end{split}
\]

For \(a=0\), equation (6) means \(K_0=D_\gamma(1)\).

**Lemma 1 (fixed-\(a\) region limit).** Fix \(a\). Suppose \(128\mid L\)
and \(B=O(\log L)\). The exact full-state region transfer for \(a\) impulses,
composed across \(B\) regions, has the same logarithmic rate in \(B\) as
\(K_a(\theta)\).

*Proof.* A collision places two impulses in one epoch with probability at
most

\[
  \binom a2\frac{127}{L-1}.
  \tag{7}
\]

This probability is \(O_a(L^{-1})\). Outside the collision event, each
impulse occupies a distinct epoch.

A live state that is not uniform can persist for one zero-input epoch. The
next field multiplier makes it uniform on the nonzero field elements. Every
transient output tilt is \(1+O(L^{-1})\), uniformly over the fixed state
space.

The distinct-epoch process therefore has the two state classes in (2).
Equation (3) gives the transform between consecutive impulses. Averaging the
limiting spacings gives (6).

The collision and boundary corrections change each nonzero entry by relative
\(1+O_a(M/L)\). Since \(B=O(\log L)\), their total logarithmic contribution
is \(o(B)\). \(\square\)

The constants hidden by \(O_a(\cdot)\) depend on \(a\). Lemma 1 is not
uniform when \(a\) grows with \(L\).

## Matrix for \(Q\) potential outer rows

For one region, each of the \(Q\) active outer rows contributes either zero
or one bit. Summing all \(2^Q\) row-bit patterns gives

\[
  T_Q(\theta):=
  \sum_{a=0}^Q\binom Qa K_a(\theta).
  \tag{8}
\]

The following representation evaluates (8) without constructing every
\(K_a\). Give each potential row an independent dummy location in \([0,1]\).
If its bit is zero, apply the identity at that location. If its bit is one,
apply \(P\). The zero- and one-bit choices therefore sum to

\[
  H:=I+P
  =
  \begin{pmatrix}
    1&1\\
    r&2-r
  \end{pmatrix}.
  \tag{9}
\]

Ordering the \(Q\) dummy locations gives \(Q+1\) spacings. Their joint law is
Dirichlet with every parameter equal to one. Hence

\[
\begin{split}
  T_Q(\theta)
  =Q!\int_{\substack{x_i\ge0\\\sum_{i=0}^Qx_i=1}}
  D_\gamma(x_0)H D_\gamma(x_1)H\cdots
  H D_\gamma(x_Q)\,d\boldsymbol x.
  \tag{10}
\end{split}
\]

An identity insertion does not change the state or the accumulated live
time. Thus the dummy locations for zero bits do not change the transfer.
Expanding every \(H=I+P\) in (10) recovers (8).

## Coefficient and beta-integral formula

Let \(D_y:=\operatorname{diag}(1,y)\), and define nonnegative matrices
\(C_{Q,k}\) by

\[
  D_y(H D_y)^Q
  =\sum_{k=0}^{Q+1}C_{Q,k}y^k.
  \tag{11}
\]

The index \(k\) counts live spacings along a state path. The sum of any fixed
\(k\) Dirichlet spacings has the beta distribution

\[
  S_{Q,k}\sim\operatorname{Beta}(k,Q+1-k).
  \tag{12}
\]

Use the endpoint conventions \(S_{Q,0}=0\) and \(S_{Q,Q+1}=1\). Then

\[
  T_Q(\theta)
  =
  \sum_{k=0}^{Q+1}
  C_{Q,k}\,
  \mathbb E[e^{-\gamma S_{Q,k}}].
  \tag{13}
\]

For \(0<k<Q+1\),

\[
  \mathbb E[e^{-\gamma S_{Q,k}}]
  ={}_1F_1(k;Q+1;-\gamma).
  \tag{14}
\]

Equations (11)--(14) are an exact finite-dimensional evaluation of the
continuum matrix. They use no sampling.

For \(Q=1\), equation (13) equals \(K_0+K_1\). For \(Q=2\), it equals
\(K_0+2K_1+K_2\). The evaluator reproduces both earlier Perron thresholds to
less than \(10^{-8}\).

## Exact outer-row reduction

Let \(Z_{D,Q}\) count bad messages supported on exactly \(Q\) outer blocks.
The probability space is the random-outer RM2Sub ensemble from
`RANDOM_OUTER_RM2SUB_RAMP.md`.

For each fixed nonzero row word, the expected message multiplicity is

\[
  \rho_B:=\frac{2^{B/2}-1}{2^B-1}.
  \tag{15}
\]

The \(Q\) active outer maps are independent. Inclusion-exclusion removes
outer-row tuples in which at least one row is zero. In the continuum, the
resulting scalar transfer is

\[
  \sum_{j=0}^Q
  (-1)^j\binom Qj
  e_0^{\mathsf T}T_{Q-j}(\theta)^B\boldsymbol1.
  \tag{16}
\]

Every term before inclusion-exclusion is nonnegative. Equation (16) is a
scalar identity; its summands have alternating signs.

Let \(\lambda_Q(\theta)\) be the Perron root of \(T_Q(\theta)\). Adding one
potential row adds all patterns in which its bit is one. Therefore

\[
  T_Q(\theta)\ge T_{Q-1}(\theta)
  \tag{17}
\]

entrywise, with a strict inequality in at least one entry. Both matrices are
irreducible for \(\theta>0\). Hence

\[
  \lambda_Q(\theta)>\lambda_{Q-1}(\theta).
  \tag{18}
\]

The first term in (16) has the unique largest Perron rate.

## Fixed-occupation exponent

For \(\delta\in(0,1/2)\), define

\[
  \Phi_{Q,\delta}(\theta)
  :=\log_2\lambda_Q(\theta)
  +\frac{\delta\theta}{\ln2}
  \tag{19}
\]

and

\[
  \eta_Q(\delta)
  :=\frac Q2-\inf_{\theta>0}\Phi_{Q,\delta}(\theta).
  \tag{20}
\]

**Theorem 2 (every fixed occupation).** Fix \(Q\ge1\). Let \(B\to\infty\)
be even, let \(L\to\infty\) be divisible by \(128\), and assume
\(B=O(\log L)\). Then

\[
  \log_2\mathbb E[Z_{\lfloor\delta LB\rfloor,Q}]
  \le
  Q\log_2L-\eta_Q(\delta)B+o_Q(B).
  \tag{21}
\]

*Proof.* Choose the active outer blocks in \(\binom LQ\) ways. Their expected
row multiplicity is \(\rho_B^Q\), and

\[
  Q\log_2\rho_B=-\frac{QB}{2}+o_Q(1).
  \tag{22}
\]

Lemma 1 and equation (16) give the inner transfer. Equation (18) shows that
the first inclusion-exclusion term determines its exponential rate. The
Chernoff factor contributes

\[
  z_L^{-\lfloor\delta LB\rfloor}
  =2^{\delta\theta B/\ln2+o(B)}.
  \tag{23}
\]

Combining these factors and minimizing over \(\theta\) proves (21).
\(\square\)

If \(B=c\log_2N+O(1)\), Markov's inequality gives

\[
  \Pr[Z_{\lfloor\delta N\rfloor,Q}>0]
  \le
  N^{Q-c\eta_Q(\delta)+o_Q(1)}.
  \tag{24}
\]

Thus the common-tilt fixed-\(Q\) condition is

\[
  c>\frac{Q}{\eta_Q(\delta)}.
  \tag{25}
\]

Theorem 2 holds for every fixed \(Q\). It does not permit \(Q=Q_N\to\infty\),
because the collision and boundary errors in Lemma 1 depend on \(Q\).

## Diagnostic at relative distance \(0.11\)

The binary64 evaluator uses (11)--(14) and optimizes one common tilt. Selected
results are:

| \(Q\) | \(\theta/Q\) | \(\eta_Q/Q\) | threshold \(Q/\eta_Q\) |
|---:|---:|---:|---:|
| 1 | 1.3889828703 | 0.2788105891 | 3.5866643489 |
| 2 | 1.3890923390 | 0.2787616464 | 3.5872940664 |
| 4 | 1.3893275056 | 0.2786565959 | 3.5886464366 |
| 8 | 1.3898708886 | 0.2784144266 | 3.5917678987 |
| 12 | 1.3905294930 | 0.2781222009 | 3.5955418047 |
| 16 | 1.3913287393 | 0.2777698552 | 3.6001026792 |
| 32 | 1.3966480144 | 0.2755077125 | 3.6296624549 |
| 64 | 1.4301403201 | 0.2648473349 | 3.7757601013 |
| 128 | 1.5474688793 | 0.2436745158 | 4.1038349726 |
| 256 | 1.6241760918 | 0.2282568179 | 4.3810301448 |
| 512 | 1.6385259487 | 0.2194122864 | 4.5576299139 |

The \(c=18/5\) common-tilt exponent remains positive through sampled
\(Q=12\). At \(Q=16\), its binary64 value is

\[
  c\eta_{16}(0.11)-16=-0.0004563392\ldots .
  \tag{26}
\]

The sign in (26) is numerical. It shows that the existing \(c=18/5\)
certificate does not extend automatically to all fixed occupations. It does
not prove that the construction has a bad codeword.

The sampled threshold reaches \(4.55763\) at \(Q=512\). A different,
nonoptimized tilt closes the complete fixed-occupation range.
RM2SUB_UNIFORM_FIXED_OCCUPATION.md proves an exact rational \(c=9\)
certificate for every fixed \(Q\ge3\) at \(\delta=0.11\). That result is
uniform in the numerical value of a fixed \(Q\), but not for a sequence
\(Q=Q_N\to\infty\).

## Status and next obligation

- **Proved here:** the ordered-integral formula (6), the matrix formula
  (10), the coefficient formula (13), and Theorem 2 for every fixed \(Q\).
- **Numerical only:** the optimized table and the first sampled failure of
  \(c=18/5\).
- **Not implied by the diagnostic:** failure of relative distance \(0.11\),
  failure of logarithmic blocks, or a lower bound on the necessary constant.
- **Closed in the companion note:** one exact logarithmic constant for the
  numerical value of every fixed \(Q\).
- **Closed elsewhere:** `RM2SUB_DENSE_OCCUPATION.md` covers occupations
  growing with \(L\) and the complete bulk bridge.
- **Open:** the structured-outer spectrum comparison.

The completed finite-\(L\) proof uses Bernoulli conditioning and finite epoch
transfers. It does not need a uniform extension of this continuum formula.
