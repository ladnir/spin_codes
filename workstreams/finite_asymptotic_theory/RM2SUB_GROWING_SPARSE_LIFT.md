# Finite-region lift for growing RM2Sub occupation

## Supersession note

This note records a valid conditional cluster-comparison route. The complete
random-outer theorem no longer depends on its open weighted-norm lemma.
`RM2SUB_DENSE_OCCUPATION.md` uses a Bernoulli conditioning argument and a
finite four-state transfer to cover every growing sparse occupation directly.

## Objective

The fixed-occupation exponent now has a constant independent of the numerical
value of \(Q\). A growing-occupation theorem additionally needs a continuum
error that is uniform when \(Q=Q_N\to\infty\).

This note separates the exact finite reductions from one remaining
construction-specific comparison. The first target range is

\[
  \frac{Q_N^2B_N}{L_N}\longrightarrow0.
  \tag{1}
\]

Condition (1) is conservative. It makes all pair-collision and discretization
errors vanish after composition across \(B_N\) regions.

## Exact outer relaxation

Use the random-outer RM2Sub-S19 ensemble. Fix a support of \(Q\) active outer
blocks. Each active message block maps to an independent uniform element of
\(\mathbb F_2^B\setminus\{0\}\).

Let \(\mathsf P_{\mathrm{nz}}\) denote this nonzero-row law, and let
\(\mathsf P_{\mathrm{bin}}\) denote the uniform law on \(\mathbb F_2^B\).
For every nonnegative function \(F\) of the \(Q\) rows,

\[
  \mathbb E_{\mathsf P_{\mathrm{nz}}^{\,Q}}[F]
  \le
  (1-2^{-B})^{-Q}
  \mathbb E_{\mathsf P_{\mathrm{bin}}^{\,Q}}[F].
  \tag{2}
\]

Under the relaxed law, all \(QB\) row bits are independent and uniform. In
particular, the \(B\) region transfers are independent before the RM2Sub
state composition.

Let \(\widehat T_{Q,L}(z)\) be the exact finite-region transfer after
averaging the \(Q\) relaxed row bits, the uniform region permutation, and the
epoch multipliers. Its state interface must retain enough boundary state to
make composition exact. If \(\widehat e_0\) is the zero initial state, then

\[
\begin{split}
  \mathbb E[Z_{D,Q}]
  \le{}&
  \binom LQ(2^{B/2}-1)^Q(1-2^{-B})^{-Q}z^{-D}\\
  &{}\cdot
  \widehat e_0^{\mathsf T}
  \widehat T_{Q,L}(z)^B\boldsymbol1.
  \tag{3}
\end{split}
\]

Equation (3) is an exact relaxation. It avoids inclusion-exclusion and
remains valid when \(Q\) grows.

## Exact placement facts

In one region, the common region permutation maps the \(Q\) active source
rows to a uniform ordered sample without replacement from \([L]\). Distinct
source rows therefore occupy distinct bit positions.

For any fixed pair of source rows, the probability that both positions lie
in one epoch is

\[
  \frac{127}{L-1}.
  \tag{4}
\]

The state produced by an impulse becomes uniform live after one intervening
zero-input epoch. A pair can violate this refresh condition only if its epoch
indices are equal or adjacent. For \(t=128\), a union bound gives

\[
  \Pr[\text{a fixed pair has no intervening zero epoch}]
  \le
  \frac{3t-1}{L-1}
  =\frac{383}{L-1}.
  \tag{5}
\]

Hence the expected number of unrefreshed pairs in one region is at most

\[
  \binom Q2\frac{383}{L-1}.
  \tag{6}
\]

A boundary transient occurs only when an impulse lies in the final epoch.
Its expected number is \(tQ/L\).

These estimates identify the scale \(Q^2/L\). They do not by themselves
bound the low-output transfer on the exceptional placements.

## Uniform zero-epoch estimate

Fix \(\tau>0\), set

\[
  \theta_Q:=\frac{Q\tau}{p},
  \qquad
  z_{Q,L}:=e^{-\theta_Q/L},
  \tag{7}
\]

and let \(W:=\operatorname{wt}(A(U))\) for
\(U\gets\mathbb F_{2^{19}}^*\). Then

\[
  0\le W\le t,
  \qquad
  \mathbb E[W]=tp.
  \tag{8}
\]

Hoeffding's lemma gives

\[
\begin{split}
  \log\mathbb E[z_{Q,L}^{W}]
  &\le
  -\frac{\theta_Qtp}{L}
  +\frac{\theta_Q^2t^2}{8L^2}\\
  &=
  -\frac{Q\tau t}{L}
  +\frac{Q^2\tau^2t^2}{8p^2L^2}.
  \tag{9}
\end{split}
\]

A region has at most \(L/t\) zero-input epochs. Thus the cumulative
zero-epoch correction in one region is at most

\[
  \exp\!\left(
    \frac{Q^2\tau^2t}{8p^2L}
  \right).
  \tag{10}
\]

Equation (10) is uniform in \(Q\). Across \(B\) regions, its logarithmic cost
is \(O(BQ^2/L)\).

## The remaining finite-kernel comparison

Let

\[
  T_Q(\theta_Q)
\]

be the continuum reset matrix from
RM2SUB_FIXED_OCCUPATION_CONTINUUM.md. Under the relaxed bit law, its
normalized version is \(2^{-Q}T_Q(\theta_Q)\).

Let \(v_*>0\) be a right Perron vector of \(H D_{1/2}\). Lift it to the full
state space by assigning its zero component to state zero and its live
component to every nonzero state. For a nonnegative full-state matrix \(A\),
define

\[
  \lVert A\rVert_{v_*}
  :=
  \max_i\frac{(Av_*)_i}{(v_*)_i}.
  \tag{11}
\]

This weighted row norm is submultiplicative. It therefore controls all \(B\)
regions without a matrix-conditioning factor.

The needed statement is the following.

**Finite-lift lemma (open).** There are constants \(C<\infty\) and
\(\alpha>0\), depending only on the frozen \(t=128,s=19\) constituent and
the selected \(\tau\), such that, for \(Q\le\alpha L\),

\[
  \left\lVert\widehat T_{Q,L}(z_{Q,L})\right\rVert_{v_*}
  \le
  2^{-Q}
  \exp\!\left(C\frac{Q^2}{L}\right)
  \left(
    e^{G_Q(1/2,\tau)}\rho(H D_{1/2})
  \right)^Q.
  \tag{12}
\]

The probability estimates (4)--(6) are insufficient to prove (12), because
discarding exceptional placements would remove their Chernoff weight. A
proof must compare their exact full-state kernels with reset paths.

The comparison can use three finite facts.

1. A separated impulse acts exactly as the reset matrix \(P\) after one
   zero-input refresh epoch.
2. Every collision cluster contains at most \(t\) impulses per epoch, and
   the state space has fixed size \(2^{19}\).
3. The reset transfer for two or more impulses has positive entries, so every
   finite collision kernel has a finite entrywise comparison constant.

The missing step is an exponential-moment bound for the product of these
cluster comparison constants. Its target cost is
\(\exp(O(Q^2/L))\), not a constant raised to \(Q\).

## Consequence of the finite-lift lemma

This subsection is conditional on (12). Choose the exact certified values

\[
  \delta:=\frac{11}{100},
  \qquad
  \tau:=\frac{11}{10},
  \qquad
  B=9\log_2N+O(1).
  \tag{13}
\]

For every \(Q\ge3\), the exact rational certificate proves

\[
  \frac{\eta_Q(\delta)}Q
  >g_*,
  \qquad
  g_*:=0.11604381273505275.
  \tag{14}
\]

Apply (12) in (3), use
\(\binom LQ\le(eL/Q)^Q\), and absorb the nonzero-row conditioning factor.
The resulting logarithmic bound is

\[
\begin{split}
  \ln\mathbb E[Z_{\lfloor\delta N\rfloor,Q}]
  \le{}&
  Q\ln(eL/Q)-g_*BQ\ln2\\
  &{}+O(BQ^2/L)+o(BQ).
  \tag{15}
\end{split}
\]

If \(Q=Q_N\to\infty\) satisfies (1), then

\[
  \ln\mathbb E[Z_{\lfloor\delta N\rfloor,Q_N}]
  \le
  -0.04439431461547\,Q_N\ln N
  -Q_N\ln Q_N
  +O(Q_N).
  \tag{16}
\]

Thus (12) would close every occupation sequence in the range (1). The result
would still leave a gap between \(Q\asymp\sqrt{L/B}\) and the all-active
endpoint.

## Status

- **Proved:** the outer relaxation (2)--(3), the exact placement estimates
  (4)--(6), the zero-epoch estimate (9)--(10), and the conditional
  implication (12) \(\Rightarrow\) (15)--(16).
- **Open:** the finite-lift lemma (12), extension beyond the first sparse
  range, and the bulk bridge.
- **Not used:** a sampled numerical grid or an assumption that exceptional
  placements can be discarded.
