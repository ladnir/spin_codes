# Asymptotic distance certificate for random-outer RM2Sub SPIN

## Certified statement

For each positive integer \(m\), set \(L_m:=128m\). Let \(B_m\) be the least
positive even integer satisfying

\[
  B_m\ge9\log_2(L_mB_m),
  \tag{1}
\]

and set \(N_m:=L_mB_m\). Then

\[
  B_m=9\log_2N_m+O(1),
  \qquad
  \frac{N_{m+1}-N_m}{N_m}=o(1).
  \tag{2}
\]

At length \(N_m\), use the random-outer RM2Sub ensemble from
`RANDOM_OUTER_RM2SUB_RAMP.md`. Every local outer is an independent uniform
linear injection

\[
  \mathbb F_2^{B_m/2}\longrightarrow\mathbb F_2^{B_m}.
\]

Use independent local coordinate permutations, independent region
permutations, and independent nonzero RM2Sub multipliers. The inner uses the
fixed \((t,s)=(128,19)\) constituent.

Let \(\mathcal C_m\) be the resulting binary linear code. Then

\[
  \Pr\!\left[
    d_{\min}(\mathcal C_m)\le\lfloor0.11N_m\rfloor
  \right]=o(1).
  \tag{3}
\]

The probability in (3) is over all sampled setup objects. The code has rate
exactly one half. In particular, a deterministic setup realization with
asymptotic relative distance at least \(0.11\) exists for every sufficiently
large admissible length.

Direct multiplication by the random local generator matrices costs

\[
  O(N_mB_m)=O(N_m\log N_m)
\]

bit operations. This certificate is a distance theorem, not a linear-time
outer construction.

## Occupation partition

For a nonzero message, let \(Q\in[L]\) be the number of nonzero local message
blocks. Let \(Z_{D,Q}\) count bad messages with occupation \(Q\), where
\(D=\lfloor0.11N\rfloor\).

Every occupation sequence has a subsequence in one of three classes:

1. \(Q\) is constant;
2. \(Q\to\infty\) and \(Q/L\le10^{-4}\);
3. \(Q/L\ge10^{-4}\).

The next three sections give uniform bounds for these classes.

## Constant occupation

`RM2SUB_ONE_ACTIVE_CONTINUUM.md` and
`RM2SUB_TWO_ACTIVE_CONTINUUM.md` give sharper bounds for \(Q=1,2\).
`RM2SUB_UNIFORM_FIXED_OCCUPATION.md` proves, for every fixed \(Q\ge3\),

\[
  \Pr[Z_{D,Q}>0]
  \le
  N^{-0.04439431461547Q+o_Q(1)}.
  \tag{4}
\]

The bound in (4) uses the schedule in (1). Hence every finite set of
constant occupations contributes \(o(1)\).

## Growing vanishing occupation

Assume \(Q\to\infty\) and put \(\alpha:=Q/L\le10^{-4}\).
`RM2SUB_DENSE_OCCUPATION.md` proves the finite transfer bound and the exponent

\[
  E(\alpha)
  \le-\kappa_0\alpha,
  \qquad
  \kappa_0:=\frac{780897}{6665600}>0.11715.
  \tag{5}
\]

Thus the principal factor is at most \(\exp(-\kappa_0QB)\). The choice of the
active positions costs

\[
  \ln\binom LQ
  \le Q\ln(eL/Q)
  \le Q(1+\ln L).
  \tag{6}
\]

Equation (1) gives

\[
  \frac{1+\ln L}{B}
  \le\frac{\ln2}{9}+o(1).
  \tag{7}
\]

The strict remaining margin is at least

\[
  \kappa_0-\frac{\ln2}{9}
  >0.04013
  \tag{8}
\]

natural units per \(QB\).

For \(p=(8/5)Q/L\), Robbins' Stirling bounds give, uniformly in this range,

\[
  -\ln b_{L,Q}(p)
  \le
  L D_{\rm KL}(Q/L\Vert p)
  +\frac12\ln Q+2.
  \tag{9}
\]

The \(B\) conditioning factors in (9) therefore contribute \(o(QB)\) when
\(Q\to\infty\). The nonzero-row relaxation contributes \(o(QB)\). The fixed
Collatz vector contributes only a constant factor.

Choose a constant \(Q_0\) large enough that the uniform Stirling remainder,
the support remainder in (7), and the row-relaxation remainder consume less
than half the margin in (8). Also require

\[
  0.02Q_0\frac9{\ln2}>2.
\]

For all sufficiently large \(N\), every
\(Q_0\le Q\le10^{-4}L\) then satisfies

\[
  \mathbb E[Z_{D,Q}]
  \le \exp(-0.02QB).
  \tag{10}
\]

The finite range \(Q<Q_0\) is covered by the constant-occupation theorem.
Summing (10) over the remaining occupation levels gives \(o(1)\).

## Positive occupation

For \(10^{-4}\le\alpha\le1\), the 31-box interval certificate proves

\[
  \Psi(\alpha)\le-\varepsilon_*,
  \qquad
  \varepsilon_*
  :=4.368933504109312\times10^{-7}.
  \tag{11}
\]

The exponent uses natural units per output bit. The support factor satisfies

\[
  \frac1N\ln\binom LQ\le\frac{\ln2}{B}=o(1).
\]

The local binomial prefactor contributes \(O(B\ln L)=o(N)\), and the
nonzero-row relaxation is also \(e^{o(N)}\). Hence

\[
  \sum_{Q/L\ge10^{-4}}\mathbb E[Z_{D,Q}]
  \le \exp[-(\varepsilon_*+o(1))N].
  \tag{12}
\]

## Completion of the first moment

The fixed-occupation theorem and a finite union give

\[
  \sum_{Q=1}^{Q_0-1}\Pr[Z_{D,Q}>0]=o(1).
\]

The additional condition on \(Q_0\) makes the sum of (10) over
\(Q_0\le Q\le10^{-4}L\) tend to zero. Equation (12) handles every remaining
occupation. A union bound, together with Markov's inequality for the two
growing ranges, proves (3).

## Arbitrary lengths

The admissible lengths satisfy the vanishing relative-gap property in (2).
The next-admissible padding and fixed-coordinate puncturing wrapper in
`ARBITRARY_LENGTH_WRAPPER.md` therefore extends the asymptotic relative
distance and rate to arbitrary target lengths, subject to its stated
dimension convention. The wrapper does not change the setup probability
space.

## Scope

This theorem certifies the random rate-half outer. It does not transfer
\(0.11\) to a structured outer family. Such a transfer requires a proved
uniform spectrum or transfer-weighted comparison.

The companion linear-time audit certifies a selected Golay--BA structured
outer at distance \(0.101\) with \(B=\Theta((\ln N)^2)\). Its larger
half-Bernoulli spectrum cost prevents the present high-occupation argument
from reaching \(0.11\).
