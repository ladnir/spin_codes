# Random-block Accumulator SPIN

## Conclusion

A random block outer followed by one uniform interleaver and one accumulator
has positive relative distance with logarithmic outer blocks.  At rate one
half, every block constant \(c>2\) gives an explicit positive distance.

For a target \(\delta\in(0,1/2)\), define

\[
  \rho(\delta):=2\sqrt{\delta(1-\delta)},
  \qquad
  \lambda(\delta):=
  \frac12-\log_2(1+\rho(\delta)).
\]

If \(\lambda(\delta)>0\), then

\[
  B_N\ge c\log_2N,
  \qquad
  c\lambda(\delta)\ge1
\]

is sufficient.  For example, \(c=3\) certifies every
\(\delta\le0.0037634\), while \(c=17\) certifies every
\(\delta\le0.0330838\).  One accumulator cannot reach the rate-half
Gilbert--Varshamov distance through this proof.  Even as \(c\to\infty\), the
certified distance approaches only \(0.0449101\).

## Ensemble and probability space

Fix a rational rate \(R\in(0,1)\).  Let \(B\) be an allowed block size for
which \(RB\) is an integer, and let \(N=LB\).  Setup performs these independent
sampling operations.

1. For each \(i\in\{1,\ldots,L\}\), sample a uniform linear injection
   \[
     C_i:\mathbb F_2^{RB}\longrightarrow\mathbb F_2^B.
   \]
2. Sample \(\Pi\gets S_N\).

For \(u=(u_1,\ldots,u_N)\), define the accumulator by

\[
  \operatorname{Acc}_N(u)_t:=\sum_{j=1}^t u_j
  \qquad(t\in\{1,\ldots,N\}).
\]

For a message \(x=(x_1,\ldots,x_L)\), define

\[
  E_N(x):=
  \operatorname{Acc}_N\!\left(
    \Pi(C_1(x_1)\Vert\cdots\Vert C_L(x_L))
  \right).
\]

The sampled setup remains fixed for all messages.  Every realized encoder is
injective because all three constituent maps are injective.  Distance failure
probabilities below are over \((C_1,\ldots,C_L,\Pi)\).

## Exact accumulator tail

Fix a nonzero outer word of weight \(h\).  The permutation makes its image
uniform on the weight-\(h\) Hamming slice.  For an integer \(D\), define

\[
  p_{N,h}(D):=
  \Pr_{\Pi\gets S_N}\!\left[
    \operatorname{wt}(\operatorname{Acc}_N(\Pi(u)))\le D
  \right].
\]

The probability depends on \(u\) only through \(h\).

Set \(k:=\lceil h/2\rceil\).  The exact accumulator enumerator gives

\[
  p_{N,h}(D)
  =
  \frac1{\binom Nh}
  \sum_{w=1}^{D}
  \binom{w-1}{k-1}
  \binom{N-w}{h-k}.
  \tag{1}
\]

**Lemma 1 (hypergeometric representation).** Let \(X\) count the marked
positions in a size-\(h\) sample drawn without replacement from a population
of \(N\) positions, of which the first \(D\) are marked.  Then

\[
  p_{N,h}(D)=\Pr[X\ge k].
  \tag{2}
\]

*Proof.* Order the \(h\) sampled positions.  The probability that the
\(k\)-th order statistic equals \(w\) is

\[
  \frac{\binom{w-1}{k-1}\binom{N-w}{h-k}}{\binom Nh}.
\]

Summing over \(w\le D\) gives (1).  The \(k\)-th order statistic is at most
\(D\) exactly when at least \(k\) sampled positions are marked. \(\square\)

**Lemma 2 (sharp uniform contraction).** Fix \(\delta\in(0,1/2)\) and set
\(D:=\lfloor\delta N\rfloor\).  For every \(1\le h\le N\),

\[
  p_{N,h}(D)
  \le \rho(\delta)^h,
  \qquad
  \rho(\delta):=2\sqrt{\delta(1-\delta)}.
  \tag{3}
\]

*Proof.* Let \(\bar\delta:=D/N\le\delta\).  For \(t>0\), Maclaurin's
inequality for elementary symmetric means gives the standard
without-replacement moment bound

\[
  \mathbb E[e^{tX}]
  \le(1-\bar\delta+\bar\delta e^t)^h
  \le(1-\delta+\delta e^t)^h.
\]

The second inequality holds because \(e^t>1\).  Lemma 1 and
\(k\ge h/2\) imply

\[
  p_{N,h}(D)
  \le
  \inf_{t>0}
  e^{-th/2}(1-\delta+\delta e^t)^h.
\]

The minimizer satisfies \(e^t=(1-\delta)/\delta\).  Substitution gives (3).
\(\square\)

Equation (3) is uniform in \(N\) and \(h\).  It improves the earlier bound
\((4e\delta)^{h/2}\) to \((4\delta(1-\delta))^{h/2}\).

## Random-block outer law

For one uniform rate-\(R\) injection, define

\[
  G_{B,R}(z):=
  \mathbb E_C\!\left[
    \sum_{v\in\mathbb F_2^{RB}\setminus\{0\}}
    z^{\operatorname{wt}(C(v))}
  \right].
\]

The image of each nonzero \(v\) is uniform on
\(\mathbb F_2^B\setminus\{0\}\).  Therefore

\[
  G_{B,R}(z)
  =
  (2^{RB}-1)\frac{(1+z)^B-1}{2^B-1}.
  \tag{4}
\]

Independence across blocks gives

\[
  \mathbb E\!\left[
    1+\sum_{h=1}^N A_{N,h}^{\mathrm{out}}z^h
  \right]
  =(1+G_{B,R}(z))^L.
  \tag{5}
\]

For \(z>0\), equation (4) implies

\[
  G_{B,R}(z)
  \le2q_R(z)^B,
  \qquad
  q_R(z):=2^{R-1}(1+z).
  \tag{6}
\]

## Asymptotic theorem

Define the sparse exponent

\[
  \lambda_A(R,\delta)
  :=-log_2q_R(\rho(\delta))
  =1-R-\log_2(1+2\sqrt{\delta(1-\delta)}).
  \tag{7}
\]

Fix \(c>0\).  For every sufficiently large integer \(L\), let \(B_L\) be
the least allowed block size satisfying

\[
  B_L\ge c\log_2(LB_L).
  \tag{8}
\]

Set \(N_L:=LB_L\) and \(K_L:=RN_L\).

**Theorem 1 (random-block Accumulator SPIN).** Suppose

\[
  \lambda_A(R,\delta)>0,
  \qquad
  c\lambda_A(R,\delta)\ge1.
  \tag{9}
\]

For the ensemble above,

\[
  \Pr\!\left[d_{\min}(E_{N_L})\le\delta N_L\right]
  =o(1)
  \qquad(L\to\infty).
  \tag{10}
\]

The admissible family has rate \(R\), block size \(\Theta(\log N_L)\), and
positive relative distance.  Dense local matrices give
\(\Theta(N_L\log N_L)\) encoding cost.  The accumulator itself costs
\(O(N_L)\) bit operations.

*Proof.* Let \(Z_D\) count nonzero messages whose encoded word has weight at
most \(D:=\lfloor\delta N_L\rfloor\).  The first-moment identity and Lemma 2
give

\[
  \mathbb E[Z_D]
  \le
  \sum_{h=1}^{N_L}
  \mathbb E[A_{N_L,h}^{\mathrm{out}}]\rho(\delta)^h
  =
  (1+G_{B_L,R}(\rho(\delta)))^L-1.
\]

Set \(\lambda:=\lambda_A(R,\delta)\).  Equations (6), (8), and (9) give

\[
  L G_{B_L,R}(\rho(\delta))
  \le
  \frac{2N_L}{B_L}2^{-\lambda B_L}
  \le
  \frac{2}{B_L}N_L^{1-c\lambda}
  =o(1).
\]

Hence \(\mathbb E[Z_D]=o(1)\).  Markov's inequality proves (10).
\(\square\)

The least-allowed rule gives \(B_L=\Theta(\log N_L)\) and
\(N_{L+1}-N_L=o(N_L)\).  The wrapper in ARBITRARY_LENGTH_WRAPPER.md therefore
extends Theorem 1 to every sufficiently large requested length with \(o(1)\)
relative rate and distance loss.

## Rate-half constants

At \(R=1/2\), condition (9) becomes

\[
  c\left[
    \frac12-
    \log_2(1+2\sqrt{\delta(1-\delta)})
  \right]
  \ge1.
  \tag{11}
\]

For fixed \(c>2\), define

\[
  r_c:=2^{1/2-1/c}-1,
  \qquad
  \delta_c:=\frac{1-\sqrt{1-r_c^2}}2.
  \tag{12}
\]

Then every \(0<\delta\le\delta_c\) satisfies (11).  No positive \(\delta\)
is obtained when \(c\le2\).

| \(c\) in \(B_N\ge c\log_2N\) | largest certified \(\delta_c\) |
|---:|---:|
| 3 | 0.00376340 |
| 4 | 0.00903140 |
| 5 | 0.01354027 |
| 6 | 0.01718506 |
| 8 | 0.02253632 |
| 10 | 0.02620820 |
| 16 | 0.03242567 |
| 17 | 0.03308384 |
| 32 | 0.03831472 |

As \(c\to\infty\), the distance ceiling is

\[
  \delta_\infty
  =
  \frac{1-\sqrt{1-(\sqrt2-1)^2}}2
  =0.0449101394\ldots.
  \tag{13}
\]

The closed formulas (11)--(13) are proved consequences of the ensemble.
The decimal table is a numerical rendering of (12), not a finite-length
certificate.

## Status

- **Proved:** the hypergeometric identity (2) and uniform contraction (3).
- **Proved:** the exact random-block law (4)--(6).
- **Proved:** the asymptotic theorem (10), including the equality case
  \(c\lambda_A=1\).
- **Proved:** the rate-half tradeoff (11)--(13) and the arbitrary-length
  extension.
- **Open optimization:** a stronger inner than one accumulator is required to
  move the distance substantially beyond \(0.04491\) at rate one half.
- **Not claimed:** an outward-rounded finite certificate for a particular
  sampled setup.
