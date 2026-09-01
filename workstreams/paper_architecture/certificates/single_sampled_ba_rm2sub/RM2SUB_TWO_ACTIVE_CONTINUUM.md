# Two-active continuum for the fixed RM2Sub-S19 inner

## Result

Use the random rate-half outer ensemble and fixed \(t=128,s=19\) inner from
`RANDOM_OUTER_RM2SUB_RAMP.md`. Let \(Z_{D,2}\) count messages supported on
exactly two outer blocks whose output weight is at most \(D\).

For \(D=\lfloor0.11002N\rfloor\), \(B=c\log_2N+O(1)\), and
\(L\to\infty\), the optimized two-active continuum gives

\[
  c>3.587809926\ldots .
  \tag{1}
\]

The optimized occupation-two threshold is \(0.00062994\ldots\) larger than
the occupation-one threshold. Thus occupation two is the more demanding
class at the optimized decimal level.

The decimal in (1) uses binary64 optimization. An exact rational inequality
certifies the same convenient schedule as occupation one:

\[
  B\ge\frac{18}{5}\log_2N+O(1).
  \tag{2}
\]

Under (2),

\[
  \Pr[Z_{\lfloor0.11002N\rfloor,2}>0]
  \le N^{-0.0049419\ldots+o(1)}.
  \tag{3}
\]

Together with the occupation-one certificate, (3) excludes bad messages
supported on one or two outer blocks. It does not control occupation
\(q\ge3\).

## Exact two-row occupancy law

Fix two distinct outer-block indices. For each active block, the random outer
map produces an independent nonzero row in \(\mathbb F_2^B\). Let

\[
  v,w\in\mathbb F_2^B\setminus\{0\}
  \tag{4}
\]

denote those rows. Region \(j\) contains \(v_j+w_j\) active input bits, where
the sum in this sentence is over the integers. Thus the region occupancy is
zero, one, or two.

The region permutation maps the two source rows to two distinct positions.
Conditioned on \(v_j=w_j=1\), the two positions form a uniform two-element
subset of \([L]\). They never cancel.

Let \(E=L/t\). Use the exact full-state epoch kernels \(W_0(z)\) and
\(W_1(z)\) from `RM2SUB_ONE_ACTIVE_CONTINUUM.md`. Define \(W_2(z)\) by
averaging the exact epoch kernel over a uniform two-element coordinate subset
of \([t]\).

The exact two-impulse region transfer is

\[
\begin{split}
  \mathcal R_{2,L}(z)
  :=\frac1{\binom L2}\Bigg(
  &t^2
  \sum_{0\le r<u<E}
  W_0^rW_1W_0^{u-r-1}W_1W_0^{E-1-u}\\
  &+\binom t2
  \sum_{r=0}^{E-1}
  W_0^rW_2W_0^{E-1-r}
  \Bigg).
  \tag{5}
\end{split}
\]

The first line of (5) covers impulses in distinct epochs. The second line
covers impulses in the same epoch. The counting identity

\[
  t^2\binom E2+E\binom t2=\binom L2
  \tag{6}
\]

shows that (5) averages every two-position set exactly once.

The exact same-epoch probability is

\[
  \frac{E\binom t2}{\binom L2}
  =\frac{t-1}{L-1}.
  \tag{7}
\]

This probability is \(O(L^{-1})\). It must remain in the finite transfer, but
it vanishes in the continuum.

## Exact two-active first moment

Let

\[
  \rho_{B,B/2}:=\frac{2^{B/2}-1}{2^B-1}.
  \tag{8}
\]

For a fixed pair of active outer blocks, every ordered pair
\((v,w)\) in (4) has expected message multiplicity
\(\rho_{B,B/2}^2\).

Let \(\mathcal R_{0,L}\) and \(\mathcal R_{1,L}\) be the exact zero- and
one-impulse region transfers. Define

\[
\begin{split}
  \mathcal S_{2,B,L}(z)
  :=\sum_{\substack{v,w\in\mathbb F_2^B\\v\ne0,\ w\ne0}}
  e_0^{\mathsf T}
  \prod_{j=1}^B
  \mathcal R_{v_j+w_j,L}(z)
  \boldsymbol1.
  \tag{9}
\end{split}
\]

The product follows the region order. Each summand in (9) is nonnegative.
Inclusion-exclusion gives the equivalent scalar identity

\[
\begin{split}
  \mathcal S_{2,B,L}(z)
  =e_0^{\mathsf T}\Big(
  &(\mathcal R_{0,L}+2\mathcal R_{1,L}+\mathcal R_{2,L})^B\\
  &-2(\mathcal R_{0,L}+\mathcal R_{1,L})^B
  +\mathcal R_{0,L}^B
  \Big)\boldsymbol1.
  \tag{10}
\end{split}
\]

Equation (10) is an identity for the complete scalar expression. It does not
assert that the matrix inside parentheses is entrywise nonnegative.

**Lemma 1 (finite two-active reduction).** The probability space includes the
random outer maps, local coordinate permutations, region permutations, and
field multipliers. For every \(z\in(0,1)\),

\[
  \mathbb E[Z_{D,2}]
  \le
  \binom L2
  \rho_{B,B/2}^2
  z^{-D}
  \mathcal S_{2,B,L}(z).
  \tag{11}
\]

*Proof.* Choose the two active outer blocks in \(\binom L2\) ways. Sum over
their two nonzero outer rows. Independence of the outer maps gives the
multiplicity \(\rho_{B,B/2}^2\).

For every fixed pair \((v,w)\), the independent region permutations give the
transfers in (9). The Chernoff inequality

\[
  \boldsymbol1\{y\le D\}\le z^{-D}z^y
  \tag{12}
\]

then gives (11). \(\square\)

## Two-impulse long-region matrix

Fix \(\theta>0\), and set \(z_L=e^{-\theta/L}\). Define

\[
  M:=2^{19}-1,
  \qquad
  r:=\frac1M,
  \qquad
  q:=1-r,
  \qquad
  p:=\frac{2^{18}}M.
  \tag{13}
\]

Here \(r\) is the live-state termination probability at an active epoch.
The symbol \(q\) in (13) denotes survival, not outer occupation.

Set

\[
  \gamma:=p\theta,
  \qquad
  a:=e^{-\gamma},
  \qquad
  f:=\frac{1-e^{-\gamma}}{\gamma}.
  \tag{14}
\]

Let \(U<V\) be the ordered locations of two independent uniform points in
\([0,1]\). Define

\[
  A(\gamma)
  :=\mathbb E[e^{-\gamma(V-U)}]
  =\frac{2(\gamma-1+e^{-\gamma})}{\gamma^2}
  \tag{15}
\]

and

\[
\begin{split}
  H(\gamma)
  &:=\mathbb E[e^{-\gamma V}]
   =\mathbb E[e^{-\gamma(1-U)}]\\
  &=\mathbb E[e^{-\gamma(U+1-V)}]
   =\frac{2(1-(1+\gamma)e^{-\gamma})}{\gamma^2}.
  \tag{16}
\end{split}
\]

The zero and one-impulse limiting matrices are

\[
  K_0=
  \begin{pmatrix}
    1&0\\
    0&a
  \end{pmatrix},
  \qquad
  K_1=
  \begin{pmatrix}
    0&f\\
    rf&qa
  \end{pmatrix}.
  \tag{17}
\]

The two-impulse limiting matrix is

\[
  K_2=
  \begin{pmatrix}
    rA&qH\\
    qrH&rH+q^2a
  \end{pmatrix}.
  \tag{18}
\]

**Lemma 2 (two-impulse region limit).** Suppose \(128\mid L\) and
\(B=O(\log L)\). The exact transfer \(\mathcal R_{2,L}(z_L)\), composed with
the zero- and one-impulse transfers across \(B\) regions, has the same
logarithmic rate in \(B\) as \(K_2\).

*Proof.* Equation (7) shows that same-epoch collisions have mass
\(O(L^{-1})\). Their output tilt is between
\(e^{-2\theta t/L}\) and one. Their state transition is also bounded because
the fixed state space has \(2^{19}\) elements. Thus these collisions change a
nonzero limiting entry by relative \(1+O(M/L)\).

It remains to analyze distinct epochs. Start from state zero. The first
impulse activates the state. If the second impulse terminates it, the state
is live only between \(U\) and \(V\). This branch gives \(rA\). If the second
impulse does not terminate it, the state is live from \(U\) through the end
of the region. This branch gives \(qH\).

Now start from a uniform nonzero state. If both impulses survive, the state
is live throughout the region. This branch gives \(q^2a\). If the first
impulse survives and the second terminates, the state is live through \(V\).
This branch gives \(qrH\) and ends in state zero.

If the first impulse terminates, the second impulse reactivates the state.
The live intervals are the prefix before \(U\) and the suffix after \(V\).
This branch gives \(rH\). These cases prove (18).

An active epoch can leave a nonuniform live state at a region boundary. A
following zero-input epoch refreshes that state. The boundary event has mass
\(O(E^{-1})\). Since \(M\) is fixed and \(B=o(E)\), all collision and
boundary errors contribute \(o(B)\) to the logarithm. \(\square\)

The distinct columns of the audited state-update map ensure that two
different coordinates have a nonzero syndrome. Therefore a same-epoch
weight-two input also activates a zero state. This finite fact prevents an
unmodeled zero-to-zero branch in the collision term.

## Spectral exponent

Define

\[
  T_2(\theta):=K_0(\theta)+2K_1(\theta)+K_2(\theta)
  \tag{19}
\]

and let \(\lambda_{2,s}(\theta)\) be its Perron root. Define

\[
  \Phi_{2,s,\delta}(\theta)
  :=\log_2\lambda_{2,s}(\theta)
  +\frac{\delta\theta}{\ln2}
  \tag{20}
\]

and

\[
  \eta_{2,s}(\delta)
  :=1-\inf_{\theta>0}\Phi_{2,s,\delta}(\theta).
  \tag{21}
\]

**Theorem 3 (occupation-two exponent).** Let \(B\to\infty\) be even, let
\(L\to\infty\) be divisible by \(128\), and assume \(B=O(\log L)\). Then

\[
  \log_2\mathbb E[Z_{\lfloor\delta LB\rfloor,2}]
  \le
  2\log_2L-\eta_{2,19}(\delta)B+o(B).
  \tag{22}
\]

*Proof.* Apply Lemma 1 with \(z=z_L\). The two random outer factors satisfy

\[
  2\log_2\rho_{B,B/2}=-B+o(1).
  \tag{23}
\]

Lemmas 1 and 2 give the leading transfer \(T_2(\theta)^B\). The two
inclusion-exclusion terms in (10) have Perron rate
\(\lambda_{1,s}(\theta)^B\). Entrywise,
\(T_2>K_0+K_1\) on a nonzero entry. Perron--Frobenius monotonicity therefore
gives

\[
  \lambda_{2,s}(\theta)>\lambda_{1,s}(\theta).
  \tag{24}
\]

Hence the subtracted terms do not change the exponential rate. The Chernoff
factor contributes

\[
  z_L^{-D}=2^{\delta\theta B/\ln2+o(B)}.
  \tag{25}
\]

Combining (23)--(25), then minimizing over \(\theta\), proves (22).
\(\square\)

For \(B=c\log_2N+O(1)\), Markov's inequality gives

\[
  \Pr[Z_{\lfloor\delta N\rfloor,2}>0]
  \le
  N^{2-c\eta_{2,19}(\delta)+o(1)}.
  \tag{26}
\]

Thus occupation two closes when

\[
  c>\frac{2}{\eta_{2,19}(\delta)}.
  \tag{27}
\]

## Constants

Binary64 optimization gives:

| \(\delta\) | \(\eta_{2,19}(\delta)\) | threshold \(2/\eta_{2,19}(\delta)\) |
|---:|---:|---:|
| 0.10000 | 0.5976083809 | 3.3466732793 |
| 0.10500 | 0.5775647292 | 3.4628153331 |
| 0.10850 | 0.5635354978 | 3.5490222136 |
| 0.11000 | 0.5575232928 | 3.5872940664 |
| 0.11002 | 0.5574431313 | 3.5878099265 |

At \(\delta=0.11002\), the numerical optimizer gives

\[
  \theta=2.778183527\ldots,
  \qquad
  \lambda_{2,19}(\theta)=1.001101901\ldots .
  \tag{28}
\]

For the exact certificate, choose

\[
  \theta_0:=\frac{2\ln2}{p}.
  \tag{29}
\]

Then \(a=1/4\). The formulas for \(f,A,H\) become expressions in
\(\ln2\). The receipt bounds \(\ln2\) with the positive atanh series, bounds
the Perron square root with exact integers, and uses
\(\ln(1+x)\le x\). It certifies

\[
  \eta_{2,19}(0.11002)>0.5569283272
  \tag{30}
\]

and

\[
  \frac{18}{5}\eta_{2,19}(0.11002)-2
  >0.0049419781.
  \tag{31}
\]

Equations (30)--(31) prove (2)--(3). The optimized decimals are not needed
for that conclusion.

## Status and next obligation

- **Proved here:** the exact two-position region law, the collision
  probability, the limiting matrix \(K_2\), the occupation-two exponent, and
  the rational \(c=18/5\) certificate.
- **Imported finite facts:** the fixed RM2Sub parameters, nonzero coordinate
  forms of \(A\), and distinct nonzero columns of the state-update map.
- **Numerical only:** the optimized threshold in (1) and the decimal table.
- **Subsequent result:** RM2SUB_FIXED_OCCUPATION_CONTINUUM.md proves the
  common-tilt theorem for every fixed occupation.
- **Open:** a bound uniform as occupation grows, the bulk bridge, and the
  structured-outer family.

The next useful task is a uniform large-occupation bound for the coefficient
formula in RM2SUB_FIXED_OCCUPATION_CONTINUUM.md.
