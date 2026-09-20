# Random SPIN: exact outer law and proof audit

## Conclusion

Random SPIN is the cleaner first asymptotic target. Independent random block
injections have an exact expected weight enumerator. Blocks of size
\(B_N=\Theta(\log N)\) therefore provide the required sparse suppression.

The repository also contains most of a random-convolution argument. It cannot
be imported as a theorem without repair. The compiled argument concerns a
sliding random outer, and two steps in its inner proof are invalid as written.
The older Expand--Convolute theorem proves a related sparse-expander ensemble,
not Random SPIN.

The exact transfer matrix below repairs the needed convolution analysis. It
gives a uniform sparse bound and a closed linear exponent. The first complete
schedule obtained from the argument used

\[
  B_N\ge42\log_2N,
  \qquad
  m_N=\lceil2\log_2N\rceil,
\]

the proof-model ensemble has rate \(1/2\), relative distance greater than
\(0.09\), and setup failure probability \(o(1)\). The constant 42 is a
conservative proof constant, not an optimized design recommendation.

Optimizing the same proof, without changing the ensemble, gives the stronger
schedule

\[
  B_N\ge17\log_2N,
  \qquad
  m_N=\left\lceil\frac{51}{50}\log_2N\right\rceil.
\]

For this schedule, Theorem 3 proves relative distance greater than
\(0.11002\) with setup failure probability \(o(1)\). This distance is within
\(7.9\times10^{-6}\) of the rate-half Gilbert--Varshamov distance. Both
statements concern the proof-model ensemble defined next; neither is yet a
finite-length certificate for a particular sampled code.

## Proof-model ensemble

Fix an even block size \(B\) dividing \(N\), and define

\[
  L:=N/B,
  \qquad
  K:=N/2.
\]

Setup performs the following independent sampling operations.

1. For each \(i\in[L]\), sample a linear injection
   \[
     C_i\gets
     \operatorname{Inj}\!\left(
       \mathbb F_2^{B/2},\mathbb F_2^B
     \right).
   \]
2. Sample one permutation \(\Pi\gets S_N\).
3. Set \(m:=\lceil\gamma\log_2N\rceil\). For every \(t\in[N]\), sample
   \(\alpha_t\gets\mathbb F_2^m\).

The proof-model inner is the non-wrapping convolution

\[
  y_t
  :=
  u_t+
  \left\langle
    \alpha_t,(y_{t-1},\ldots,y_{t-m})
  \right\rangle,
  \tag{I}
\]

where \(y_j:=0\) for \(j\le0\). Define \(I_N(u):=y\). For every convolution
setup, \(I_N\) is an invertible linear map because each \(y_t\) contains \(u_t\)
with coefficient one.

For a message \(x=(x_1,\ldots,x_L)\), define

\[
  O_N(x):=C_1(x_1)\Vert\cdots\Vert C_L(x_L)
\]

and

\[
  E_N(x):=I_N(\Pi(O_N(x))).
\]

The complete setup remains fixed for all messages. Every realized outer map
is injective. The probability of distance failure is over the outer
injections, the global permutation, and the convolution setup.

This ensemble is a proof model. Dense random \(B/2\)-by-\(B\) local maps cost
\(\Theta(NB)=\Theta(N\log N)\) bit operations when \(B=\Theta(\log N)\).

## Exact expected outer enumerator

For \(z\in(0,1)\), define the expected nonzero local generating mass

\[
  G_B(z):=
  \mathbb E_C\!\left[
    \sum_{u\in\mathbb F_2^{B/2}\setminus\{0\}}
    z^{\operatorname{wt}(C(u))}
  \right].
\]

**Lemma 1 (exact local law).** For a uniform linear injection \(C\),

\[
  G_B(z)
  =
  (2^{B/2}-1)
  \frac{(1+z)^B-1}{2^B-1}.
  \tag{1}
\]

*Proof.* Fix a nonzero \(u\). The image \(C(u)\) is uniform on
\(\mathbb F_2^B\setminus\{0\}\). Hence

\[
  \mathbb E_C[z^{\operatorname{wt}(C(u))}]
  =
  \frac{(1+z)^B-1}{2^B-1}.
\]

Summing over the \(2^{B/2}-1\) nonzero inputs proves (1). \(\square\)

Let \(A_{N,h}^{\mathrm{out}}\) denote the number of nonzero outer words of
weight \(h\). Independence across blocks gives the exact global identity

\[
  \mathbb E\!\left[
    1+\sum_{h=1}^N A_{N,h}^{\mathrm{out}}z^h
  \right]
  =
  (1+G_B(z))^L.
  \tag{2}
\]

No independence among codewords inside one block is needed. Equation (2)
uses only linearity of expectation inside a block and independence across
blocks.

## Logarithmic blocks suffice on the outer side

Fix \(z<\sqrt2-1\), and define

\[
  q(z):=\frac{1+z}{\sqrt2}<1.
\]

Equation (1) gives

\[
  G_B(z)\le 2q(z)^B.
  \tag{3}
\]

Choose an allowed block size

\[
  B_N=\lceil c\log_2N\rceil_{\mathcal B}
\]

with

\[
  c>\frac{1}{\log_2(1/q(z))}.
  \tag{4}
\]

Then

\[
  \frac{N}{B_N}G_{B_N}(z)=o(1).
\]

Equations (2) and (3) therefore imply

\[
  \mathbb E\!\left[
    \sum_{h=1}^N A_{N,h}^{\mathrm{out}}z^h
  \right]
  =
  (1+G_{B_N}(z))^{N/B_N}-1
  =o(1).
  \tag{5}
\]

This is a proved exponential sparse bound. It is qualitatively different
from the polynomial boundary of a repeated fixed-constituent BA outer.

## Linear-weight outer exponent

For fixed \(z>0\), (1) and (2) give

\[
  \frac1N\log_2
  \mathbb E\!\left[
    1+\sum_h A_{N,h}^{\mathrm{out}}z^h
  \right]
  \le
  \max\!\left\{
    0,\log_2(1+z)-\frac12
  \right\}
  +o(1).
  \tag{6}
\]

Coefficient extraction yields the exponent

\[
  \Phi_{\mathrm{blk}}(\eta)
  :=
  \inf_{z>0}
  \left[
    \max\!\left\{
      0,\log_2(1+z)-\frac12
    \right\}
    -\eta\log_2z
  \right].
  \tag{7}
\]

The minimization gives

\[
  \Phi_{\mathrm{blk}}(\eta)
  =
  \begin{cases}
    \eta\log_2(1+\sqrt2),
      &0<\eta\le1-2^{-1/2},\\[2mm]
    H_2(\eta)-\frac12,
      &1-2^{-1/2}\le\eta<1.
  \end{cases}
  \tag{8}
\]

The second branch is the rate-half random-code exponent. Formula (8) is the
same envelope currently used for the sliding dense outer, but the random-block
derivation is exact and has no boundary-termination issue.

## What the old repository proof establishes

The compiled dense source uses:

- a systematic sliding random outer with logarithmic memory;
- one uniform interleaver; and
- a scalar time-varying recursive inner with logarithmic memory.

The files framework.tex, outerDense.tex, innerDenseScalar.tex, and
integrationDense.tex contain the reduction, outer envelope, inner envelopes,
and integration argument.

These sources contain useful proved ingredients:

1. the exact first-moment reduction through Hamming slices;
2. the exact run-count distribution of a uniform weight-\(h\) word;
3. fair-output domination while the recursive state is noncritical;
4. a linear-weight run lower-tail exponent; and
5. a binomial output-tail exponent.

The advertised \(0.109\) theorem is not currently complete. Its
moderate-weight step invokes the constant-factor geometric envelope described
below. Its linear-weight step also uses a sampled-grid maximization rather
than an interval or analytic certificate.

## Two proof defects that must not be imported

### Reversed monotonicity in the geometric envelope

The compiled corollary claims

\[
  N2^{-m}\le \rho^h
  \qquad
  (h\ge\kappa\log_2N).
\]

For \(0<\rho<1\), the function \(\rho^h\) decreases with \(h\). The cutoff
condition gives

\[
  \rho^h
  \le
  \rho^{\kappa\log_2N},
\]

which is the opposite direction from the proof. A polynomial term cannot be
bounded by \(\rho^h\) uniformly through \(h=\Theta(N)\).

The same issue affects the comparison of a fixed exponential-in-\(N\)
binomial tail with \(\rho^h\) unless their exponential rates are compared.

### Input runs do not imply distinct live episodes

The forced-termination proof selects \(r\) input runs and treats them as
\(r\) distinct live-state episodes. A new input run starts a new episode only
when the convolution state is zero. If the previous episode survives the
intervening zero gap, several input runs belong to one episode.

For example, the two runs in an input prefix \(101\) need not produce two
episodes. Therefore, the product of \(r\) termination probabilities is not
justified by the number of input runs alone.

The preserved scalar-dense refinement notes recognize this episode-clustering
problem and develop gap-aware machinery. That machinery is not integrated
into the compiled asymptotic theorem.

## Relation to Expand--Convolute

Expand--Convolute already proves high minimum distance for a related
randomized construction. Its outer matrix is a sparse Bernoulli expander,
and its convolution analysis uses a finite-state Markov chain. The expander
randomness directly determines the transition law for each fixed message.
The source is Raghuraman, Rindal, and Tanguy,
[Expand--Convolute Codes for Pseudorandom Correlation Generators from
LPN](https://eprint.iacr.org/2023/882).

The published asymptotic theorem uses

\[
  w=C\ln N,
  \qquad
  m=C_m\log N,
  \qquad
  C>2,
  \qquad
  C_m>1.
\]

Its analytic rate condition includes

\[
  R<
  \frac1{\ln2}\frac{e-1}{e+1}
  \left(\frac{20}{41}-\delta\right)^2.
  \tag{9}
\]

Even at \(\delta=0\), the right side of (9) is less than \(0.16\).
Consequently, the published analytic theorem does not prove the rate-half
case. The paper's rate-half tables evaluate a tighter finite formula.

The reusable part is the non-wrapping convolution proof. All \(m\) feedback
coefficients are random. For a fixed message under the Bernoulli expander,
the convolution is an \((m+1)\)-state homogeneous Markov chain. The proof
couples this chain to a reversible two-state chain when \(m\ge\log N+2\).
This route avoids the fixed-tap run-to-episode defect identified above.

Random SPIN instead uses independent random block injections followed by a
separate global permutation. Conditioning on the outer weight makes the
interleaved word uniform on a Hamming slice. The convolution proof must
therefore control a without-replacement support process. The
Expand--Convolute theorem is a strong proof template, but it does not directly
imply the Random SPIN theorem.

## Exact without-replacement convolution chain

The non-wrapping inner gives a finite chain that retains the global
interleaver exactly. Fix an outer weight \(h\), and let \(U\) be uniform on the
weight-\(h\) Hamming slice.

Before time \(t\), let \(R_t\) be the number of input ones that remain in
positions \(t,\ldots,N\). Let \(Z_t\in\{0,\ldots,m\}\) be the number of
trailing output zeros, clipped at \(m\). Thus \(Z_t=m\) exactly when the
convolution state is zero.

Conditioned on the complete history before time \(t\),

\[
  \Pr[U_t=1]=\frac{R_t}{N-t+1}.
  \tag{10}
\]

After sampling \(U_t\), the output law is

\[
  \Pr[Y_t=1\mid U_t,Z_t]
  =
  \begin{cases}
    1/2,&Z_t<m,\\
    U_t,&Z_t=m.
  \end{cases}
  \tag{11}
\]

The updates are

\[
  R_{t+1}:=R_t-U_t,
  \qquad
  Z_{t+1}:=
  \begin{cases}
    0,&Y_t=1,\\
    \min\{m,Z_t+1\},&Y_t=0.
  \end{cases}
  \tag{12}
\]

Equations (10)--(12) define an exact time-inhomogeneous Markov chain on
\((R_t,Z_t)\). The proof is direct. Sampling without replacement gives (10).
When \(Z_t<m\), the convolution state is nonzero, so
\(\langle\alpha_t,s_t\rangle\) is a fresh fair bit. When \(Z_t=m\), the state
is zero and (I) gives \(Y_t=U_t\).

This chain is the preferred proof interface. It can be evaluated exactly at
finite lengths and compared with the homogeneous Expand--Convolute chain
asymptotically. It does not count input runs or assume that they create
distinct episodes.

## Exact transfer matrix

The same chain admits a smaller coefficient-generating representation. This
representation removes the remaining-input coordinate from the state and
records the input weight with a formal variable.

Let \(X\) mark input weight, let \(s\) mark output weight, and index rows and
columns by \(z\in\{0,\ldots,m\}\). Define the nonnegative matrix
\(T_m(X,s)\) by

\[
  T_m(X,s)[z,z']
  :=
  \begin{cases}
    (1+X)s/2,&z<m,\ z'=0,\\
    (1+X)/2,&z<m,\ z'=z+1,\\
    Xs,&z=m,\ z'=0,\\
    1,&z=m,\ z'=m,\\
    0,&\text{otherwise}.
  \end{cases}
  \tag{13}
\]

Let \(e_m\) be the unit vector at the initial zero-state index \(m\), and let
\(\boldsymbol 1\) be the all-ones vector.

**Lemma 2 (exact bivariate transform).** For the inner map in (I),

\[
  \mathbb E_{I_N}\!\left[
    \sum_{u\in\mathbb F_2^N}
    X^{\operatorname{wt}(u)}s^{\operatorname{wt}(I_N(u))}
  \right]
  =
  e_m^{\mathsf T}T_m(X,s)^N\boldsymbol 1.
  \tag{14}
\]

*Proof.* Suppose first that \(z<m\). The convolution state is nonzero, so the
fresh random row \(\alpha_t\) makes the output bit fair for either input bit.
Summing the two input choices with weights \(1\) and \(X\) gives the first two
rows of cases in (13). If \(z=m\), the state is zero and \(Y_t=U_t\); input
zero stays at \(m\) with weight one, while input one moves to zero with weight
\(Xs\). Multiplying the one-step transforms from the initial state and summing
over the terminal state proves (14). \(\square\)

Consequently, if \(U\) is uniform on the weight-\(h\) slice, then

\[
  \mathbb E_{U,I_N}\!\left[
    s^{\operatorname{wt}(I_N(U))}
  \right]
  =
  \frac{
    [X^h]e_m^{\mathsf T}T_m(X,s)^N\boldsymbol 1
  }{\binom Nh}.
  \tag{15}
\]

For \(0<s<1\), an integer threshold \(D\), and any \(a>0\), the lower-tail
probability therefore satisfies

\[
  \begin{aligned}
  \Pr_{U,I_N}\!\left[
    \operatorname{wt}(I_N(U))\le D
  \right]
  &\le
  s^{-D}
  \frac{
    [X^h]e_m^{\mathsf T}T_m(X,s)^N\boldsymbol 1
  }{\binom Nh}\\
  &\le
  s^{-D}a^{-h}
  \frac{
    e_m^{\mathsf T}T_m(a,s)^N\boldsymbol 1
  }{\binom Nh}.
  \end{aligned}
  \tag{16}
\]

The first inequality is the exponential-moment bound
\(\mathbf 1\{W\le D\}\le s^{W-D}\). The second uses nonnegative
coefficients. Equations (14)--(16) are exact or rigorous finite bounds; they
do not invoke an independence approximation for the interleaved support.

The sparse saddle regime already admits an analytic spectral bound. Fix
\(a\in(0,1)\) and \(s\in(0,1)\), and set

\[
  b:=\frac{1+a}{2}.
\]

**Lemma 3 (Perron-root perturbation).** Suppose \(b(1+s)<1\), and let
\(\lambda_m(a,s)\) be the Perron root of \(T_m(a,s)\). Then
\(\lambda_m(a,s)>1\) and

\[
  1
  =
  \frac{bs}{\lambda_m}
  \frac{1-(b/\lambda_m)^m}{1-b/\lambda_m}
  +
  \left(\frac{b}{\lambda_m}\right)^m
  \frac{as}{\lambda_m-1}.
  \tag{17}
\]

In particular,

\[
  0<\lambda_m-1
  \le
  \frac{as(1-b)}{1-b(1+s)}\,b^m.
  \tag{18}
\]

*Proof.* The matrix is irreducible, has a self-loop of weight one at state
\(m\), and has a positive cycle leaving and returning to that state. Hence
its Perron root is greater than one. Normalize a positive right Perron vector
by \(v_0=1\). Its eigenvector equations are

\[
  \lambda_m v_z=bs+bv_{z+1}
  \quad (z<m),
  \qquad
  (\lambda_m-1)v_m=as.
\]

Substitution from \(z=m-1\) back to \(z=0\) gives (17). The first term on the
right side of (17) is at most
\(bs/(1-b)\). Therefore the last term is at least

\[
  1-\frac{bs}{1-b}
  =
  \frac{1-b(1+s)}{1-b}.
\]

Rearranging and using
\((b/\lambda_m)^m\le b^m\) proves (18). \(\square\)

For \(m=\lceil\gamma\log_2N\rceil\), (18) gives

\[
  N(\lambda_m-1)=o(1)
  \quad\text{whenever}\quad
  \gamma>\frac{1}{\log_2(1/b)}.
  \tag{19}
\]

Equation (19) controls the exponential growth of the matrix power in the
sparse saddle regime. A complete finite bound must also control the
Perron-vector or matrix-norm prefactor uniformly as \(m\) grows.

## Linear-weight exponent

The transfer matrix closes the linear-weight range analytically. For fixed
\(a>0\) and \(s\in(0,1)\), define

\[
  \lambda_\infty(a,s)
  :=
  \max\!\left\{
    1,\frac{(1+a)(1+s)}2
  \right\}.
  \tag{20}
\]

**Lemma 4 (limiting transfer growth).** Let
\(\lambda_m(a,s)\) be the Perron root of \(T_m(a,s)\). Then

\[
  \lim_{m\to\infty}\lambda_m(a,s)
  =
  \lambda_\infty(a,s).
  \tag{21}
\]

If \(m_N=\Theta(\log N)\), then

\[
  \limsup_{N\to\infty}
  \frac1N\log_2
  \left(
    e_{m_N}^{\mathsf T}
    T_{m_N}(a,s)^N
    \boldsymbol 1
  \right)
  \le
  \log_2\lambda_\infty(a,s).
  \tag{22}
\]

*Proof.* Set \(b=(1+a)/2\). If \(b(1+s)<1\), Lemma 3 gives
\(\lambda_m\to1\). If \(b(1+s)>1\), the principal active-state submatrix gives
\(\liminf_m\lambda_m\ge b(1+s)\). Equation (17) then gives the matching upper
limit. If \(b(1+s)=1\), any subsequential limit greater than one makes the
last term in (17) vanish and the first term strictly less than one. This
contradicts (17), so the critical limit is also one.

For (22), normalize a positive right Perron vector by \(v_0=1\). The
eigenvector equations used in Lemma 3 give

\[
  v_z
  =
  \frac{bs}{\lambda_m}
  \sum_{j=0}^{m-z-1}
  \left(\frac b{\lambda_m}\right)^j
  +
  \left(\frac b{\lambda_m}\right)^{m-z}v_m,
  \qquad
  v_m=\frac{as}{\lambda_m-1}.
\]

Equation (17) also gives
\(\lambda_m-1\ge as(b/\lambda_m)^m\). Hence the ratio between the largest and
smallest coordinates of \(v\) is \(\exp(O(m))\). Comparing
\(\boldsymbol 1\) with coordinatewise multiples of \(v\) bounds the matrix
power by \(\lambda_m^N\exp(O(m))\). Now use \(m_N=O(\log N)\) and (21).
\(\square\)

Define the binary entropy function

\[
  H_2(x):=-x\log_2x-(1-x)\log_2(1-x)
  \qquad (0<x<1).
\]

Fix \(\delta\in(0,1/2)\) and \(\eta\in(0,1)\). Equations (16), (20), and
(22) give the convolution exponent

\[
  F_\delta(\eta)
  :=
  \inf_{\substack{a>0\\0<s<1}}
  \left[
    \log_2\lambda_\infty(a,s)
    -\eta\log_2a
    -\delta\log_2s
    -H_2(\eta)
  \right].
  \tag{23}
\]

Set

\[
  \eta_{\mathrm{conv}}(\delta)
  :=
  1-\frac1{2(1-\delta)}.
\]

Direct differentiation gives

\[
  F_\delta(\eta)
  =
  \begin{cases}
    -\eta\log_2a_\eta
    -\delta\log_2s_\eta
    -H_2(\eta),
      &0<\eta\le\eta_{\mathrm{conv}}(\delta),\\[1mm]
    H_2(\delta)-1,
      &\eta_{\mathrm{conv}}(\delta)\le\eta<1,
  \end{cases}
  \tag{24}
\]

where

\[
  a_\eta
  :=
  \frac{\eta}{\sqrt{\eta^2+\delta^2}+\delta},
  \qquad
  s_\eta
  :=
  \frac{1-a_\eta}{1+a_\eta}.
  \tag{25}
\]

The first branch lies on
\((1+a_\eta)(1+s_\eta)=2\). For the second branch, the unconstrained
coefficient saddles are

\[
  a=\frac{\eta}{1-\eta},
  \qquad
  s=\frac{\delta}{1-\delta}.
\]

For an outer word \(u\) of weight \(h\), define

\[
  p_{N,h}(\delta)
  :=
  \Pr_{\Pi,I_N}\!\left[
    \operatorname{wt}(I_N(\Pi(u)))\le\delta N
  \right].
  \tag{26}
\]

The probability does not depend on the selected weight-\(h\) word.

**Theorem 1 (linear-weight closure at relative distance \(0.09\)).**
Let \(B_N=\Theta(\log N)\) be even and divide \(N\). Let
\(m_N=\Theta(\log N)\). For every sequence \(h_N\) with
\(h_N/N\to\eta\in(0,1)\),

\[
  \limsup_{N\to\infty}
  \frac1N\log_2
  \left(
    \mathbb E[A_{N,h_N}^{\mathrm{out}}]\,
    p_{N,h_N}(0.09)
  \right)
  \le
  J_{0.09}(\eta)
  <0,
  \tag{27}
\]

where

\[
  J_\delta(\eta)
  :=
  \Phi_{\mathrm{blk}}(\eta)+F_\delta(\eta).
  \tag{28}
\]

*Proof.* The first inequality combines (8), (16), and (22). It remains to
check the sign. Define

\[
  \eta_{\mathrm{blk}}:=1-2^{-1/2},
  \qquad
  \eta_{\mathrm{conv}}:=\eta_{\mathrm{conv}}(0.09).
\]

On \((0,\eta_{\mathrm{blk}}]\), the derivative of \(J_{0.09}\) is

\[
  \log_2\!\left(
    (1+\sqrt2)
    \frac{\sqrt{\eta^2+0.09^2}+0.09}{1-\eta}
  \right).
\]

This derivative is increasing, so \(J_{0.09}\) is convex. Also,
\(J_{0.09}(0)=0\). Interval arithmetic gives

\[
  J_{0.09}(\eta_{\mathrm{blk}})
  \in[-0.125943,-0.125942].
\]

Convexity therefore gives \(J_{0.09}(\eta)<0\) throughout the first region.
On
\([\eta_{\mathrm{blk}},\eta_{\mathrm{conv}}]\), the entropy terms cancel and
the remaining expression is increasing. Its right endpoint satisfies

\[
  J_{0.09}(\eta_{\mathrm{conv}})
  \in[-0.070598,-0.070597].
\]

For \(\eta\ge\eta_{\mathrm{conv}}\), equation (24) gives

\[
  J_{0.09}(\eta)
  =
  H_2(\eta)+H_2(0.09)-\frac32
  \le
  H_2(0.09)-\frac12
  \in[-0.063531,-0.063530].
\]

The script certify_random_spin_linear_exponent.py computes the displayed
intervals with mpmath interval arithmetic. \(\square\)

The theorem is uniform on every compact interval
\(\eta\in[\varepsilon,1-\varepsilon]\). The exponent tends to zero as
\(\eta\downarrow0\). More precisely,

\[
  J_{0.09}(\eta)
  =
  \eta\log_2\!\left(0.18(1+\sqrt2)\right)
  +o(\eta),
\]

and the linear coefficient lies in
\([-1.202378,-1.202377]\). Thus the unresolved range is sublinear input
weight, not linear input weight.

## Uniform sparse bound

The small-\(\eta\) saddle can be moved inside the stable Perron region at a
controlled cost. This gives one bound for every weight through \(N/25\).

**Lemma 5 (uniform sparse convolution bound).** Set

\[
  \delta:=\frac9{100},
  \qquad
  m_N:=\lceil2\log_2N\rceil.
\]

There are constants \(C>0\) and \(N_0\) such that, for all \(N\ge N_0\) and
all \(1\le h\le N/25\),

\[
  p_{N,h}(\delta)
  \le
  C N^4\sqrt h\,\rho_*^h,
  \qquad
  \rho_*:=
  \frac9{40}\exp\!\left(\frac{631}{2640}\right).
  \tag{29}
\]

Interval arithmetic gives

\[
  \rho_*
  \in[0.28574950,0.28574951]
  <
  \sqrt2-1.
  \tag{30}
\]

*Proof.* Fix \(h\) in the stated range and set

\[
  \eta:=\frac hN,
  \qquad
  \chi:=\frac32,
  \qquad
  a:=\frac{\eta}{\delta(\chi+1)},
  \qquad
  s:=\frac{1-\chi a}{1+a}.
\]

The range on \(h\) gives \(a\le8/45\) and \(s>0\). With
\(b=(1+a)/2\), the active-state row sum is

\[
  b(1+s)=1-\frac{\chi-1}{2}a<1.
\]

Lemma 3 therefore gives

\[
  \lambda_m-1
  \le
  \frac{2s(1-b)}{\chi-1}b^m
  \le
  2b^m.
\]

At the endpoint, \(b\le53/90\). The choice of \(m_N\) gives
\(Nb^{m_N}=o(1)\). The Perron-vector formula in Lemma 4 also gives

\[
  \frac{\max_zv_z}{\min_zv_z}\le C N^4.
\]

Indeed, equation (17) implies
\(\lambda_m-1\ge as(b/\lambda_m)^m\). For sufficiently large \(N\),
\(\lambda_m\le2\), so \(v_m\le4^m\le4N^4\). The remaining coordinates have
a constant lower bound. The schedule inequality

\[
  2\log_2(90/53)-1>0.5278
\]

also gives a constant lower bound for \(v_m\), including at \(h=1\).
Consequently,

\[
  e_m^{\mathsf T}T_m(a,s)^N\boldsymbol 1
  \le C N^4.
\]

For \(1\le h<N\), Stirling's upper bound on \(h!\) and
\(\log(1-x)\ge-x/(1-x)\) give

\[
  \binom Nh
  \ge
  \frac1{3\sqrt h}
  \left(\frac{eN}{h}\right)^h
  \exp\!\left(
    -\frac{h(h-1)}{2(N-h+1)}
  \right).
  \tag{31}
\]

For the inclusive integer threshold \(D=\lfloor\delta N\rfloor\),
\(s^{-D}\le s^{-\delta N}\).
Apply (16), (31), and the inequalities
\(\log(1+a)\le a\) and
\(-\log(1-\chi a)\le\chi a/(1-\chi a)\). The result is

\[
  p_{N,h}(\delta)
  \le
  C N^4\sqrt h\,
  [\delta(\chi+1)]^h
  \exp\!\left\{
    h\left[
      \frac{\chi^2a}{(\chi+1)(1-\chi a)}
      +
      \frac{\eta}{2(1-\eta)}
    \right]
  \right\}.
\]

Both correction terms increase with \(\eta\). At \(\eta=1/25\), their sum is

\[
  \frac{12}{55}+\frac1{48}
  =
  \frac{631}{2640}.
\]

Since \(\delta(\chi+1)=9/40\), equation (29) follows. \(\square\)

**Corollary 2 (sparse first-moment sum).** Suppose \(B_N\) is even, divides
\(N\), and satisfies

\[
  B_N\ge42\log_2N.
\]

Then

\[
  \sum_{h=1}^{\lfloor N/25\rfloor}
  \mathbb E[A_{N,h}^{\mathrm{out}}]\,
  p_{N,h}(0.09)
  =o(1).
  \tag{32}
\]

*Proof.* Apply Lemma 5 and use \(\sqrt h\le\sqrt N\). Equation (2) gives

\[
  \sum_{h\le N/25}
  \mathbb E[A_{N,h}^{\mathrm{out}}]\,
  p_{N,h}(0.09)
  \le
  C N^{9/2}
  \left[
    (1+G_{B_N}(\rho_*))^{N/B_N}-1
  \right].
\]

Set \(q_*:=(1+\rho_*)/\sqrt2\). Equations (3) and (30) give \(q_*<1\).
For all sufficiently large \(N\), the right side is at most a constant times

\[
  \frac{N^{11/2}}{B_N}q_*^{B_N}.
\]

The interval check gives

\[
  \log_2(1/q_*)>0.1373904,
  \qquad
  42\log_2(1/q_*)-\frac{11}{2}>0.27039.
\]

Thus the last display tends to zero. \(\square\)

## Complete admissible-length theorem

For each sufficiently large integer \(L\), let \(B_L\) be the least positive
even integer satisfying

\[
  B_L\ge42\log_2(LB_L).
\]

Set

\[
  N_L:=LB_L,
  \qquad
  K_L:=N_L/2,
  \qquad
  m_L:=\lceil2\log_2N_L\rceil.
\]

Sample the proof-model ensemble at these parameters.

**Theorem 2 (Random SPIN at relative distance \(0.09\)).** Over the independent
outer injections, global permutation, and convolution setup,

\[
  \Pr\!\left[
    d_{\min}(E_{N_L})\le0.09N_L
  \right]
  =o(1)
  \qquad (L\to\infty).
  \tag{33}
\]

The admissible family has rate \(1/2\). Its outer blocks and convolution
memory are both \(\Theta(\log N_L)\). Direct encoding uses
\(\Theta(N_L\log N_L)\) bit operations.

*Proof.* Split the first-moment sum into three ranges. Corollary 2 handles
\(h\le N_L/25\). Theorem 1 is uniform on
\(h/N_L\in[1/25,0.9]\), so that range contributes \(2^{-\Omega(N_L)}\).
For the remaining range, use \(p_{N,h}\le1\) and extract the outer tail at
\(z=9\). Equations (2) and (6) give

\[
  \sum_{h\ge0.9N_L}
  \mathbb E[A_{N_L,h}^{\mathrm{out}}]
  \le
  2^{N_L(H_2(0.9)-1/2+o(1))}
  =
  2^{-\Omega(N_L)}.
\]

Interval arithmetic gives
\(H_2(0.9)-1/2<-0.03100\). Hence the complete first moment tends to zero.
Markov's inequality proves (33). \(\square\)

The least-even definition gives \(B_L=\Theta(\log N_L)\) and
\(B_{L+1}-B_L=O(1)\). Therefore

\[
  N_{L+1}-N_L
  =
  O(N_L/\log N_L)
  =
  o(N_L).
\]

The padding and puncturing wrapper in ARBITRARY_LENGTH_WRAPPER.md extends the
family to every sufficiently large length with \(o(1)\) relative rate and
distance loss.

## Finite diagnostics

The script evaluate_random_spin_finite.py propagates the joint law of input
weight, output weight, and clipped zero-state. It then combines that law with
the exact outer spectrum. A separate integer recurrence matches exhaustive
enumeration for all lengths through five at memory two. The check compares 57
input-weight and output-weight coefficients.

The binary64 evaluator gives the following first-moment values at inclusive
distance \(\lfloor0.09N\rfloor\):

| \(N\) | \(B\) | \(m\) | outer blocks | \(\log_2\mu\) |
|---:|---:|---:|---:|---:|
| 512 | 16 | 24 | 32 | \(2.02265\) |
| 512 | 32 | 24 | 16 | \(-4.30265\) |
| 1024 | 32 | 28 | 32 | \(-3.28969\) |

The \(B=16\) control fails because weights one through eight remain too
large. At \(B=32\), weights four through six dominate both successful
points. These values validate the expected sparse obstruction and the exact
dynamic program. They are not outward-rounded finite certificates.

## Conservative first-moment conclusion

Fix a target \(\delta_R\in(0,1/2)\). The first-moment target is

\[
  \sum_{h=1}^N
  \mathbb E[A_{N,h}^{\mathrm{out}}]\,
  p_{N,h}(\delta_R)
  =o(1).
  \tag{34}
\]

Theorem 2 proves (34) at \(\delta_R=0.09\). The proof uses three ranges.

1. **Sparse and sublinear:** \(1\le h\le N/25\). Apply Lemma 5 and the exact
   outer generating function.
2. **Linear:** \(N/25<h<0.9N\). Apply Theorem 1 uniformly.
3. **Top:** \(0.9N\le h\le N\). Use the outer tail with \(p_{N,h}\le1\).

The conservative theorem leaves three optimization targets: reduce the outer
block constant, obtain an outward-rounded finite certificate, and move the
distance toward the rate-half Gilbert--Varshamov point

\[
  H_2(\delta_{\mathrm{GV}})=\frac12.
\]

## Optimized near-GV schedule

This section resolves the two asymptotic optimization targets. It keeps the
same probability space and changes only the numerical saddle and schedule
parameters.

Set

\[
  \delta_*:=\frac{5501}{50000}=0.11002.
\]

**Lemma 6 (linear closure at \(\delta_*\)).** Let
\(B_N=\Theta(\log N)\) be even and divide \(N\), and let
\(m_N=\Theta(\log N)\). For every sequence \(h_N/N\to\eta\in(0,1)\),

\[
  \limsup_{N\to\infty}\frac1N\log_2\left(
    \mathbb E[A_{N,h_N}^{\mathrm{out}}]\,
    p_{N,h_N}(\delta_*)
  \right)
  \le J_{\delta_*}(\eta)<0.
  \tag{35}
\]

The bound is uniform when \(\eta\) ranges over a compact subset of
\((0,1)\).

*Proof.* Repeat the proof of Theorem 1 with \(\delta_*\) in place of
\(0.09\). The analytic reductions are unchanged. Outward-rounded interval
arithmetic gives

\[
\begin{aligned}
  J_{\delta_*}(1-2^{-1/2})
    &\in[-0.074030,-0.074029],\\
  J_{\delta_*}(\eta_{\mathrm{conv}}(\delta_*))
    &\in[-0.011076,-0.011075],\\
  H_2(\delta_*)-\frac12
    &\in[-0.000023719,-0.000023718],\\
  \log_2(2\delta_*(1+\sqrt2))
    &\in[-0.912609,-0.912608].
\end{aligned}
\]

The first three strict inequalities close the three linear regions. The last
one is the negative right derivative at zero. The script
`certify_random_spin_linear_exponent.py` verifies these intervals. \(\square\)

The following parameterized statement records which constants are sufficient
for the sparse range. It distinguishes the mathematical conditions from the
particular numerical choice below.

**Lemma 7 (parameterized sparse closure).** Fix
\(\delta\in(0,1/2)\), a cutoff \(\eta_*\in(0,1)\), an interior parameter
\(\chi>1\), and a memory constant \(\gamma>0\). Define

\[
\begin{aligned}
  a_*&:=\frac{\eta_*}{\delta(\chi+1)},
  &b_*&:=\frac{1+a_*}{2},\\
  \epsilon_*&:=
  \frac{\chi^2a_*}{(\chi+1)(1-\chi a_*)}
  +\frac{\eta_*}{2(1-\eta_*)},
  &\rho_*&:=\delta(\chi+1)e^{\epsilon_*}.
\end{aligned}
\]

Suppose

\[
  \chi a_*<1,
  \qquad
  \gamma\log_2(1/b_*)>1.
  \tag{36}
\]

For \(m_N=\lceil\gamma\log_2N\rceil\), there are constants \(C,N_0\)
such that, for every \(N\ge N_0\) and \(1\le h\le\eta_*N\),

\[
  p_{N,h}(\delta)
  \le C N^{2\gamma}\sqrt h\,\rho_*^h.
  \tag{37}
\]

If also \(\rho_*<\sqrt2-1\), then the sparse first-moment sum is \(o(1)\)
for every even divisor schedule satisfying \(B_N\ge c\log_2N\), provided

\[
  c\log_2\!\left(\frac{\sqrt2}{1+\rho_*}\right)
  >2\gamma+\frac32.
  \tag{38}
\]

*Proof.* In Lemma 5, use
\(a=\eta/[\delta(\chi+1)]\) and
\(s=(1-\chi a)/(1+a)\). The active-state row sum is
\(1-(\chi-1)a/2<1\). Condition (36) makes
\(Nb_*^{m_N}=o(1)\) and bounds the Perron-vector ratio by
\(CN^{2\gamma}\). The same binomial lower bound then gives (37), because
the exponent correction is increasing in \(\eta\) and is at most
\(\epsilon_*\). Combining (37) with the exact outer generating function
leaves at most a constant times

\[
  \frac{N^{2\gamma+3/2}}{B_N}
  \left(\frac{1+\rho_*}{\sqrt2}\right)^{B_N}.
\]

Condition (38) makes this expression tend to zero. \(\square\)

For the optimized theorem, choose

\[
  \eta_*:=\frac1{5000},
  \qquad
  \chi:=\frac{101}{100},
  \qquad
  \gamma:=\frac{51}{50},
  \qquad
  c:=17.
\]

Here

\[
  a_*=\frac{1000}{1105701},
  \quad
  b_*=\frac{1106701}{2211402},
  \quad
  \epsilon_*=\frac{1241938871}{2219984824218}.
\]

Interval arithmetic certifies

\[
\begin{aligned}
  \rho_*&\in[0.2212639483446276,0.2212639483446278],\\
  \gamma\log_2(1/b_*)-1&>0.0186697,\\
  17\log_2\!\left(\frac{\sqrt2}{1+\rho_*}\right)
    -\left(2\gamma+\frac32\right)&>0.0576243.
\end{aligned}
\tag{39}
\]

Thus every inequality in Lemma 7 is strict.

For each sufficiently large integer \(L\), let \(B_L^*\) be the least
positive even integer satisfying

\[
  B_L^*\ge17\log_2(LB_L^*).
\]

Set

\[
  N_L^*:=LB_L^*,
  \qquad
  K_L^*:=N_L^*/2,
  \qquad
  m_L^*:=\left\lceil\frac{51}{50}\log_2N_L^*\right\rceil.
\]

**Theorem 3 (optimized Random SPIN).** Sample the proof-model ensemble using
the preceding admissible parameters. Over the independent outer injections,
global permutation, and convolution setup,

\[
  \Pr\!\left[d_{\min}(E_{N_L^*})\le\delta_*N_L^*\right]
  =o(1)
  \qquad(L\to\infty).
  \tag{40}
\]

The family has rate \(1/2\), relative distance greater than \(0.11002\),
outer block size \(\Theta(\log N_L^*)\), and convolution memory
\(\Theta(\log N_L^*)\). Direct encoding uses
\(\Theta(N_L^*\log N_L^*)\) bit operations.

*Proof.* Lemma 7 and (39) close
\(1\le h\le N_L^*/5000\). Lemma 6 is uniform on
\(h/N_L^*\in[1/5000,0.9]\), so the middle range contributes
\(2^{-\Omega(N_L^*)}\). The outer-tail argument in Theorem 2 closes
\(h\ge0.9N_L^*\), independently of \(\delta_*\). Markov's inequality proves
(40). \(\square\)

The least-even rule gives \(B_L^*=\Theta(\log N_L^*)\) and relative gaps
\((N_{L+1}^*-N_L^*)/N_L^*=o(1)\). The same padding and puncturing wrapper
therefore extends Theorem 3 to every sufficiently large requested length,
with \(o(1)\) relative rate and distance loss.

## Status

- **Proved:** the exact random-block identities (1) and (2), the logarithmic
  outer schedule (3)--(5), the outer exponent (6)--(8), and the exact
  without-replacement chain (10)--(12), including its bivariate transfer
  identity and coefficient bounds (13)--(16), and the Perron-root bound
  (17)--(19).
- **Proved:** the limiting transfer growth (20)--(22), the optimized
  convolution exponent (23)--(25), and linear-weight closure at relative
  distances \(0.09\) and \(0.11002\).
- **Proved:** the uniform sparse bound (29)--(32), the complete
  conservative theorem (33), the parameterized sparse conditions
  (36)--(38), the optimized theorem (40), and the arbitrary-length
  extensions.
- **Imported proof ingredients:** the Hamming-slice reduction and local
  convolution lemmas listed above.
- **Invalid as written:** the constant geometric envelope and the
  run-to-episode product step.
- **External neighboring theorem:** Expand--Convolute proves a sparse-expander
  plus convolution ensemble.
- **Open optimization:** close the final \(7.9\times10^{-6}\) to the
  rate-half Gilbert--Varshamov distance, further reduce the block constant,
  and replace the binary64 finite diagnostics with outward-rounded
  certificates.
