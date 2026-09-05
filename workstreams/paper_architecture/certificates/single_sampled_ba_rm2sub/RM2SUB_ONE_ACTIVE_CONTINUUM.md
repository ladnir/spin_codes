# One-active continuum for the fixed RM2Sub-S19 inner

## Result

Consider the random rate-half outer ensemble and fixed \(t=128,s=19\)
RM2Sub inner in `RANDOM_OUTER_RM2SUB_RAMP.md`. Let \(Z_{D,1}\) count messages
supported on one outer block whose output weight is at most \(D\).

For \(D=\lfloor0.11002N\rfloor\) and \(B=c\log_2N+O(1)\), the optimized
one-active continuum gives

\[
  c>3.587179987\ldots .
  \tag{1}
\]

The decimal in (1) uses binary64 optimization. An exact rational inequality
certifies the convenient choice

\[
  B\ge\frac{18}{5}\log_2N+O(1).
  \tag{2}
\]

Under (2),

\[
  \Pr[Z_{\lfloor0.11002N\rfloor,1}>0]
  \le N^{-0.0026858\ldots+o(1)}.
  \tag{3}
\]

Equation (3) controls occupation one. The general fixed-occupation theorem is
in `RM2SUB_FIXED_OCCUPATION_CONTINUUM.md`.

## Why the finite two-state envelope is insufficient

The frozen finite evaluator uses a zero/live envelope. Its zero-input live
entry pays

\[
  \kappa:=\frac{2^{19}-1}{2^{19}-2}>1.
  \tag{4}
\]

A region contains \(E=L/128\) epochs. Direct iteration therefore pays
\(\kappa^E\). This factor is negligible at the frozen \(E=64\), but it grows
exponentially when \(L\to\infty\).

The construction does not have this repeated loss. A nonuniform live state
persists for at most one zero-input epoch. The next fresh field multiplier
makes the state uniform on the nonzero field elements. The asymptotic proof
uses the exact \(2^{19}\)-state kernel before reducing to two states.

## Exact epoch kernel

Set

\[
  M:=2^{19}-1.
  \tag{5}
\]

Let \(A:\mathbb F_2^{19}\to\mathbb F_2^{128}\) and
\(C:\mathbb F_2^{128}\to\mathbb F_2^{19}\) be the audited RM2Sub maps. Every
coordinate form of \(A\) is nonzero, and every column \(C(e_j)\) is nonzero.

For \(a\in\{0,1\}\), let \(X_0:=0\), and let \(X_1:=e_J\) for
\(J\gets[128]\). For \(z\in(0,1)\), define

\[
  W_a(z)[q,q']
  :=
  \mathbb E_{J,\alpha}\!\left[
    z^{\operatorname{wt}(X_a+A(q))}
    \boldsymbol1\{q'=\alpha q+C(X_a)\}
  \right],
  \tag{6}
\]

where \(\alpha\gets\mathbb F_{2^{19}}^*\). Equation (6) is the exact
\(2^{19}\)-state epoch kernel.

Let \(E=L/128\). Define

\[
  \mathcal R_{0,L}(z):=W_0(z)^E
  \tag{7}
\]

and

\[
  \mathcal R_{1,L}(z)
  :=
  \frac1E\sum_{j=0}^{E-1}
  W_0(z)^jW_1(z)W_0(z)^{E-1-j}.
  \tag{8}
\]

The region permutation makes (8) exact for one active outer block.

## Long-region limit

Fix \(\theta>0\), and set

\[
  z_L:=e^{-\theta/L}.
  \tag{9}
\]

For a uniform nonzero state, each coordinate of \(A(q)\) is one with
probability

\[
  p:=\frac{2^{18}}{2^{19}-1}.
  \tag{10}
\]

Define

\[
  \gamma:=p\theta,
  \qquad
  a:=e^{-\gamma},
  \qquad
  f:=\frac{1-e^{-\gamma}}{\gamma}.
  \tag{11}
\]

The limiting inactive and marked-region matrices are

\[
  K_0(\theta)=
  \begin{pmatrix}
    1&0\\
    0&a
  \end{pmatrix}
  \tag{12}
\]

and

\[
  K_1(\theta)=
  \begin{pmatrix}
    0&f\\
    f/M&(1-1/M)a
  \end{pmatrix}.
  \tag{13}
\]

**Lemma 1 (marked-region limit).** Suppose \(128\mid L\) and
\(B=O(\log L)\). The exact full-state transfers (7)--(8), composed across
\(B\) regions and started from state zero, have the same logarithmic rate in
\(B\) as (12)--(13).

*Proof.* A zero state remains zero under a zero-input epoch. A uniform live
state remains uniform live, and

\[
  \mathbb E[z_L^{\operatorname{wt}(A(q))}]
  =1-\frac{128p\theta}{L}+O(L^{-2}).
  \tag{14}
\]

Raising (14) to \(E=L/128\) gives \(a\).

An impulse entering state zero creates the nonzero state \(C(e_J)\). One
following zero-input epoch refreshes that state. The transient output factor
is \(1+O(L^{-1})\). Averaging the remaining live suffix gives

\[
  \int_0^1e^{-\gamma x}\,dx=f.
  \tag{15}
\]

For a uniform live entering state, the impulse terminates it with probability
\(1/M\). The live prefix gives \(f/M\). On survival, the prefix and suffix
cover all but one transient epoch, which gives \((1-1/M)a\).

An impulse in the final epoch can leave a nonuniform boundary state. This
event has probability \(1/E\). Resolving the boundary state changes a
nonzero limiting entry by relative \(1+O(M/E)\). Since \(B=o(E)\), the total
logarithmic correction is \(o(B)\). \(\square\)

## Spectral exponent

Let

\[
  T_1(\theta):=K_0(\theta)+K_1(\theta)
  =
  \begin{pmatrix}
    1&f\\
    f/M&(2-1/M)a
  \end{pmatrix}.
  \tag{16}
\]

Its Perron root is

\[
  \lambda_1(\theta)
  =
  \frac{
    1+(2-1/M)a+
    \sqrt{\big(1-(2-1/M)a\big)^2+4f^2/M}
  }2.
  \tag{17}
\]

For \(\delta\in(0,1/2)\), define

\[
  \eta_1(\delta)
  :=
  \frac12-
  \inf_{\theta>0}
  \left(
    \log_2\lambda_1(\theta)+
    \frac{\delta\theta}{\ln2}
  \right).
  \tag{18}
\]

**Theorem 2 (occupation-one exponent).** Let \(B\to\infty\) be even. Let
\(L\to\infty\) be divisible by \(128\), and assume \(B=O(\log L)\). Then

\[
  \log_2\mathbb E[Z_{\lfloor\delta LB\rfloor,1}]
  \le
  \log_2L-\eta_1(\delta)B+o(B).
  \tag{19}
\]

*Proof.* The exact random-outer multiplicity is

\[
  \rho_B:=\frac{2^{B/2}-1}{2^B-1},
  \qquad
  \log_2\rho_B=-\frac B2+o(1).
  \tag{20}
\]

Lemma 1 gives the region rate \(\lambda_1(\theta)^B\). The Chernoff factor is

\[
  z_L^{-\lfloor\delta LB\rfloor}
  =2^{\delta\theta B/\ln2+o(B)}.
  \tag{21}
\]

Multiplying these factors and minimizing over \(\theta\) proves (19).
\(\square\)

For \(B=c\log_2N+O(1)\), Markov's inequality gives

\[
  \Pr[Z_{\lfloor\delta N\rfloor,1}>0]
  \le N^{1-c\eta_1(\delta)+o(1)}.
  \tag{22}
\]

Thus occupation one closes when \(c>1/\eta_1(\delta)\).

## Constants and exact certificate

At \(\delta=0.11002\), binary64 optimization gives

\[
  \theta=1.388982314\ldots,
  \qquad
  \eta_1=0.2787705115\ldots,
  \qquad
  \frac1{\eta_1}=3.587179987\ldots .
  \tag{23}
\]

For an exact certificate, choose

\[
  \theta_0:=\frac{\ln2}{p}.
  \tag{24}
\]

Then \(a=1/2\) and \(f=1/(2\ln2)\). The receipt bounds

\[
  \ln2=2\sum_{k\ge0}\frac{3^{-(2k+1)}}{2k+1}
  \tag{25}
\]

with exact rational partial sums. It also bounds the square root in (17) with
integers and uses \(\ln(1+x)\le x\). The resulting exact inequality is

\[
  \eta_1(0.11002)>0.2785238501.
  \tag{26}
\]

Therefore

\[
  \frac{18}{5}\eta_1(0.11002)-1
  >0.0026858604,
  \tag{27}
\]

which proves (2)--(3).

## Admissible lengths and status

Set \(L_m:=128m\). Let \(B_m\) be the least positive even integer satisfying

\[
  B_m\ge\frac{18}{5}\log_2(L_mB_m).
  \tag{28}
\]

Then \(N_m:=L_mB_m\), \(B_m=\Theta(\log N_m)\), and
\((N_{m+1}-N_m)/N_m=o(1)\).

- **Proved here:** the exact epoch kernel, the long-region matrices, the
  occupation-one exponent, and the rational \(c=18/5\) certificate.
- **Numerical only:** the optimized threshold in (1).
- **Open beyond this note:** a uniform growing-occupation bound, the bulk
  bridge, implementation equivalence, and the structured-outer comparison.

The arbitrary-length wrapper cannot be invoked for full minimum distance
until the remaining occupation regimes are closed.
