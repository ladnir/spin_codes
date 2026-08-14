# Audit of a rank-two conditioned-row outer lemma

## Scope

The current implementation accepts `conditioned_rows=2`, but then replaces
the two retained bits by their XOR parity. It maximizes the even states
$00,11$, maximizes the odd states $01,10$, and frees the other retained
codeword. That path is an XOR-collapse relaxation, not a rank-two lemma.

The algebra below retains the three possible sums of two ordered bits. It then
uses an entrywise rank-one envelope. The algebra is not a bound for the frozen
sloped construction with independent puncture lanes. It becomes applicable
after the construction adopts the global-lane puncture rule described in
`conditioned_row_sloped_layout_audit.md`.

## Probability space and local factors

Fix positive packet fugacities $t=(t_0,\ldots,t_g)$, with $t_0=1$. Let
$R_w(t)$ be the normalized packet moment of a 64-bit column of weight $w$.
Fix linear-BL coefficients $p_b\in(0,1]$ whose pairwise sums are at least
one.

In each tile, choose two distinct data rows. Their independent messages are
uniform in $\mathbb F_2^{64}$, so their encoded rows are independent uniform
codewords $C,D$ in the extended BCH code $\mathcal C$. The other 62
messages remain independent and uniform. Conditional on retained bits
$x,y\in\{0,1\}$ at a coordinate in band $b$, linear BL gives

\[
f_{b,x+y}(t):=
\left(
2^{-62}\sum_{u=0}^{62}{62\choose u}
R_{u+x+y}(t)^{1/p_b}
\right)^{p_b}.
\tag{B1}
\]

The BL dimension condition is unchanged because fixing two rows removes two
direct-sum copies of the same message maps. The free-message factor is
$2^{62\cdot64}$ per tile.

This probability space requires two global lane matchings. After the
entrywise separation below, each matching must regroup into complete EBCH
codewords. The frozen independent-lane puncture sampler does not provide the
first required matching through all punctures. Thus the hole formula below is
inapplicable to that sampler.

If one global lane `L*` contains all punctured blocks, choose another lane
`L' != L*` for the second row. Each lane class is a perfect matching in every
band. The separated factors then regroup globally even though one codeword's
three band cells occupy shifted physical tiles.

## Rank-one envelope

For one band, write $\ell_s:=\log f_{b,s}$ for $s=0,1,2$, and form

\[
G_b:=
\begin{pmatrix}
f_{b,0}&f_{b,1}\\
f_{b,1}&f_{b,2}
\end{pmatrix}.
\]

Choose any $d_b\in\mathbb R$, and define

\[
a_b:=\max\left\{
\frac{\ell_0}{2},
\frac{\ell_2}{2}-d_b,
\frac{\ell_1-d_b}{2}
\right\},
\qquad
u_{b,0}:=e^{a_b},\qquad u_{b,1}:=e^{a_b+d_b}.
\tag{B2}
\]

The three terms in the maximum are exactly the constraints

\[
f_{b,0}\le u_{b,0}^2,\qquad
f_{b,1}\le u_{b,0}u_{b,1},\qquad
f_{b,2}\le u_{b,1}^2.
\tag{B3}
\]

Therefore $G_b(x,y)\le u_{b,x}u_{b,y}$ entrywise. Multiplying (B3) over
coordinates separates the retained codewords $C$ and $D$.

## One-row enumerators and Cauchy

Put $c=(42,43,43)$ and $r_b:=u_{b,1}/u_{b,0}=e^{d_b}$. Define

\[
\mathcal B_{\rm n}(u;\theta):=
\left(\prod_bu_{b,0}^{c_b}\right)
\sqrt{
A_{01}(r_0^2,r_1^{2\theta})
\,A_{12}(r_1^{2(1-\theta)},r_2^2)
},
\tag{B4}
\]

where $A_{01}$ and $A_{12}$ are the existing exact split enumerators.
For every $\theta\in[0,1]$, Cauchy--Schwarz gives

\[
\sum_{c\in\mathcal C}\prod_{b,k\in b}u_{b,c_k}
\le \mathcal B_{\rm n}(u;\theta).
\]

Consequently, a normal tile contributes at most

\[
T_{\rm n}:=
2^{62\cdot64}\mathcal B_{\rm n}(u;\theta_{\rm n})^2.
\tag{B5}
\]

The square is the product of two independently summed one-row enumerators.
It is not an extra free-message factor.

## Puncture and graph holes

In a hole tile, choose the punctured data row as the first retained row. Choose
one of the other 63 rows as the second retained row. The puncture removes the
first row's band-zero bit. If the graph bit is $g\in\{0,1\}$ and the second
row's bit is $y$, the hole factor is $f_{0,g+y}$. Equation (B3) gives

\[
f_{0,g+y}\le u_{0,g}u_{0,y}.
\tag{B6}
\]

Thus the first codeword uses a punctured one-row enumerator, the second uses a
normal one-row enumerator, and the graph bit contributes $u_{0,g}$.

Let $P_{01}$ be the existing puncture-averaged band-$(0,1)$ enumerator,
including its divisor 42. Define

\[
\mathcal B_{\rm p}(u;\theta):=
u_{0,0}^{41}u_{1,0}^{43}u_{2,0}^{43}
\sqrt{
P_{01}(r_0^2,r_1^{2\theta})
\,A_{12}(r_1^{2(1-\theta)},r_2^2)
}.
\tag{B7}
\]

For graph bit $g$, the hole-tile contribution is at most

\[
T_{{\rm h},g}:=
2^{62\cdot64}
\mathcal B_{\rm p}(u;\theta_{\rm p})
\mathcal B_{\rm n}(u;\theta_{\rm n})u_{0,g}.
\tag{B8}
\]

Let $A_w^{\rm graph}$ be the stored graph-code weight spectrum. The complete
outer moment is bounded by

\[
T_{\rm n}^{128}
2^{-24}\sum_{w=0}^{128}A_w^{\rm graph}
T_{{\rm h},0}^{128-w}T_{{\rm h},1}^{w}.
\tag{B9}
\]

Equation (B9) uses the same graph-syndrome uniformity, distinct-hole sampling,
and graph-coordinate bijection as the one-row proof. A verifier must bind
those assumptions rather than treating 128 graph bits independently.

## Normalization and coefficient extraction

At unit fugacity, every factor in (B1) equals one. Choosing $d_b=0$ gives
$a_b=0$ and $u_{b,0}=u_{b,1}=1$. Both one-row enumerators then have mass
$2^{64}$. Hence every normal and hole tile has mass

\[
2^{62\cdot64}2^{64}2^{64}=2^{4096}.
\]

The graph spectrum in (B9) has normalized mass one. The global mass is
$2^{256\cdot4096}=2^K$. This identity checks that the two retained rows were
summed exactly once. It also forbids an additional $2^{128}$ factor.

For a packet profile $v$, subtract
$\langle v,\log_2t\rangle$ once from the logarithm of (B9). The packet-profile
normalization is not part of either retained-row enumerator.

## Strength and required artifacts

The bound is genuinely different from the current parity collapse. It uses
all three local values $f_{b,0},f_{b,1},f_{b,2}$ before separating the two
codewords. The parity collapse instead replaces the diagonal pair by
$\max(f_{b,0},f_{b,2})$ and retains only $C\oplus D$. Neither relaxation
uniformly dominates the other. A proof system should retain the minimum of
the one-row, parity-collapse, and rank-one-envelope rows.

A strict improvement is possible when the fitted matrix $G_b$ is close to
rank one and the even-state maximum is wasteful. Equality in the rank-one
step occurs when $f_{b,1}^2=f_{b,0}f_{b,2}$ and the chosen $d_b$ realizes
that factorization. No theorem guarantees a strict improvement after the two
Cauchy bounds.

The rank-one-envelope lemma needs no new spectra. It reuses the exact
band-$(0,1)$, band-$(1,2)$, puncture-averaged band-$(0,1)$, and graph
weight spectra.

An exact rank-two tile sum without (B3) would require new second-order split
complete weight enumerators. For each ordered pair $(C,D)\in\mathcal C^2$,
those artifacts must record the four symbol counts
$n_{b,00},n_{b,01},n_{b,10},n_{b,11}$ in both bands of each Cauchy pair.
The punctured version must also mark the removed band-zero symbol. Ordinary
split weight spectra do not determine these second-order enumerators.

Before any diagnostic, static gates should verify (B3), all source masses,
the divisor-42 puncture identity, the graph mass $2^{24}$, and the unit-
fugacity mass $2^K$.
