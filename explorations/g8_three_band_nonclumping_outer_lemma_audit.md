# Audit of the three-band nonclumping outer lemma

## Verdict

Lemma 1 and Corollary 2 are valid after the local-factor interface is made
precise. The proposed diagnostic must not run from the current plan as
written. A packet-profile point does not determine message support size, and
the existing star and motif tables do not contain the required local moments.

The following formulation removes both ambiguities.

## Local factor and exact envelope

Fix a band of size $c\in\{42,43\}$, one tile, and all coordinate
permutations. For the 64 messages incident to that tile, let

\[
G_t(x_1,\ldots,x_{64})
 :=\prod_{k=1}^{c}R_{w_k}(t),
\]

where $w_k$ is the number of ones in coordinate $k$ among the 64 projected
codewords. Here

\[
R_w(t):=
\frac{[z^w](\sum_{j=0}^{8}{8\choose j}t_jz^j)^8}{{64\choose w}}.
\]

This is the tile monomial after averaging the independent lane permutation.
The global unpunctured outer moment factors as a product of these tile
monomials before summing over the independent messages. If the construction
also averages coordinate permutations, those random choices are local to one
tile factor and may be included in the product probability space.

For $0\le d\le r\le64$, define

\[
A_d(t):=\sum_{w=0}^{d}{d\choose w}R_w(t)^3
\]

and

\[
\overline m_{c,r}(t):=
\frac{1}{(2^{64}-1)^r}
\sum_{d=0}^{r}(-1)^{r-d}{r\choose d}
2^{(64-c)d}A_d(t)^c.
\tag{A1}
\]

Every projection from one 64-bit EBCH message to a fixed band is surjective.
For $d$ unrestricted messages, each projected $d$-tuple has
$2^{(64-c)d}$ preimages. The $c$ projected coordinates then separate, and
their common sum is $A_d(t)$. Inclusion-exclusion over the events $x_i=0$
proves

\[
\overline m_{c,r}(t)
=\mathbb E[G_t(X_1,\ldots,X_r,0,\ldots,0)^3],
\]

where the $X_i$ are independent uniform nonzero messages. The value is
independent of the tile and of the selected $r$ incident blocks. Thus the
subset maximum in the original definition is unnecessary for fixed tile
monomials.

The note instead defines $F_t:=\mathbb E_\pi[G_{t,\pi}]$ by averaging
coordinate permutations before taking the cube. In that convention, (A1)
need not equal $\mathbb E[F_t^3]$. Convexity gives the sufficient bound

\[
\mathbb E[F_t^3]
\le \mathbb E_{X,\pi}[G_{t,\pi}^3]
=\overline m_{c,r}(t).
\tag{A2}
\]

Accordingly, the theorem should either use fixed $G_t$, or replace every
$m_{b,r}$ by the envelope $\overline m_{c_b,r}$. Equating (A1) with the cube
moment of the pre-averaged factor is not justified.

## Finner and normalization checks

After conditioning on an activity set $S$, the active message variables are
independent. Each message occurs in one tile factor from each of the three
bands. Generalized Hölder therefore applies with exponent three to every tile
factor. Variables local to one factor, such as a tile permutation, do not
violate the fractional-cover condition.

The support decomposition also has the correct multiplier. There are
$(2^{64}-1)^{|S|}$ messages with activity set $S$, while the conditional
expectation uses uniform nonzero active values.

At $t=(1,\ldots,1)$, every $R_w(t)=1$. Hence $A_d(t)=2^d$, and (A1) becomes

\[
\frac{1}{(2^{64}-1)^r}
\sum_{d=0}^{r}(-1)^{r-d}{r\choose d}2^{64d}=1.
\]

Summing the support decomposition gives
$(1+2^{64}-1)^{|I|}=2^{64|I|}=2^K$. This check also detects an accidental
second message-count factor.

## Corollary 2

The coefficient algebra in Corollary 2 is correct. Pair capacity supplies

\[
1\le \exp\!\left(\eta\left[{s\choose2}-\sum_bT_b(S)\right]\right).
\]

To apply Hölder across the three bands, place $\exp(-3\eta T_b(S))$ inside
the band-$b$ sum. Taking its cube root recovers
$\exp(-\eta\sum_bT_b(S))$. Therefore the polynomial must contain the
unrooted local moment and the factor $\exp(-3\eta{r\choose2})$, exactly as
stated. An implementation may use base-two exponentials only after rescaling
$\eta$ consistently.

## Support-size summation

The nine packet counts of a leaf vertex describe the encoded word after the
outer map. They do not determine

\[
s:=|\{i\in I:x_i\ne0\}|,
\]

the number of nonzero 64-bit outer-message blocks. Different messages can
have the same packet-weight profile and different values of $s$. A bound for
a packet-profile coefficient must therefore sum the right-hand side of
Corollary 2 over every $s\in\{0,\ldots,16384\}$. Restriction to selected
support sizes is sound only after a separate certificate proves that the
target coefficient has no contribution from the omitted sizes. A dual
packet-profile barycenter is not such a certificate.

## Artifact and integration gates

The current artifacts establish useful inputs but do not yet implement this
interface.

* `certify_three_band_tile_map.py` establishes tile size, pair capacity, and
  the four-block incidence counts.
* The fixed-band projection artifacts establish the surjectivity used in
  (A1). Their detailed weight counts are not needed for (A1).
* `certify_three_band_star.py`, `certify_three_band_small_motifs.py`, and
  `three_band_exact_moments.py` certify scalar support/OR moments. They do not
  serialize the nine-fugacity values $\overline m_{c,r}(t)$.
* `probe_packet8_three_band_profile_enumerator.py` computes $R_w(t)$ in
  binary64, but it neither evaluates (A1) outward nor performs the required
  sum over all support sizes.

Before a numerical result is treated as evidence, a new evaluator must bind
the frozen fugacities, evaluate (A1) without cancellation error, and sum all
support sizes. It must reproduce $\overline m_{c,r}(1)=1$ and the global
$2^K$ normalization. A proof artifact also needs outward coefficient
evaluation and the exact graph/puncture replacement. The existing scalar
star or motif tables cannot substitute for these gates.

## Rejected conditional linear-BL extension

A proposed extension replaces the three cube roots by linear-BL coefficients
$p_b\in(0,1]$. It is not valid after conditioning every active message to be
nonzero. The linear-BL inequality uses uniform measure on the full message
space. The punctured set $\mathbb F_2^{64}\setminus\{0\}$ is not a vector
space, and the same constant-one inequality need not hold there.

The actual band projections give a direct counterexample. Let
$P_b:\mathbb F_2^{64}\to\mathbb F_2^{c_b}$ be the three surjective band
maps, with $(c_0,c_1,c_2)=(42,43,43)$. Choose a nonzero $x_0$ outside all
three kernels and define

\[
f_b(y):=\boldsymbol 1[y=P_bx_0].
\]

The joined projections determine $x_0$. For uniform
$X\in\mathbb F_2^{64}\setminus\{0\}$,

\[
\mathbb E\!\left[\prod_b f_b(P_bX)\right]=\frac1{2^{64}-1},
\qquad
\mathbb E[f_b(P_bX)^{1/p_b}]
=\frac{2^{64-c_b}}{2^{64}-1}.
\]

The admissible linear-BL point $(p_0,p_1,p_2)=(1/2,1/2,1)$ would assert

\[
\frac1{2^{64}-1}
\le
\frac{2^{42.5}}{(2^{64}-1)^2},
\]

which is false. Thus the conditional-BL premise fails even though every pair
of coefficients sums to at least one.

If the conditional inequality were available, the later algebra would be
correct. Writing $P:=\sum_bp_b$ and $q_b:=p_b/P$ gives
$p_b/q_b=P$. Hölder across bands would therefore place
$\overline m_{b,r}^{P}\exp(-\eta{r\choose2}/q_b)$ in the band polynomial and
raise its coefficient to $q_b$. Real powers
$h_b=1/p_b\ge1$ also cause no analytic problem: inclusion-exclusion remains
an identity for the fixed monomial, and Jensen remains valid.

That correct downstream algebra does not repair the false premise. The
original support-conditioned argument may use ordinary fractional Hölder
weights whose sum is at most one. A genuine linear-BL improvement must be
applied before support conditioning, or it must include a newly proved
constant for the punctured measure. Neither route currently preserves the
pair-capacity support sum in the proposed form.

## Motif convergence gate

The test $e\Delta\max_b\rho_b<1$ is not, by itself, a complete polymer
certificate for the expansion after (7). It bounds an edge-animal branching
series. A polymer proof must also control incompatibility between connected
edge sets. In the uncentered expansion, the polymer activities additionally
depend on the one-vertex fugacity selected by coefficient extraction.

A fugacity-independent formulation first converts the binary occupation
variables to Ising spins. For an edge of color $b$, set

\[
J_b:=\log(1+\rho_b),
\qquad
\tau_b:=\tanh(J_b/4).
\]

The constant and one-spin terms can be absorbed into the vertex weights.
The remaining high-temperature expansion assigns each edge the activity
$\tau_b$. Summing the spins gives factors of absolute value at most one.
Thus every connected edge polymer $\gamma$ has

\[
|w(\gamma)|\le\prod_{e\in E(\gamma)}\tau_{c(e)}.
\]

The following Kotecký--Preiss condition is a sound convergence gate:

\[
\sup_v\sum_{\gamma\ni v}
|w(\gamma)|e^{|V(\gamma)|}\le1.
\tag{A3}
\]

Condition (A3) controls both connected-polymer weights and their overlap
incompatibility. It also supplies an absolute tail bound after a finite motif
truncation.

For a scalar gate, let $\Delta$ be the maximum degree, let
$D:=2(\Delta-1)$ be the maximum degree of the line graph, and let
$\tau:=\max_b\tau_b$. The standard connected-subset bound in the line graph
gives at most $\Delta(eD)^{k-1}$ rooted polymers with $k$ edges. Therefore

\[
e^2(\Delta+D)\tau\le1
\tag{A4}
\]

is a simple sufficient condition for (A3). This condition is conservative;
weighted colored counting can improve it.

For the fitted values, $e\Delta\max_b\rho_b\approx9.6$, so even the proposed
weaker test fails. After centering,
$\max_b\tau_b\approx0.00461$. With $\Delta=189$ and $D=376$, the left side
of (A4) is about $19.3$. The current motif branch must therefore stop under a
sound elementary polymer gate. This verdict does not exclude a sharper
colored cluster expansion. Dobrushin-type uniqueness alone does not certify
the six-edge truncation and its remainder.
