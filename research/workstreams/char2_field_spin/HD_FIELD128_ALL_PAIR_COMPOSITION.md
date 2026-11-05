# Composing every packet pair before the fractional power

The gamma full-RS construction admits composition before the fractional power for every consecutive packet pair. The earlier restriction to at most one singleton per packet reduced coefficient cost; it was not a hypothesis of the probabilistic argument. This note justifies the stronger bound without changing the encoder.

Fix the construction and local caps from [the gamma return proof](HD_FIELD128_GAMMA_RETURN.md). In particular, work over \(F=\operatorname{GF}(2^{128})\), use width 64 and state dimension \(s\in\{16,24,28,30\}\), and fix a degree-16 scalar \(\gamma\) over GF256. All outer maps, regional permutations, pair swaps, and packet refreshes have the stated independent distributions.

## From the outer comparison to a fixed-route moment

Fix a nonempty set \(U\) of \(r\) active outer groups, a mark \(0<z<1\), and an exponent \(0<\alpha\le1\). Put \(q=|F|\), and let \(\beta(q)=q^{160}/(q-1)^{32}\) be the per-active-group factor from [the outer counting domination](HD_FIELD128_LIFT_PROOF.md#outer-counting-measure). Let \(L\) be the forbidden output weight. Write \(\Pi\) for the complete route, including swaps, and \(D\) for all refresh diagonals. Define

\[
Z_U(\Pi,D)
:=\mathbb E_X\bigl[z^{\operatorname{wt}(I_{\Pi,D}X)}\bigr],
\qquad
X\gets F^{256r},
\]

where the remaining outer coordinates are zero and \(I_{\Pi,D}\) is the routed inner encoder. This \(X\) belongs to the outer comparison measure; it is not asserted to be an actual outer codeword.

Let \(B_U\) be the event that a nonzero message with active group set exactly \(U\) has output weight at most \(L\). Conditional on \(\Pi,D\), the outer counting domination and projective scaling give

\[
\Pr_{\rm outer}[B_U\mid\Pi,D]
\le\min\{1,a_r Z_U(\Pi,D)\},
\qquad
a_r=\frac{\beta(q)^r}{(q-1)z^L}.
\]

The comparison may include all-zero input values. No such mass is removed here. Since \(\min(1,x)\le x^\alpha\) for \(x\ge0\), conditional Jensen gives

\[
\Pr[B_U]
\le a_r^\alpha\,
\mathbb E_\Pi\left[
\left(\mathbb E_D[Z_U(\Pi,D)\mid\Pi]\right)^\alpha
\right].
\tag{1}
\]

This Jensen step averages the refresh randomness after fixing the whole route. It does not average route histograms or swaps inside the fractional power.

## Why the raw local operators compose

Fix \(\Pi\). Its packet \(t\) has a fixed potential input set \(S_t\) and occupancy \(j_t=|S_t|\). Under the product comparison measure, the input values in packet \(t\) are uniform on \(F^{S_t}\). They are independent of all previous input values and refresh coins. Thus they are independent of the incoming state, including the measure weighted by preceding output marks.

The local proof uses a cone of nonnegative measures on the state space. A measure is represented by \((Z,W)\) when its mass at zero is at most \(Z\) and each nonzero atom is at most \(W/(q^s-1)\). The gamma local operator \(P_j(z)\) maps such an incoming bound to

\[
(Z,W)P_j(z).
\]

This domination includes the current packet's output mark and fresh refresh. It holds for every position set of size \(j\), not an average over position sets. The MDS rank formulas and deterministic principal-nonsingularity lemma establish this uniformity for every \(j=0,\ldots,64\).

Applying the domination successively therefore yields

\[
\mathbb E_D[Z_U(\Pi,D)\mid\Pi]
\le e_0\left(\prod_{t=1}^{32768}P_{j_t}(z)\right)\mathbf1,
\qquad e_0=(1,0),\quad \mathbf1=(1,1)^\mathsf T.
\tag{2}
\]

The terminal vector bounds total mass in both state classes. No state reset or independence of consecutive occupancies is assumed.

Group the factors in (2) into chronological pairs,
\(B_i=P_{j_{2i-1}}P_{j_{2i}}\). Expanding only the state classes at pair boundaries and applying fractional subadditivity gives

\[
\left(e_0\prod_iB_i\mathbf1\right)^\alpha
\le e_0\prod_i B_i^{\circ\alpha}\mathbf1.
\tag{3}
\]

Equations (1)--(3) justify every paired factor \((P_jP_k)^{\circ\alpha}\). There is no singleton-count restriction. In particular, the argument does not infer validity merely because a new coefficient table is smaller than an old upper bound.

## Averaging the route after composition

A packet contains 32 routed pairs. Let \(h=(u,v)\) count its singleton and full pairs, respectively. There are

\[
m(h)=\binom{32}{u}\binom{32-u}{v},\qquad u+v\le32,
\]

assignments of these pair types to its slots. The local cap is \(P_{u+2v}\) for every singleton orientation. Consequently the fair-swap probabilities sum to one without a Jensen step. The factors \(2^u\) counting different source supports, when applicable, remain in the separate support count.

Define the complete paired polynomial

\[
T(x,y)=\sum_{h=(u,v)}\sum_{k=(u',v')}
m(h)m(k)x^{u+u'}y^{v+v'}
\bigl(P_{u+2v}P_{u'+2v'}\bigr)^{\circ\alpha}.
\tag{4}
\]

The order of the matrix product in (4) is chronological. Multiplicities lie outside the fractional power. All \(561^2\) ordered histogram pairs occur, although only \(65^2\) distinct occupancy products are needed.

For a regional histogram with \(n_1\) singleton and \(n_2\) full pairs, the averaged regional cap is

\[
R(n_1,n_2)=
\frac{[x^{n_1}y^{n_2}]T(x,y)^{64}}
{\binom{4096}{n_1}\binom{4096-n_1}{n_2}}.
\tag{5}
\]

Coefficient extraction retains sampling without replacement inside the region. It does not replace the packet sequence by independent draws. Independent regional routes allow the usual ordered even/odd product across 256 regions, with the state carried through every boundary. The union over \(U\) remains outside the fractional power.

For nonnegative matrices, \((AB)^{\circ\alpha}\le A^{\circ\alpha}B^{\circ\alpha}\). Thus (4) is coefficientwise no larger than the earlier good-pair polynomial, with equality at \(\alpha=1\). This comparison measures the improvement after validity has been established by (1)--(3).

## Implementation scope

The [all-pair diagnostic](theory/hd_field128_gamma_all_pair_screen.py) stores the complete paired polynomial, with maximum singleton degree 64, full-pair degree 64, and input degree 128. Its adapter sets the single-packet components to zero and uses the existing pair power 64 and regional normalizer. It does not discard the former dirty cases.

The implementation uses FLOAT arithmetic and remains a proposal diagnostic. Neither coefficientwise improvement nor a negative FLOAT result supplies a directed certificate. A directed implementation and complete support coverage remain separate obligations. The original [good-pair note](HD_SHIFTED_GOOD_PAIRS.md) and all frozen gamma receipts are unchanged.
