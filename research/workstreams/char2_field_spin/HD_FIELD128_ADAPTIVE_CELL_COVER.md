# Adaptive profile cells for the all-pair bound

The selected-weight verifier evaluates every even/odd profile separately. Repeating that quadratic calculation at every weight is unnecessary when nearby profiles admit a common upper bound. This note gives such a bound using directed alias numerators, convex scalar normalization, and exact projected support counts. It retains the chronological state transfer and uses no common state potential.

This is a coverage lemma, not a completed cover or a claim about unsampled weights. No coefficient provider or frozen receipt is changed. A new implementation requires independent review and directed arithmetic.

## Fixed construction and profile domain

Fix one gamma full-RS construction from [the gamma return note](HD_FIELD128_GAMMA_RETURN.md), including its state dimension. Use the [all-pair composition](HD_FIELD128_ALL_PAIR_COMPOSITION.md). Put

\[
q=2^{128},\quad P=4096,\quad G=2P=8192,\quad n=2^{21},
\quad \beta(q)=\frac{q^{160}}{(q-1)^{32}}.
\]

For an active source-group set \(U\subseteq\{0,\ldots,G-1\}\), let \(r=|U|\), and let \(a,b\) count its full pairs under the even and odd cyclic matchings. For \(0<r<G\), the number of active runs in the group cycle is

\[
k=r-a-b\ge1.
\]

Every such integer profile belongs to the real polytope

\[
\mathcal D=\left\{(r,a,b):
\begin{array}{l}
1\le r\le G-1,\quad 0\le a,b\le P,\\
2a,2b\le r,\quad a,b\ge r-P,\quad r-a-b\ge1
\end{array}\right\}.
\tag{1}
\]

We only require this containment. We do not assert that every integer point in \(\mathcal D\) is realized by a support. Impossible profiles can be retained as a conservative relaxation. The all-active profile \((G,P,P)\) is handled separately or by a terminal certificate; the empty support is excluded.

A cell \(B\) is a bounded rational polytope contained in \(\mathcal D\). It represents the supports whose integer profiles lie in \(B\). Each cell has fixed exact witnesses \(z,\alpha\), where \(0<z<1\) and \(0<\alpha\le1\), and a fixed forbidden weight \(L=\lfloor\delta n\rfloor\). Different cells may use different witnesses for the same construction and target.

## Regional matrices over one cell

Let \(T(x,y)\) be the positive paired polynomial before regional extraction, with the fractional power already placed as in the all-pair proof. For a regional weight \(r\) and full-pair count \(f\), define

\[
u=r-2f,\qquad e=P-r+f,\qquad
H(r,f)=\frac{P!}{e!u!f!},\qquad
R(r,f)=\frac{[x^uy^f]T(x,y)^{64}}{H(r,f)}.
\tag{2}
\]

The entries of \(R\) are regional caps, not transition probabilities. They already average the route after the all-pair fractional bound. For a fixed support of profile \((r,a,b)\), its failure probability is at most

\[
\left(\frac{\beta(q)^r}{(q-1)z^L}\right)^\alpha
e_0\bigl(R(r,a)R(r,b)\bigr)^{128}\mathbf1,
\quad e_0=(1,0),\quad\mathbf1=(1,1)^\mathsf T.
\tag{3}
\]

Choose one fixed positive marker pair \((x_i,y_i)\) for each parity \(i=0,1\) within the cell. Let \(N_i(B)\) be a nonnegative matrix that bounds, entrywise, every marked regional numerator

\[
x_i^{r-2f_i}y_i^{f_i}
[x^{r-2f_i}y^{f_i}]T(x,y)^{64},
\qquad (r,a,b)\in B\cap\mathbb Z^3,
\quad f_0=a,\ f_1=b.
\tag{4}
\]

The superscripted marker factor in (4) multiplies the coefficient; it is not part of the coefficient-extraction operator. A directed alias atlas supplies \(N_i(B)\): take componentwise maxima over all residues in an enclosing histogram rectangle. The rectangle may include additional histograms. Positive alias sums still dominate every requested coefficient. A directed positive-marker fallback may be combined by taking componentwise minima of independently valid upper bounds.

Set

\[
D_i(r,f)=H(r,f)x_i^{r-2f}y_i^f,
\qquad
Q_B=e_0\bigl(N_0(B)N_1(B)\bigr)^{128}\mathbf1.
\tag{5}
\]

Monotonicity of products of nonnegative matrices gives a bound on (3) by

\[
Q_B\,2^{\Phi(r,a,b)},
\]

where

\[
\Phi(r,a,b)=
\alpha\bigl(r\log_2\beta(q)-\log_2(q-1)-L\log_2z\bigr)
-128\bigl(\log_2D_0(r,a)+\log_2D_1(r,b)\bigr).
\tag{6}
\]

The matrix multiplication order is even then odd. The state is carried through all 256 regions. No commutation or state reset occurs.

## A vertex bound for the scalar factor

Extend \(H\) to real nonnegative histogram counts using

\[
H(r,f)=\frac{\Gamma(P+1)}{\Gamma(P-r+f+1)\Gamma(r-2f+1)\Gamma(f+1)}.
\tag{7}
\]

This agrees with (2) at integer histograms. The function \(\log\Gamma(t+1)\) is convex for \(t\ge0\). Thus \(-\log H\) is a sum of convex functions of affine counts, plus a constant. All remaining terms of (6) are affine. Consequently \(\Phi\) is convex on \(\mathcal D\), including lower-dimensional faces.

Let \(V(B)\) be the complete vertex set of the nonempty cell. Every point of \(B\) is a convex combination of its vertices, so

\[
\Phi(r,a,b)\le v_B:=\max_{v\in V(B)}\Phi(v)
\quad\text{for every }(r,a,b)\in B.
\tag{8}
\]

A verifier must evaluate (8) with directed upper rounding. Rational cell vertices need not be integer profiles. Their histogram counts are still nonnegative by (1), so (7) is defined. Endpoint checks along the weight axis alone are not sufficient; all cell vertices are required.

## Counting supports in a cell without enumerating profile pairs

Let \(C_B\) be any upper bound on the number of actual supports represented by \(B\). Summing (3) over those supports and applying (8) gives

\[
\boxed{\qquad
\Pr[\text{a bad message has its active support in }B]
\le C_B Q_B 2^{v_B}.
\qquad}
\tag{9}
\]

There is no fractional power on \(C_B\).

Useful bounds for \(C_B\) come from projections. The exact number of supports with weight \(r\) and full-pair count \(f\) in one fixed matching is

\[
M_P(r,f)=\binom Pf\binom{P-f}{r-2f}2^{r-2f}.
\tag{10}
\]

The exact number with weight \(r\) and \(k\) active runs in the group cycle is

\[
A_G(r,k)=\frac Gk\binom{r-1}{k-1}\binom{G-r-1}{k-1},
\quad 1\le k\le\min(r,G-r).
\tag{11}
\]

Set both counts to zero outside their feasible integer ranges. For each integer \(r\) in the cell's weight projection, choose integer sets \(E_r,O_r,K_r\) containing, respectively, all projected values of \(a,b,r-a-b\) from integer profiles in the cell. Integer interval envelopes of the corresponding real projections suffice. A valid choice in (9) is

\[
C_B:=\min\left\{
\sum_r\sum_{a\in E_r}M_P(r,a),\quad
\sum_r\sum_{b\in O_r}M_P(r,b),\quad
\sum_r\sum_{k\in K_r}A_G(r,k)
\right\}.
\tag{12}
\]

To justify the first sum, every support in the cell appears once in its even-matching marginal class. The other sums use the odd-matching and run partitions. They may include supports outside the cell, which only enlarges the bound. No independence between the matching counts is assumed. There is no extra factor for the number of values of the omitted profile coordinate.

For each fixed weight, exact marginal and run arrays have length \(O(G)\). Prefix sums answer interval queries without an \(O(G^2)\) profile-pair loop. A bounded implementation can stream weights and accumulate counts for a batch of cells instead of retaining every weight's arrays. This is an arithmetic operation-count observation, not a runtime measurement.

## Preserve the singleton-free faces

Before adaptive subdivision, partition the integer profile domain into

\[
\begin{aligned}
\mathcal F_0&=\mathcal D\cap\{2a=r\},\\
\mathcal F_1&=\mathcal D\cap\{2b=r\},\\
\mathcal I&=\mathcal D\cap\{r-2a\ge1,\ r-2b\ge1\}.
\end{aligned}
\tag{13}
\]

The faces do not intersect inside \(\mathcal D\), since simultaneous equality would give \(k=0\). Equation (13) covers all integer profiles in (1). It is a partition of integer profiles, not of every real point of \(\mathcal D\).

On \(\mathcal F_0\), write \(r=2u\), \(a=u\), and \(b=u-k\). Contract the selected even pairs to a cycle of length \(P\). The exact support count at \((u,k)\) is

\[
A_P(u,k)=\frac Pk\binom{u-1}{k-1}\binom{P-u-1}{k-1},
\quad 1\le u\le P-1,\quad 1\le k\le\min(u,P-u).
\tag{14}
\]

For a cell on that face, sum (14) over its integer \((u,k)\) domain to obtain the exact number of represented supports. At fixed integer \(u\), a convex cell gives an interval of feasible \(k\); a prefix sum evaluates it. The same count applies on \(\mathcal F_1\), with \(a=u-k\) and \(b=u\). The chronological matrix factors must still be kept in their respective orders.

Convexity of (6) is preserved under these affine face parameterizations. Face vertices, alias numerator bounds, and (9) therefore apply without changing the transfer argument. Keeping the faces separate retains the contracted-cycle count at profiles such as \((r,0,r/2)\); replacing it by a general marginal count would lose that improvement.

## Coverage and receipt obligations

A passing cell covers its complete integer domain, not merely its sampled center. Cells that fail a requested probability budget may be subdivided; they remain unresolved until their children pass. A passing selected weight does not guarantee that any nontrivial neighboring cell passes.

A complete certificate must verify the following obligations.

1. Every receipt uses the same field, fixed maps, gamma construction, state dimension, outer ensemble, route, and target. Witnesses and valid upper-bound refinements may differ by cell; constructions may not be mixed.
2. Each alias bank has directed coefficients and transforms, immutable source pins, and the full all-pair degree range. A FLOAT proposal or source hash alone does not establish a coefficient enclosure.
3. Each cell's marker pair is fixed per parity. Its numerator matrices cover every relevant histogram, its complete rational vertex set is checked, and all scalar and matrix operations are rounded upward.
4. The integer leaf domains are disjoint and cover the requested root domains. A checked recursive split at integer \(m\), into coordinate \(\le m\) and coordinate \(\ge m+1\), proves lattice coverage without enumerating all profiles. Face splits use their integer coordinates \((u,k)\). Discarded empty cells need a valid emptiness justification.
5. The root domains cover every required nonzero weight. Separately certified sparse ranges and the terminal event may replace those roots only when their weight domains and construction obligations match. The all-active support is included exactly once.
6. Sum the accepted cell ceilings and all other required event ceilings once with directed arithmetic. Events need not be probabilistically independent. A timeout leaves an explicit pending domain and cannot produce a whole-code claim.

The support-count upper bounds in (12) may count external supports conservatively; that does not alter the disjoint domains of the actual events. A final verifier must distinguish this internal counting relaxation from accidentally repeating a receipt or omitting part of the profile domain.

The next implementation step is a small directed cell checker with exact vertex and lattice-partition tests. Only after independent review should a bounded coverage driver reuse pinned all-pair banks and adaptively split unresolved cells. No such driver or completed cover is supplied by this note.
