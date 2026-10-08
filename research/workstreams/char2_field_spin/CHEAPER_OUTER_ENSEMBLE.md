# Cheaper outer ensembles: necessary gates and exact local laws

The current higher-rate target uses a GF256 outer of length 256 and dimension 176, followed by one coordinate-regional route. The inner has width 16 and state dimension 8. Replacing the large outer by independent small blocks fails deterministically. Two inexpensive column-dependent mixers also fail, for different reasons. A two-shear butterfly remains a candidate, but it has no whole-code certificate.

This note records proof screens, not implementation timings. The probability statements below concern explicitly stated ideal setup distributions.

## A supported subcode forces a window obstruction

Fix an outer code \(C\subseteq\mathbb F^{256}\). Suppose \(C\) contains an \(r\)-dimensional subcode supported in \(h\) coordinate regions. Apply this subcode independently to \(L\) rows. Fix any realization of the regional route and of a linear causal inner with state dimension \(s\), zero initial state, and identity feedthrough.

Choose one packet-aligned contiguous window in each of the \(h\) regions. Let their total length be \(E\). Requiring all outer coordinates outside these windows to vanish imposes at most \(hL-E\) equations. The remaining source dimension is therefore at least

\[
D\ge rL-(hL-E)=E-(h-r)L.
\]

Next impose zero state at every window end. These conditions impose at most \(sh\) additional equations. If

\[
E-(h-r)L>sh,
\tag{1}
\]

a nonzero source remains. Zero input and zero state force zero output throughout every gap. The output is nonzero, because identity feedthrough makes the causal inner invertible. Its weight is at most \(E\).

Thus (1) is a deterministic failure criterion. It holds for every route, every allowed window placement, and every inner coefficient realization. Input diagonals do not change this argument because they preserve zero coordinates.

At requested \(K=2^{20}\), the current padded geometry has

\[
L=6016,\quad N=256L=1540096,\quad e=154009,\quad s=8.
\]

A nonzero output of weight at most \(e\) violates the strict 10% target. The smallest aligned budget used here is

\[
E=16\left\lceil\frac{(h-r)L+sh+1}{16}\right\rceil.
\]

| Supported subcode \([h,r]\) | Aligned \(E\) | Dimension lower bound | Exit equations | Result |
|---|---:|---:|---:|---|
| \([16,11]\) | 30224 | 144 | 128 | Deterministic failure |
| \([32,22]\) | 60432 | 272 | 256 | Deterministic failure |
| \([64,44]\) | 120848 | 528 | 512 | Deterministic failure |
| \([128,88]\) | 241680 | 1040 | 1024 | Budget exceeds \(e\); not excluded |
| \([26,1]\) | 150624 | 224 | 208 | Deterministic failure |
| \([27,1]\) | 156640 | 224 | 216 | Budget exceeds \(e\); not excluded |

These failures do not rely on padding. Discard the final 59 complete source rows, leaving \(176(6016-59)=1048432\le2^{20}\) source coordinates. The dimension bound becomes

\[
D\ge E-(h-r)L-59r.
\]

The four failing cases then require budgets 30864, 61728, 123440, and 150672. All remain below \(e\). In the weight-26 case, \(D\ge213>208\).

Consequently, a shared outer with a word of weight at most 26 cannot meet this geometry's distance target. Outer distance at least 27 is necessary, not sufficient. Independent small-block outers fail even when every constituent is MDS. An invertible source premix cannot repair them: it changes the parametrization, not the outer code.

## Sparse systematic parity is already expensive at the first gate

Consider a systematic generator \([I;P]\), with 176 source coordinates and 80 parity coordinates. If a column of \(P\) has at most 25 nonzeros, its source basis word has weight at most 26. The preceding obstruction applies.

A fixed-degree design therefore needs degree at least 26 in every column. A direct GF256 implementation already uses \(176\cdot26=4576\) scalar-coefficient products per four outer rows, before any sharing optimization.

For an independent-entry ensemble, let each parity entry be nonzero with probability \(\rho\). The exact probability that some source basis word triggers the obstruction is

\[
1-\left(1-\sum_{j=0}^{25}\binom{80}{j}\rho^j(1-\rho)^{80-j}\right)^{176}.
\tag{2}
\]

At \(\rho=0.73\), (2) is about \(2^{-39.2272}>2^{-40}\). At \(\rho=0.74\), this particular lower bound falls below \(2^{-40}\). Passing it proves nothing about other words. Thus independently sparse entries cannot remain sparse while meeting even this necessary setup margin.

There is an exact rank-moment interface for this ensemble. If \(u\) systematic coordinates are allowed and \(v\) parity coordinates are forbidden, write \(D=\dim\ker P_{v\times u}\). Independent sparse-uniform entries give

\[
\mathbb E[q^D]
=\sum_{a=0}^{u}\binom ua(q-1)^a
\left[\frac1q+\frac{q-1}{q}
\left(1-\frac{q\rho}{q-1}\right)^a\right]^v.
\tag{3}
\]

Equation (3) is useful for a rank-aware proof; it is not a whole-code bound. A setup matrix shared across outer rows does not make their shortened dimensions independent.

## Column-dependent mixing changes the code

Start with eight copies of GRS\([32,22]\). At each coordinate position \(j\), collect the eight component symbols in \(x_j\in\mathbb F_{256}^8\) and emit \(M_jx_j\). The resulting code has length 256 and dimension 176 when every \(M_j\) is invertible.

Using the same \(M_j=M\) at every position merely recombines component codewords. The direct-sum obstruction remains. Coordinate-dependent mixing is genuinely different. It uses fixed-size local maps before the existing regional route, not another global routing stage.

One layer of two-by-two switches still decomposes the code into four supported \([64,44]\) subcodes. It fails (1), regardless of the switch coefficients. Two layers are the first connectivity pattern not excluded by that argument.

A particularly cheap three-stage switch also fails. Consider

\[
(x,y)\longmapsto(u,v)=(x+y,\ y+a(x+y)).
\tag{4}
\]

Use (4) on the standard eight-wire butterfly, with distances 1, 2, and 4. Every coefficient may vary independently with the position and switch. Nevertheless,

\[
(z,z,z,z,z,z,z,z)\longmapsto(0,0,0,0,0,0,0,z)
\]

for every coefficient assignment. Taking eight identical constituent codewords retains a supported \([32,22]\) subcode. Random coefficients do not repair this common invariant subspace.

## A fixed binary postmix with fresh diagonals fails the 40-bit target

Another inexpensive proposal sets

\[
M_j=(I+J_8)D_j,
\]

where \(J_8\) is the all-ones matrix and every diagonal entry of every \(D_j\) is sampled independently from \(\mathbb F_{256}^{*}\). The matrices are sampled once and shared across the \(L\) rows. The binary matrix is invertible: \((I+J_8)^2=I\).

The following lower bound uses only one fixed pair of components. Let \(S\) be a 13-element subset of the 32 evaluation points, and let

\[
g_S(X)=\prod_{a\notin S}(X-a).
\]

There are \(m=(q^2-q)/2=32640\) monic irreducible quadratics over \(\mathbb F_q\), for \(q=256\). For every ordered distinct pair \((f,g)\), the words \(g_Sf\) and \(g_Sg\) have degree 21 and support exactly \(S\). The ratio \(f/g\) is well-defined at every evaluation point.

Let \(E_{S,f/g}\) be the event that the relative diagonal entries on \(S\) equal \(\alpha f/g\) for some \(\alpha\in\mathbb F_q^*\). It has probability \((q-1)^{-12}\). On this event, choose proportional component words so that their diagonally scaled values agree. Their sum is zero. The map \(I+J_8\) therefore leaves these two components unchanged and produces an outer word of weight exactly 26.

The window argument then yields whole-code failure for every route and every inner coefficient realization. It also applies after discarding the 59 padded rows.

Here is an exact Bonferroni bound on the setup event. Write \(R=m(m-1)\) and

\[
\lambda=\binom{32}{13}R\,255^{-12}.
\]

For two supports with intersection size \(j\), the number of ordered support pairs is

\[
c_j=\binom{32}{13}\binom{13}{j}\binom{19}{13-j}.
\]

For the same ratio and distinct supports, the joint probability is

\[
p_j=\begin{cases}255^{-24}&j=0,\\255^{-(25-j)}&j\ge1.\end{cases}
\]

For distinct ratios, this is an upper bound when \(j\le4\). Their intersection is empty when \(j\ge5\): equality up to a scalar at five points forces a degree-at-most-four cross-product polynomial to vanish identically. Distinct monic irreducible pairs give distinct ratios, even up to a scalar.

It follows that

\[
\Pr[\text{bad setup}]
\ge\lambda-\frac12\sum_{j=0}^{12}c_jp_j
\left(R+\mathbf1_{j\le4}R(R-1)\right)
=1.979890647388062\ldots\times10^{-12}
>2^{-40}.
\tag{5}
\]

The displayed decimal is only a presentation value. The replay compares exact rational numbers. The lower bound is about \(2^{-38.8777}\); its pair term is at most \(0.595574\lambda\). Thus even this fully connected cheap mixer does not meet the requested margin under the stated shared setup law.

## A two-shear butterfly and a full-GL reference

A two-product switch avoids the invariant in (4):

\[
u=x+b y,\qquad v=y+a u,
\qquad a,b\gets\mathbb F_q^*.
\]

Its support law does not depend on the nonzero input values. Encode the first and second wire by bits 1 and 2, and put \(r=1/(q-1)\). The exact transition probabilities are

| Input mask | Output probabilities |
|---|---|
| 0 | \(0:1\) |
| 1 | \(3:1\) |
| 2 | \(1:r,\ 3:1-r\) |
| 3 | \(2:r,\ 1:r(1-r),\ 3:(1-r)^2\) |

Composing twelve independent switches gives an exact 256-state support kernel for one eight-wire column. The implementation exposes this finite kernel. It is not yet a shortened-rank bound across 32 columns or across shared outer rows.

For comparison, let the \(M_j\) instead be independent uniform elements of \(\mathrm{GL}(8,q)\). Regard eight constituent words as one GRS\([32,22]\) word over an alphabet of size \(Q=q^8\). For a fixed support of size \(u\), its number of nonzero words is

\[
A_u(Q)=(Q-1)\sum_{i=0}^{u-11}(-1)^i\binom{u-1}{i}Q^{u-11-i},
\qquad u\ge11.
\]

Every nonzero column becomes uniform among the \(Q-1\) nonzero eight-vectors. Therefore the expected scalar weight enumerator is

\[
\sum_{u=11}^{32}\binom{32}{u}A_u(Q)
\left(\frac{(1+(q-1)z)^8-1}{Q-1}\right)^u.
\tag{6}
\]

Sum coefficients of degrees at most 26 and divide by \(q-1\). Exact arithmetic bounds the probability of outer distance at most 26 by \(2^{-322}\). This passes only the necessary outer-distance gate, not the whole-code target.

There is also an exact kernel-moment interface. Suppose a support permits \(t_j\) of the eight output coordinates at column \(j\). Put \(w_j=(q^{t_j}-1)/(Q-1)\), and let \(e_u\) denote the elementary symmetric polynomial. For the shortened outer dimension \(D\),

\[
\mathbb E[q^D]=1+\sum_{u=11}^{32}A_u(Q)e_u(w_1,\ldots,w_{32}).
\tag{7}
\]

Equation (7) is for one outer row. Raising its expectation to the number of rows would incorrectly treat the shared matrices as fresh randomness. A usable replacement must control this common setup or explicitly change the setup law.

## Cost and next action

The following counts are algebraic GFNI512 product counts per four outer rows. They omit Boolean operations, movement, and systematic normalization. Eight existing small GRS kernels contribute approximately \(8\cdot41=328\) products; see [PARITY_STATE_COST.md](PARITY_STATE_COST.md).

| Column postmix | Additional products | With eight small outers | Current result |
|---|---:|---:|---|
| Three-stage one-product switch | 384 | 712 | Deterministic failure |
| \((I+J_8)D_j\) | 256 | 584 | Failure probability exceeds \(2^{-40}\) |
| Three-stage two-shear switch | 768 | 1096 | Exact local support law; no global bound |
| Full \(\mathrm{GL}(8,q)\) | 2048 | 2376 | Strong outer-distance gate; expensive |

For context, the current large-outer FFT transpose uses 1203 products per four rows. Its lower-product basis circuit has different Boolean and data-movement costs; product counts alone do not predict speed. The column-postmixed generator is not systematic in its original source coordinates. Restoring a systematic interface can require additional inverse column maps. The table is therefore optimistic for that interface.

The next useful proof task is the two-shear butterfly's joint shortened-rank spectrum under a setup shared across rows. Its local support kernel is available, but neither an ordinary word first moment nor an independent-row substitution suffices. Do not implement the rejected one-product or fixed-binary-diagonal alternatives. No cheaper outer in this note has a whole-code certificate, and the existing certified GRS256/176 construction remains unchanged.

Replay:

```text
python -B research/workstreams/char2_field_spin/theory/cheaper_outer_ensemble_tests.py
python -B research/workstreams/char2_field_spin/theory/cheaper_outer_ensemble_screen.py
```

Eight tests pass. They include exhaustive GF4 switch and sparse-rank checks, an independent tiny causal confinement example, and an exhaustive GF4 check of the Bonferroni event multiplicities. The production-size comparisons use exact integers and fractions. No native build or benchmark was run.
