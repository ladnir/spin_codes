# Exact dense-class certificate for GRS[128,87]

The K20 dense class passes with at least 13,690 bits of margin. This is a separate proof for the proposed GRS[128,87] outer. It does not alter the GRS[256,176] certificate or establish the new construction's whole-code bound.

## Construction and integer scope

Work over GF256. The outer consists of \(L\) rows of a fixed GRS[128,87] code. Each coordinate region independently permutes its row labels. Independent nonzero coordinate diagonals are sampled independently of this route and of the inner.

The intended inner is the unchanged symmetric Cauchy recurrence with width 16, state dimension 8, eta 2, and fresh nonzero state diagonals per packet. It starts with zero state and discards its final state. Identity feedthrough makes the inner invertible. The dense proof only needs invertibility and independence from the route and input diagonals.

Pad the requested \(K=2^{20}\) to a multiple of \(87\cdot64\). This gives

\[
L=12096,\qquad K_{\mathrm{pad}}=1052352,\qquad
N=128L=1548288.
\]

The design rate is \(87/128>2/3\). The forbidden output weights are 1 through

\[
e=\lceil N/10\rceil-1=154828.
\]

For an output \(y\), let \(x=B^{-1}y\), where \(B\) is the inner map. This certificate covers exactly the dense inverse-input class

\[
371590=\lceil0.24N\rceil\le\operatorname{wt}(x)\le N.
\]

## Membership bound and exact length lift

Condition on \(B\) and fix an inverse input \(x\) of density \(\alpha\). The regional membership theorem in [REGIONAL_GLOBAL_ROUTE.md](REGIONAL_GLOBAL_ROUTE.md) gives

\[
\Pr[x\in PD(C^{\oplus L})]
\le (L+1)^{128}
\left[
\lambda^{-128\alpha}
\left((1-\alpha)^{128}
+256^{-41}(1-\alpha+\tfrac{256}{255}\lambda\alpha)^{128}\right)
\right]^L
\tag{1}
\]

for every \(\lambda>0\). The bound is uniform over the vector's regional weight profile. Conditioning on \(B\) preserves the route and input-diagonal laws. There is no union over regional profiles.

The old GRS[32,22] rational witnesses can be reused exactly. For a fixed density and tilt, define

\[
X=(1-\alpha)^{32},\qquad
Y=256^{-10}(1-\alpha+\tfrac{256}{255}\lambda\alpha)^{32}.
\]

The GRS[128,87] bracket is

\[
X^4+\frac{Y^4}{256}\le(X+Y)^4.
\tag{2}
\]

The additional parity coordinate strengthens the fourth-power comparison. This is a scalar inequality between MDS membership bounds; it does not assert a block or tensor structure for the new outer.

The replay freezes 152 rational tilts, all at least one. For cell \([j/200,(j+1)/200]\), it bounds the decreasing first term at the left endpoint and the increasing second term at the right endpoint. It bounds \(\lambda^{-128\alpha}\) at the left endpoint. Each cell is accepted only after checking the small witness, (2), and the target inequality using exact fractions:

\[
\operatorname{base}_{128,87}^{200}\,2^{256\cdot128}
\le\lambda^{128j}.
\tag{3}
\]

Equation (3) proves a per-coordinate membership exponent of \(32/25\). The 152 cells cover every real density in \([0.24,1]\). The replay also checks a disjoint integer partition of weights 371590 through 1548288, including both endpoints.

## From membership to bad-output existence

The GF256 ball through radius \(e\) contains at most \(2^{1.27N}\) vectors. The exact Chernoff comparison is

\[
\sum_{j=0}^{e}\binom Nj255^j
\le(10/9)^N2295^{N/10},
\qquad (10/9)^{100}2295^{10}<2^{127}.
\]

Combine this with (1)--(3). For every fixed independent invertible inner, the probability of any forbidden output in the dense class is at most

\[
(L+1)^{128}2^{-N/100}
\le2^{-(\lfloor N/100\rfloor-128\cdot14)}
=2^{-13690}.
\]

The output union already covers all dense inverse weights. There is no factor of 152: those cells establish one uniform membership bound, rather than 152 independently summed bad events.

The [new middle replay](FIELD128_DENSITY_CERTIFICATE.md) ends at weight 371589, so its upper boundary joins this dense class exactly. Sparse coverage and the final sum must be checked separately. For that join, the new outer's minimum nonzero row weight is 42; hence 1106 active rows imply weight at least 46452, above the middle lower endpoint 46449.

This partial certificate is over GF256. A complete base-field distance theorem can later lift to compatible extension fields. The dense class alone does not lift: a nonzero base-field component of an extension-field word can have smaller inverse-input weight.

## Replay

```text
python -B research/workstreams/char2_field_spin/theory/field128_dense_certificate_tests.py
python -B research/workstreams/char2_field_spin/theory/field128_dense_certificate.py
```

Four tests pass, including all 152 exact cell checks, integer coverage, endpoint inequalities, the extra-parity identity, and rejection of an invalid tilt. The independent CLI replay reports margin 13690. No floating optimizer, native encoder, benchmark, or remote execution is used. Next, join matching sparse and middle receipts and sum their rigorous failure bounds before making a whole-code claim.
