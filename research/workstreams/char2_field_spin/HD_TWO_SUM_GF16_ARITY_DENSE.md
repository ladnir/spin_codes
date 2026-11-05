# G16 dense arities through a weight coefficient DP

The frozen Group16 native candidate uses sixteen RS[16,8] rows and the
two-sum mixer with labels `06 + GF16`. Its measured implementation is
recorded separately in [HD_TWO_SUM_GF16_NATIVE.md](HD_TWO_SUM_GF16_NATIVE.md).
This note records a proof screen, not another native experiment or a
whole-code certificate.

At one inherited witness, all dense arities 3–16 pass the FLOAT screen for
17 active source groups with ordered matching counts `(8,0)`. The cohort
has sixteen two-row groups and one dense group. The union over its 65,399
fixed dense row sets is **−585.5421219556642 bits**. Only this source class
and cohort are covered by the screen; no directed rounding is claimed.

## Fixed row sets and the counting measure

Fix the complete active-group support and one dense row set
\(A\subseteq[16]\), with \(|A|=a\). The row set remains fixed across all
16 outer columns. For an output mask \(T\), let

\[
N_A(T)=\#\{x\in\mathbb F_{256}^{A}:\operatorname{supp}(Hx)=T\}.
\]

Coordinates in this ambient comparison space may be zero. These counts
are not the conditional law of an input with every coordinate nonzero.
The established shortened-RS comparison assigns the column-mask sequence
\((T_1,\ldots,T_{16})\) mass

\[
255^{-8a}\prod_{j=1}^{16}N_A(T_j).
\]

The exact helper
[hd_two_sum_closed_kernel.py](theory/hd_two_sum_closed_kernel.py)
shows that \(\sum_{|A|=a}N_A(T)^{16}\) depends only on \(w=|T|\).
Write this integer as \(M_a(w)\), and define the per-mask weight

\[
W_a(w)=\left(\frac{M_a(w)}{255^{8a}}\right)^{1/16}.
\]

Hölder gives the pointwise comparison

\[
\sum_{|A|=a}255^{-8a}\prod_{j=1}^{16}N_A(T_j)
\leq \prod_{j=1}^{16}W_a(|T_j|).
\]

The sum on the left holds each row set fixed. The right side is a product
of positive comparison weights, not a claim that the construction samples
a new row set in each column. The multiplicity \(\binom{16}{a}\) is already
inside \(M_a\); it is not charged again.

## Ordered regional actions without a mask table

For the fixed source class, \(R_i(m)\) is the existing two-state upper
matrix for region \(i\), uniformly over its possible source subsets of
size \(m\). The matching alternates with \(i\). These matrices need not
commute, and the reduction does not assert exchangeability of the actual
subsets of the fixed source support.

Within one outer column, let \(s\) sparse groups have one active row and
\(d\) have two active rows; put \(t=s+d\). Each two-row group gives full
output support with probability \(239/255\), or a hole at each of the 16
labelled regions with probability \(1/255\). The pair choices remain fixed
across columns; this local law is identical for all 120 choices.

The suffix DP keeps a dense-weight index \(w\), a labelled-hole count
\(k\), and the two inner states. Its terminal coefficient is the common
potential \(v=(1,p)\) at \((w,k)=(0,0)\). At region \(i\), the recurrence
sums over a dense bit \(b\in\{0,1\}\) and a local hole count \(j\):

\[
V_i(w,k)=\sum_{b=0}^{1}\sum_{j=0}^{k}
 \binom{k}{j}R_i(t-j+b)V_{i+1}(w-b,k-j).
\]

Out-of-range coefficients are zero. Thus the weight coefficient sums the
actions of all labelled masks of weight \(w\), in order
\(R_0\cdots R_{15}\). It already counts each target mask once; there is no
additional factor \(\binom{16}{w}\).

For each \((s,d)\), the DP then mixes hole counts with
\(\binom{d}{k}239^{d-k}/255^d\). Denote the resulting weight-indexed vectors
by \(D_{s,d}(w)\). The empty and nonempty dense-column caps are

\[
c_0(s,d)=\max_{j\in\{0,1\}}
 \frac{W_a(0)D_{s,d}(0)_j}{v_j},\qquad
c_1(s,d)=\max_{j\in\{0,1\}}
 \frac{\sum_{w=1}^{16}W_a(w)D_{s,d}(w)_j}{v_j}.
\]

Both mixtures precede the potential quotient. The implementation uses 17
weight coefficients rather than 65,536 labelled masks. It reuses the same
DP vectors for all 14 arity classes.

## Whole-column bound and scope

For \(e\in\{0,1\}\), form

\[
P_e(X)=\sum_{s+d\leq16}
 \binom{16}{s}\binom{16-s}{d}2^s
 255^{s+2d}c_e(s,d)X^{s+2d}.
\]

Raise \(P_0(X)+YP_1(X)\) to the 16th power and retain degrees
\(X\geq288\), \(Y\geq9\). The 32 nonzero sparse RS rows each have weight
at least 9, giving the first floor. The dense group has at least one
nonzero RS row, so it occupies at least 9 columns, giving the second.

Let \(K=15708305994709246428578365440\) be the exact number of cyclic source
supports of weight 17 with ordered matching counts `(8,0)`. The retained
positive tail is multiplied by

\[
\frac{K\,17\,120^{16}}
 {255\,p\,255^{256}\,z^{734003}}.
\]

The factors respectively count source supports, the dense-group choice,
and fixed sparse pair choices. The denominator includes one projective
normalization, one endpoint-potential factor, and the 32 sparse row
normalizers. Dense row normalizers occur inside \(W_a\), not in this
prefix. The screen fixes

\[
z=\frac{8970988817957597}{9007199254740992},\qquad
\log_2p=-\frac{2266954815}{33554432}.
\]

No new witness was fitted. This is an \(\alpha=1\) comparison. It does not
reuse the older g8 certificate as a certificate for the g16 construction.

## Recorded FLOAT results and checks

| Dense arity | Failure log2 screen |
|---:|---:|
| 3 | −681.007807586 |
| 4 | −671.134474394 |
| 5 | −661.582250467 |
| 6 | −652.332001720 |
| 7 | −643.471807802 |
| 8 | −635.129666227 |
| 9 | −627.368421858 |
| 10 | −620.133755755 |
| 11 | −613.299565576 |
| 12 | −606.772533668 |
| 13 | −600.566040142 |
| 14 | −594.812802386 |
| 15 | −589.674520476 |
| 16 | −585.629347971 |
| Union | **−585.542121956** |

Every arity's largest retained term is at \((X,Y)=(288,9)\). The run took
2.776 seconds, including an independent g8 reproduction. The original
labelled g8 helper and the new weight DP both give
`−184.77272755943716` for the six-arity g8 union at the same literal.

The preceding full-arity-only control independently gave h1
`−585.6293479711762` and h2 `−476.80214895720337`, in 1.216 seconds combined.
The h2 result covers two full-arity groups only, not all dense arity pairs.
Its aggregate nonempty-column floor is 18. These results are saved in
`hd_two_sum_gf16_full_dense_8_0_h1_h2_float.json`.

The focused suite passes 21 tests: seven new arity-DP tests, eight
full-arity tests, and six closed-kernel tests. Controls include independent
labelled hole enumeration with noncommuting matrices, scalar mask counts,
all g8 arity caps, the separate full-arity seen-count DP, zero coefficients,
normalization, and resource guards. Independent read-only reviews of both
new helpers found no issue. No native run accompanied these screens.

Reproduce from `theory`:

```text
python -B -m unittest hd_two_sum_gf16_arity_dense_tests hd_two_sum_gf16_full_dense_tests hd_two_sum_closed_kernel_tests
python -B hd_two_sum_gf16_arity_dense.py --output NEW_OUTPUT.json --max-seconds 120
```

The CLI refuses to overwrite an artifact. Frozen SHA-256 pins:

| File | SHA-256 |
|---|---|
| `hd_two_sum_gf16_arity_dense.py` | `539e8c819edc36a351c21e75ca41e91513c0de271d8c23419503ee68c770142a` |
| `hd_two_sum_gf16_arity_dense_tests.py` | `905de1c23697578c2442ccca308b6868965476c4f4854659dd791d9600793421` |
| `hd_two_sum_closed_kernel.py` | `2ee0fe57de750a19466155317d595589695252de450518b436eab6c5c0fbf979` |
| `hd_two_sum_gf16_arity_dense_8_0_all14_float.json` | `c9de91df117fea070552fb4a2b5ab75ef1e7ef265c78b9adf6ec69f77f7c4f99` |
| `hd_two_sum_gf16_full_dense_8_0_h1_h2_float.json` | `7d07826d73e83b465ed985b433fb09933afac5cb01071765f122acfefc80e258` |

Next, test the two-weight coefficient DP for the 105 unordered dense-arity
pairs. Broader source/weight coverage and fresh directed replay remain
separate obligations. No such extension was launched in this iteration.
