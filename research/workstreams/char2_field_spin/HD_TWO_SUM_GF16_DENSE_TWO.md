# G16 two-dense-group screen by two weight coefficients

This extends the fixed-source h1 method in
[HD_TWO_SUM_GF16_ARITY_DENSE.md](HD_TWO_SUM_GF16_ARITY_DENSE.md) to two dense
groups. The construction remains the frozen Group16 candidate: sixteen
RS[16,8] rows, the `06 + GF16` two-sum mixer, fair-swapped paired routing,
and the existing h4 inner code. No native code or frozen proof helper is
changed.

The scope is 17 active source groups, ordered matching counts `(8,0)`, and
the cohort of fifteen two-row groups plus two dense groups. Each dense
group has a fixed row set of size 3–16. The screen partitions their
65,399² ordered row-set assignments into 105 unordered arity pairs.
It uses one inherited witness, with no fitting:

\[
z=\frac{8970988817957597}{9007199254740992},\qquad
\log_2p=-\frac{2266954815}{33554432}.
\]

The implementation is FLOAT-only. A passing value is proposal evidence,
not a directed receipt or a whole-code certificate.

All 105 separate arity classes pass at this literal. Their union is
**−476.62769692618843 bits**. The merged-arity control fails at
**+258.67213841832563 bits**; both results are retained.

## Two weights preserve labelled regional order

For a dense group of arity \(a\), the h1 note defines the per-mask Hölder
weight

\[
W_a(w)=\left(\frac{\sum_{|A|=a}N_A(T)^{16}}{255^{8a}}\right)^{1/16},
\qquad |T|=w.
\]

The closed kernel makes this depend only on \(w\). For two fixed dense
groups, apply the two bounds separately. Their per-column comparison
weight is \(W_a(w_1)W_b(w_2)\); their row sets remain fixed across columns.

Fix a sparse category \((s,d)\), and put \(t=s+d\). Let \(R_i(m)\) be the
source-conditioned regional matrix envelope from the h1 interface. The
suffix coefficient \(V_i(w_1,w_2,k)\) records two dense weights and a count
of labelled sparse holes. Its terminal value is \(v=(1,p)\) at
\((w_1,w_2,k)=(0,0,0)\), and zero elsewhere. The recurrence is

\[
V_i(w_1,w_2,k)=
\sum_{b_1,b_2\in\{0,1\}}\sum_{j=0}^{k}
\binom{k}{j}R_i(t-j+b_1+b_2)
V_{i+1}(w_1-b_1,w_2-b_2,k-j).
\]

The `01` and `10` branches use the same matrix but increment different
indices. Both must be added. The DP applies the matrices in their original
order, even when they do not commute. Its coefficients sum every labelled
mask pair once, so neither weight index receives an extra binomial factor.

Mix the hole count with \(\binom{d}{k}239^{d-k}/255^d\), producing
\(D_{s,d}(w_1,w_2)\). For \(e\in\{0,1,2\}\), form

\[
c_e^{a,b}(s,d)=\max_{j\in\{0,1\}}\frac{1}{v_j}
\sum_{\mathbf1_{w_1>0}+\mathbf1_{w_2>0}=e}
W_a(w_1)W_b(w_2)D_{s,d}(w_1,w_2)_j.
\]

Thus \(e\) counts nonempty dense group-columns. Both the sparse-hole and
dense-mask mixtures precede the potential quotient. There are 17² weight
coefficients rather than a table of \(4^{16}\) labelled mask pairs.

## Column powers and once-only counts

For each \(e\), define

\[
P_e^{a,b}(X)=\sum_{s+d\leq15}
\binom{15}{s}\binom{15-s}{d}2^s255^{s+2d}
c_e^{a,b}(s,d)X^{s+2d}.
\]

Raise \(P_0+YP_1+Y^2P_2\) to power 16 and retain \(X\geq270\), \(Y\geq18\).
The first floor sums the distance-9 bounds for 30 nonzero sparse rows.
Each dense group occupies at least 9 columns. Their individual constraints
imply the retained aggregate floor 18; the converse is not assumed.

Write \(K=15708305994709246428578365440\) for the exact source-support
count at `(8,0)`. Multiply the positive tail by

\[
\frac{K\binom{17}{2}120^{15}}
 {255\,p\,255^{240}\,z^{734003}}\,m_{a,b},
\qquad
m_{a,b}=\begin{cases}1&a=b,\\2&a<b.\end{cases}
\]

The two dense source positions are ordered canonically. When \(a<b\),
both assignments of the arities to those positions occur. The universal
regional envelopes and the two-weight sum give the same bound for them,
so their factor 2 is paid once after the column power. It is not an
assumption that the actual fixed source subsets are exchangeable.

The factors \(\binom{16}{a}\) and \(\binom{16}{b}\) are already inside
the Hölder moments. The identity

\[
\sum_{3\leq a\leq b\leq16}
m_{a,b}\binom{16}{a}\binom{16}{b}=65399^2
\]

checks the complete ordered row-set partition.

## Merged-arity control

The same vector bank also supports one envelope for all dense row sets:

\[
W_H(w)^{16}=
\sum_{a=3}^{16}\frac{\sum_{|A|=a}N_A(T)^{16}}{255^{8a}},
\qquad |T|=w.
\]

The helper constructs this sum as exact integer numerators over the common
denominator \(255^{128}\), then evaluates its roots in FLOAT arithmetic.
Using \(W_H(w_1)W_H(w_2)\) covers all ordered row-set assignments in one
contraction. It receives no arity multiplicity or other row-set factor.
It may be looser than the sum of the 105 separate bounds because Hölder
now also groups different arities.

## Results and validation

[hd_two_sum_gf16_dense_two.py](theory/hd_two_sum_gf16_dense_two.py) and its
[tests](theory/hd_two_sum_gf16_dense_two_tests.py) are isolated new files.
Nine tests pass in 0.890 seconds. They include exact labelled mask/hole enumeration with
noncommuting matrices, scalar coefficient counts, shared-potential order,
the 105-class partition, aggregate moments against all g8 row sets, zero
operators, and scope/resource guards. The full-arity specialization agrees
with the separate seen-count h2 DP.

The fixed-literal run completed in 16.188 seconds. It built one regional
bank and one two-weight bank; the latter stage completed at 1.261 seconds.
The independent full-arity seen-count control agrees with every full-arity
column cap within \(5.12\cdot10^{-13}\) bits.

| Separate arity class or control | Failure log2 screen |
|---|---:|
| `(16,16)` | −476.8021489572043 |
| `(15,16)`, including its factor 2 | −479.8473214623341 |
| `(15,15)` | −484.89249396746345 |
| `(14,16)`, including its factor 2 | −484.98560337185427 |
| All 105 separate classes | **−476.62769692618843** |
| One merged-arity envelope | **+258.67213841832563** |

Every class has its largest retained term at \((X,Y)=(270,18)\). The
full-arity result also reproduces the earlier independent screen
`−476.80214895720337` to FLOAT precision. The complete raw control is
`hd_two_sum_gf16_dense_two_8_0_all105_float.json`.

The merged envelope loses about 735.30 bits relative to the separate
union. Its empty-mask weights explain the scale of this loss:

\[
W_a(0)=\binom{16}{a}^{1/16}255^{-a/2}.
\]

The merged empty mass is dominated by arity 3, whereas large output masks
prefer arity 16. At the peak \(Y=18\), there are 14 empty dense
group-column slots among 32. Replacing the arity-16 empty weight by the
arity-3 weight changes their contribution by approximately

\[
14\left(\frac{13}{2}\log_2 255+
          \frac{\log_2 560}{16}\right)=735.47\text{ bits}.
\]

This calculation diagnoses the comparison loss; it is not an equality
for the whole tail. Grouping all arities lets different columns favor
different arities in the upper envelope. It does not establish a failure
of the construction. Preserve the separate-arity result and the failed
merged control instead of treating the latter as a replacement.

The CLI retains a partial result on its explicit 120-second budget and
refuses to overwrite an artifact:

```text
python -B -m unittest hd_two_sum_gf16_dense_two_tests
python -B hd_two_sum_gf16_dense_two.py --output NEW_OUTPUT.json --max-seconds 120
```

Frozen SHA-256 pins:

| File | SHA-256 |
|---|---|
| `hd_two_sum_gf16_dense_two.py` | `d8a4976db2c7b69def936ec7855791b6efc52a065318835b37bee79e9cab04d4` |
| `hd_two_sum_gf16_dense_two_tests.py` | `0b31c9df07524bd3db7c061b6ebea20d7e8565fa62be378b83436fb688b54080` |
| `hd_two_sum_gf16_dense_two_8_0_all105_float.json` | `fea001af482281e7f6f7372e621be230501c858a2b299c2aa3a5a7426e22486e` |

After this gate, the useful next step is broad source/cohort coverage and
directed replay of a passing comparison, retaining arity information or
testing small arity classes. No optimization or native run accompanied
this screen.
