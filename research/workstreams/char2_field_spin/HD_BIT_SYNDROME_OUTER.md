# Four additional XORs for five outer syndromes

This proposal strengthens the experimental G16 outer without adding field
products, keys, or a data pass. Exact local checks pass. A subsequent
[native experiment](HD_BIT_SYNDROME_NATIVE.md) measures about 3.4 ms,
and a [positive column majorant](HD_BIT_SYNDROME_KERNEL.md) covers both
quiet and nonquiet words. No whole-code certificate is established.
The certified Toeplitz encoder and every earlier receipt remain unchanged.

The motivating failure is the selected arity-eight class in
[the closure record](HD_34MS_CLOSURE.md). Its latest coherent-row-set bound
remains positive. The candidate below changes the quiet subspace itself;
it does not fit another envelope to the same two-sum construction.

## Map and setup distribution

Keep sixteen RS[16,8] constituents and the existing ordered AES-field
labels

\[
 u_i=c+\sum_{j=0}^3 b_j a_j(i),\qquad
 c=\texttt{06},\quad b=(\texttt{01},\texttt{0c},\texttt{50},\texttt{b0}),
\]

where `a_j(i)` is bit `j` of `i`, embedded as zero or one in GF(256).
These are the existing labels `06+GF16`; they are neither zero nor one,
and the coset is disjoint from its inverse.

Replace the two-sum matrix by

\[
 H_*:=I+\mathbf1\mathbf1^T+uu^T+\sum_{j=0}^3 a_j a_j^T.
\]

The ideal per-column mixer is `D2 H* D1`. Both diagonals have independent
uniform nonzero field entries, sampled once at setup, independently
between groups and columns.
The existing common-scalar normalization of D1 can be retained: multiplying
all input scales by one scalar and dividing all output scales by that
scalar preserves this linear map. The route and inner map are not changed
by this proposal. In particular, no post-refresh-coupled timing is assumed.

Put `V=[1,a0,a1,a2,a3]` and `w=(c,b0,b1,b2,b3)^T`. Then

\[
 H_*-I=V B V^T,\qquad B=I_5+ww^T.
\]

Every entry of `V^T V` is zero: the relevant binary counts are 16, 8, or 4.
Consequently `(H*-I)^2=0` and `H*^2=I`. Moreover,

\[
 \det B=1+w^Tw=(c+b_1+b_2+b_3)^2=\texttt{19}\ne0.
\]

The matrix V has rank five. Thus a column vector x is quiet, meaning
`H*x=x`, exactly when all five affine-coordinate sums `V^T x` vanish.
This equivalence concerns the pre-output-diagonal vector, not equality
after applying D2.

Every matrix entry remains nonzero. Off the diagonal, a binary dot product
of one leaves `u_i u_j`; a dot product of zero leaves `1+u_i u_j`, which
is nonzero by inverse-disjointness. A diagonal entry is either `u_i^2` or
`u_i^2+1`. Neither vanishes. Every singleton therefore still has full
sixteen-coordinate output support.

## Exact fixed-support quiet counts

Fix a set S of active input coordinates in one column. Its affine rank is
the rank of the augmented binary labels `(1,i)` for `i in S`. Binary-matrix
rank is unchanged on extension from GF(2) to GF(256). The number of vectors
supported within S whose five syndromes vanish is

\[
 256^{|S|-\operatorname{rank} V_S}.
\]

Möbius inversion on subsets removes zero coordinates. The exact number
with support equal to S is

\[
 Z_*(S)=\sum_{T\subseteq S}(-1)^{|S|-|T|}
             256^{|T|-\operatorname{rank} V_T}.
\]

The helper computes all 65,536 integers once, using a positive kernel-size
table followed by exact integer inversion. These are counts, not
probabilities. Any use in a distance bound must retain the appropriate
input-diagonal law, RS counting normalization, and fixed row set across
all sixteen columns.

Any three distinct augmented labels are independent. Quiet nonzero
columns therefore have weight at least four. Four positions forming an
affine parallelogram admit 255 quiet full-support vectors, so this minimum
is attained. This replaces the old quiet triples; it does not eliminate
all quiet columns.

At arity eight, exactly 30 sets have affine rank four: the affine
hyperplanes. The other 12,840 sets have rank five. Exact full-support
quiet counts refine this rank classification:

| Fixed eight-row sets | Quiet full-support count per set |
|---:|---:|
| 1,920 | 0 |
| 10,080 | 16,386,810 |
| 840 | 16,516,605 |
| 30 | 4,162,573,845 |

The old two-sum count is 272,798,207,821,845 for every eight-row set. The
reduction is about 16 bits even for the largest new count, and about
24 bits for the two nonzero rank-five classes. These are local counting
gains, not distance-margin bits. Whole-word constraints and live-state
paths still require a proof.

For every eight-row set A, output supported within A is equivalent to
quietness. If A is not an affine hyperplane, both `V_A` and `V_(A^c)` have
rank five, so `H*_(A^c,A)` also has rank five. In each of the thirty
hyperplane cases, an exact independent field-rank check gives rank four,
equal to `rank V_A`. The kernel of the cross restriction therefore equals
the quiet kernel in every case. This excludes an additional nonquiet
branch confined to A; it does not give the old nonquiet output-weight
floor of nine.

## Pair regression and receipt boundary

The universal one-hole law does not survive. The table lists output-hole
histograms over the 255 nonzero input ratios for each fixed pair. The last
column counts the source pairs with that histogram.

| Ratios by number of holes | Source pairs |
|---|---:|
| 0:239, 1:16 | 12 |
| 0:240, 1:14, 2:1 | 13 |
| 0:241, 1:12, 2:2 | 2 |
| 0:242, 1:12, 4:1 | 68 |
| 0:243, 1:10, 2:1, 4:1 | 12 |
| 0:244, 1:10, 6:1 | 8 |
| 0:245, 1:6, 2:3, 4:1 | 1 |
| 0:245, 1:8, 2:1, 6:1 | 4 |

All 30,600 projective pair inputs were checked exactly. Six holes really
occur, for example for source pair `(3,5)`. The helper returns complete
labelled mask counts, not just these weight summaries. A fixed pair must
remain fixed across columns in a new bound.

Independent D2 entries still make output amplitudes independent and
nonzero conditional on the complete support mask. The invertibility-based
RS counting comparison also remains applicable. Neither fact transfers
an old distance receipt: the old one-hole recurrence, rank-two closed
kernel, and quiet/nonquiet decomposition describe a different outer map.
Inner caps and source-route counts may be reused only through their
unchanged input interfaces. No bound for the certified Toeplitz ensemble
is relabelled.

## Explicit circuit delta and performance risk

The native two-sum circuit already computes

\[
 S=\sum_i x_i,\quad A_j=\sum_{i:a_j(i)=1}x_i,\quad
 T=cS+\sum_j b_j A_j.
\]

Its correction at output i is `S+cT+sum_j a_j(i)b_jT`. Replace its four
shared output deltas by

\[
 T+A_0,\quad b_1T+A_1,\quad b_2T+A_2,\quad b_3T+A_3.
\]

Keep the common term `S+cT` and the same subset-XOR output expansion.
This is exactly H*. An instrumented scalar circuit counts eight fixed
field products and 66 binary XORs, versus eight and 62 for the current
two-sum circuit. Thus the algebraic total remains 896 field products per
128-input batch. The additional 64 binary XORs per batch give 131,072 at
the current K20 batching count. No diagonal key, store, or payload changes.

These are circuit counts, not timing predictions. In particular, the four
input bit sums must remain live until the output deltas are formed. The
old compiler could release them after computing T. A future implementation
must inspect register allocation and preserve the existing inlined store
structure. This candidate is more performance-plausible than adding a
third Vandermonde sum, which would add eight fixed field products per
column, but it has not been benchmarked.

## Checks and next gate

Eight exact tests pass in 0.596 seconds. They verify the factorization,
Gram matrix, determinant, involution, singleton entries, all 128 binary
basis inputs, additional deterministic random inputs, every projective
pair, the affine ranks, and the thirty exceptional cross restrictions.
Independent GF(2)/GF(4) toy enumeration checks
the forced-nonzero count inversion. The complete quiet count table is
nonnegative and sums to `256^11`, the size of the full quiet subspace.

Reproduce from `theory`:

```text
python -m unittest hd_bit_syndrome_outer_tests
python hd_bit_syndrome_outer.py
```

Source SHA256:
`a2c2d0fee734fad2c4b74da7850589f8e84180798dcf5f79beab30e06939292d`.
Test SHA256:
`d66f36f2172f5a1f1bcd6a02f0bdba514700eb102c93ab24287ddb2c47fbfa05`.

The next decisive proof gate should use these exact quiet counts and
labelled pair masks before attempting a full provider. It must also
construct a new positive bound for nonquiet columns. A successful quiet
gate alone is insufficient. No numerical distance screen, native edit,
benchmark, or certificate run was performed for this proposal.
