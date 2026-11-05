# Two dense groups: a bounded feature-envelope control

The positive feature interface is sound, but its first control is too loose
to support a certificate. For the fixed 36-row-set type described below,
discarding labelled dense-mask placement loses 467.35 bits already with one
dense group. The two-group control also fails. These are FLOAT diagnostics,
not directed bounds or evidence that the construction fails.

No encoder, old witness, helper, or receipt changed. Only two frozen controls
were evaluated; there was no fitting, type-histogram sweep, or native benchmark.

## Fixed source class and comparison measure

The construction is XOR2 over GF(256), with the existing h4 inner map and
independently swapped alternating pair routes. Fix 17 active source groups
with source matching counts `(8,0)`. Of these, `h` groups are dense and
`n=17-h` have two active constituent rows. This experiment uses `h=1` as a
control and `h=2` as its target. Every dense group has its fixed row set in
the same deterministic type `B`: the 36 arity-five row sets assigned to
`a4_mask240_fit` in the frozen twelve-type proposal. Row sets and sparse
constituent pairs remain fixed through all 32 columns.

For a fixed dense row set `A`, let `N_A(T)` count input vectors supported
inside `A` whose image under `I+11^T` has output support `T`. This comparison
allows zero input coordinates. The shortened-RS counting argument gives
the comparison factor `256^(32|A|)/255^(16|A|)`. Consequently define

\[
 W_B(T)^{32}=\sum_{A\in B}\frac{N_A(T)^{32}}{255^{16|A|}}.
\]

Hölder bounds the sum over one fixed `A` across 32 columns by the product
of these column weights. Applying this separately to two fixed dense
groups gives `W_B(T_1)W_B(T_2)` in each column. The ordered choices in `B^2`
are already counted: there is no additional `36^2`, `2!`, or dense comparison
factor outside the column polynomial.

## Positive feature and chronological resource envelope

For a tuple of dense masks, retain three quantities: the number `e` of
nonempty masks, their total weight `t`, and the number `f` of full masks.
The feature mass is the sum of the products of the actual labelled `W_B`
weights. It is a positive convolution, not a symmetric replacement by
binomial mask counts. The two-group control has 26 nonzero features.

Consider a sparse column category with `s` singleton groups and `d` double
groups, of which `k` cancel. Before dense offsets are added, each regional
occupancy is between `d-k` and `s+d`, and its sum over eight regions is
`7s+8d-6k`. Every dense mask tuple in feature `(e,t,f)` therefore satisfies

\[
 d-k+f\le m_i\le s+d+e,\qquad
 \sum_{i=0}^7(m_i-(d-k+f))=7s+2k+t-8f.
\]

The resource DP maximizes componentwise over this enlarged set of occupancy
sequences. It prepends the corresponding regional matrix to each suffix
vector, preserving the eight labelled regions in their original order.
Different vector coordinates may use different maximizing suffixes; positive
matrix multiplication makes this a valid upper envelope.

For each fixed nonempty count `e`, sum the positive feature masses and the
exact cancellation probabilities

\[
 \binom dk\frac{254^{d-k}}{255^d}
\]

before taking the shared-potential quotient. With one global potential
`v=(1,p)` this gives a scalar column cap `c_{s,d,e}`. Form

\[
 P(X,Y)=\sum_{s+d\le n}\sum_e
 \binom ns\binom{n-s}d2^s255^{s+2d}
 c_{s,d,e}X^{s+2d}Y^e.
\]

The implementation uses positive ordered-column convolution to retain the
tail `X>=34n`, `Y>=17h` in `P^32`. The latter is a relaxation of the two
individual dense-column floors when `h=2`, not an assertion that the
individual floors follow from their sum. For this source class, multiply
the tail by

\[
 \frac{K_{17}(8,0)\binom{17}{h}28^n z^{-734003}}
 {255\,p\,255^{32n}},\qquad
 K_{17}(8,0)=15708305994709246428578365440.
\]

Thus the target has floors `(510,34)`, sparse normalization exponent `480`,
and source-color count `136`. The projective and endpoint factors occur
once. The resource cap is universal over the fixed sparse pair assignments;
the `28^n` choices are counted once, not redrawn each column.

## Bounded observations

Both controls use exactly
`z=8973939342378351/9007199254740992` and
`log2(p)=-2386334719/33554432`.

| Control | FLOAT log2 bound | Dominant `(X,Y)` | Elapsed |
|---|---:|---:|---:|
| One dense group, feature envelope | +411.170497929167 | `(1004,32)` | 1.17 s |
| Same one-group type, prior labelled-mask envelope | -56.177151934930 | — | existing artifact |
| Two dense groups, feature envelope | +856.091648792149 | `(938,64)` | 1.80 s |

Each new control used 171 cached occupancy intervals. The two-group operator
build took 0.80 seconds. The sparse floor is not the dominant coefficient
in either failed control. The one-group comparison holds the type, witness,
counts, and positive floors fixed; it exposes loss from enlarging the
labelled dense-placement/sparse-occupancy set, not a new Hölder loss. It
does not establish that any particular retuning could or could not recover
this loss. The third authorized frozen control was not run.

Six tests pass: exact nonuniform labelled feature convolution; exact tiny
Hölder tensorization; noncommuting ordered-resource dominance; exhaustive
tiny fixed sparse geometries with whole dense masks and cancellation coins;
exact integer two-floor column enumeration; and the frozen type/count scope.
The independent source audit found no material mathematical issue.

The next useful control should retain the labelled dense occupancy vector,
whose entries lie in `{0,1,2}`. There are only `3^8=6561` such vectors.
Combining the two dense mask weights by this vector and nonempty count
preserves the required information before a sparse-offset DP. A maximum
width-15 suffix output has `6561*121*2` doubles, about 12.1 MiB; intermediate
storage and runtime still require a bounded implementation check. This is
a proposal for the next turn, not an implemented or measured improvement.

### Proposed labelled-vector algorithm (not implemented)

This section records the proposal made after the coarse control. The later
[labelled-control note](HD_XOR_DENSE_TWO_LABELLED.md) records its isolated
implementation and one bounded two-group result; the coarse artifacts
below remain unchanged.

Compute the positive convolution `M[e,ell]` from the `256^2` mask pairs once,
where `ell_i=1_{i in T_1}+1_{i in T_2}` and `e` counts nonempty masks. The
regional operator depends on this labelled occupancy vector, not on the
names of the two masks that produced it. The aggregate dense floor still
requires retaining `e` separately.

For sparse interval `[l,u]=[d-k,s+d]`, let `w=u-l`. Build one suffix DP over
ternary occupancy codes and shifted sparse sums. Prepending labelled region
`i` and offset `b` uses

\[
 V_{\rm new}[b+3c,j]
 =\max_{0\le\delta\le w}
 R_i[l+\delta+b]V_{\rm old}[c,j-\delta],
\]

with a componentwise maximum. All prefixes share the already computed
suffix indexed by `c`. There are 136 intervals at `n=15`; this is not one
DP per dense mask pair. For a fixed interval the required final targets are
`j=7w-5k`, `0<=k<=w`. At a partial stage, discard a sum only if no target is
reachable using the remaining sparse coordinate range. This pruning is
independent of matrix values and therefore exact for the enlarged resource
set.

After the final stage, multiply by `M[e,ell]` and sum over `ell`; accumulate
the exact binomial cancellation weight into the corresponding `(s,d,e)`
vector. These two positive sums commute. Neither may be replaced by a
potential maximum until both have finished. Then use the same scalar
polynomial and tail as above. This avoids storing a full category-by-offset
bank and never interchanges a sum with the resource maximum.

Counting a two-state matrix action as one DP candidate, the reachable-sum
version requires 150,969,492 candidates across 136 intervals, versus
399,538,080 without pruning. The existing unpruned one-dense-group offset
bank has 24,339,342 candidates. With pruning, the largest single DP array is
398,034 doubles, about 3.04 MiB; temporary arrays and indexing add overhead.
The existing saved one-group FLOAT banks took about 2.1--2.7 seconds each,
including their other work. A 10--30 second estimate for one new FLOAT bank
is only a planning estimate, not a measurement. A first implementation
should enforce a bounded one-bank control before attempting any type sweep.

## Reproduction and pins

From the `theory` directory, use fresh output filenames:

```powershell
C:/Python314/python.exe -m unittest hd_xor_dense_two_feature_tests -v
C:/Python314/python.exe hd_xor_dense_two_feature.py --control h1-olda4 --output fresh_h1.json
C:/Python314/python.exe hd_xor_dense_two_feature.py --control h2-olda4 --output fresh_h2.json
```

SHA-256 pins for this bounded record:

| File in `theory/` | SHA-256 |
|---|---|
| `hd_xor_dense_two_feature.py` | `80faccc10a252513b8ad189a68d376f34c47d7715a88d5496ba3340966dd1455` |
| `hd_xor_dense_two_feature_tests.py` | `c6be0e6f5641f04ca33f5b95367ed9c80c1552f9c4c281d488f2f37b117c41b1` |
| `hd_xor_dense_two_feature_h1_olda4.json` | `0cc10b3f8fee95b05e367edaf9347c652c7d5bb045550818dfee50180fa35bd8` |
| `hd_xor_dense_two_feature_h2_olda4.json` | `84a71526752354ed4794864bc8a228c4778fb0eb14e3829893729c21ebe3d87e` |
| Prior `hd_xor_dense_holder_8_0_h1_arity12.json` | `4577166d2aefbe1912ac6bc096c8f08a3604439c03f72ff0a2dd7ad4fef74a0e` |
