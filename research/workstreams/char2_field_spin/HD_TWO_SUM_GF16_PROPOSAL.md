# A sixteen-row affine two-sum proposal

A suitable sixteen-label affine coset exists in the AES field. It gives an
involutory two-sum mixer with full singleton output and the same one-hole
law for every input pair. The proposal now has a default-off implementation
and a bounded sparse proof screen. Neither supplies a whole-code distance
certificate; the existing checkpoints remain unchanged.

The proposed outer group has sixteen RS[16,8] constituents instead of eight
RS[32,16] constituents. Both have 128 input and 256 output field coordinates.
The shorter RS transforms offer a circuit saving, but their changed support
floors require a new proof. The existing checkpoints remain unchanged.

## One explicit AES-byte coset

Work in `F=GF(256)` with polynomial `0x11b`. Its subfield `K=GF(16)` is

```text
00 01 0c 0d 50 51 5c 5d b0 b1 bc bd e0 e1 ec ed
```

The ordered basis `(01,0c,50,b0)` gives this list by binary subset sums.
Choose `c=06`. The sixteen labels `u_i` are the corresponding elements of
the affine coset `C=c+K`:

```text
06 07 0a 0b 56 57 5a 5b b6 b7 ba bb e6 e7 ea eb
```

Define `H=I+11^T+uu^T`. Thus
`H_ij=delta_ij+1+u_i*u_j`, with addition and multiplication in `F`.

Here is a criterion that explains the choice, rather than merely checking
the list. Write `Tr` for the absolute trace from `K` to GF(2), and put
`t=c+c^16`, `n=c^17`. Every element of `C` has relative trace `t`. For this
choice,

```text
t = bc,  n = b1,  Tr(1/t) = 0,  Tr(n/t^2) = 1.
```

For any `a` in `C`, the relative trace of `a^-1` is `t/N(a)`. Hence its
inverse belongs to `C` exactly when `N(a)=1`. Writing `a=c+v`, this equation
becomes `v^2+t*v+n=1`. After division by `t^2`, its solvability in `K` is
equivalent to `Tr((n+1)/t^2)=0`. But

\[
 \operatorname{Tr}((n+1)/t^2)
 =\operatorname{Tr}(n/t^2)+\operatorname{Tr}(1/t)
 =1.
\]

Thus `C` and `C^-1` are disjoint. The identity
`Tr(1/t^2)=Tr(1/t)` uses invariance of the absolute trace under squaring.
More generally, every nonzero relative trace `t` satisfying `Tr(1/t)=0`
gives such a coset. The concrete check found 112 eligible choices of `c`,
representing seven cosets; only the one above is proposed here.

## Involution and sparse support law

The labels satisfy `sum(u_i)=sum(u_i^2)=0`, and 16 is zero in characteristic
two. Therefore, for the matrix `L=[1,u]`, we have `L^T L=0`, and
`H^2=(I+LL^T)^2=I`.

Every entry is nonzero. A diagonal entry is `u_i^2`; an off-diagonal entry
could vanish only if two coset elements were inverses. Fix two distinct
input columns `i,j`. Their minors have three forms:

| Selected output rows | Determinant |
|---|---|
| `i,j` | `1` |
| `k,l` both outside `i,j` | `(u_i+u_j)(u_k+u_l)` |
| `i,k`, with `k` outside `i,j` | `1+u_i(u_i+u_j+u_k)` |

The second form is nonzero because the labels are distinct. In the third
form, the triple sum remains in `C`; inverse-disjointness therefore excludes
zero. The exchanged mixed case is identical. Thus every two-row minor of
every two-column restriction is nonzero.

A singleton input consequently has all sixteen outputs nonzero. For a fixed
pair, each output row has one nonzero input ratio that makes it zero. The
minor condition makes these sixteen ratios distinct. Independent nonzero
input diagonals make the ratio uniform over 255 values. The exact pair law is

\[
 \Pr[\text{full support}]=239/255,\qquad
 \Pr[\text{only output }i\text{ is zero}]=1/255\quad(0\le i<16).
\]

Independent nonzero output diagonals restore independent nonzero amplitudes
conditional on the complete support mask. Fixed input row sets remain fixed
across all sixteen columns. No row-set resampling is needed for this law.

An exact local check found zero zero-entries among 256 entries, zero zero-minors
among 14,400 minors, and `H^2=I`. It independently enumerated all 30,600
projective pair inputs and obtained the stated histogram for all 120 pairs.
The inverse-coset and circuit-count arguments also received an independent
read-only review.

## Circuit opportunity, not a timing prediction

Let `b=(01,0c,50,b0)`, and index the inputs by four binary bits. Compute
`S=sum_i x_i` and the four bit sums `A_j=sum_{i:bit_j(i)=1} x_i`. Then

\[
 T=cS+A_0+b_1A_1+b_2A_2+b_3A_3,
 \qquad (Hx)_i=x_i+S+cT+\sum_{j:bit_j(i)=1}b_jT.
\]

Computing `T` uses four fixed field products. Computing `cT,b_1T,b_2T,b_3T`
uses four more. Thus the proposed mixer needs eight fixed products per
sixteen-output column. The common-scalar diagonal gauge still gives 31
variable products per column.

The five input sums admit 26 binary XORs by recursive pair sharing: eight
pair sums, seven XORs for the odd-position sum, and eleven XORs for the
remaining sums on the eight pair sums. Combining `T` takes four more;
forming `S+cT` takes one; the sixteen corrections take fifteen; applying
them takes sixteen. This is a 62-binary-XOR circuit, not a compiled instruction
count. It does not include register moves, loads, spills, or setup work.

Using the existing zero-coefficient FFT simplifications, an RS transform with
dimension `k` and length `2k` has `k*log2(k)-k+1` fixed products. The resulting
product counts for equal 128-input layouts are:

| Work | Eight RS[32,16] constituents | Sixteen RS[16,8] constituents |
|---|---:|---:|
| RS transforms | `8*49 = 392` | `16*17 = 272` |
| Normalized input/output diagonals | `32*15 = 480` | `16*31 = 496` |
| Affine two-sum mixer | `32*6 = 192` | `16*8 = 128` |
| Total field products | `1064` | `896` |

The difference is 168 products per such circuit batch, or 344,064 products
over 2,048 batches under the current K20 batching convention. This does not
predict a wall-clock saving. The sampled diagonal key count increases by
16 per group, about 3.33%; payload size is unchanged. Larger live columns
and a different RS schedule may alter register pressure and memory traffic.

For comparison, the same binary-XOR counting convention gives 113 XORs per
RS[32,16] constituent and 41 per RS[16,8] constituent. The mixer alone adds
32 XORs across the equal-size layout (`16*62` versus `32*30`), while the RS
circuits save 248. The net algebraic count saves 216 binary XORs; compiler
instruction combining can change this comparison.

## The proof risk from shorter constituents

An active RS[16,8] row has support at least nine, not seventeen. For a sparse
cohort with `A` active constituent rows, the counting-measure normalization
becomes `255^(-8A)` and the aggregate support floor becomes `9A`. There are
sixteen column factors rather than thirty-two. The hole recurrence can use
sixteen labelled regional factors and the new `239/255` full-support mass,
but its witnesses and final bounds must be recomputed.

Dense cohorts are essential. In fact, for every setup, four constituent
rows already give a deterministic outer-weight upper bound of 76. Their
message space has dimension 32. Requiring both column sums `S,T` to vanish
on fifteen selected columns imposes at most thirty homogeneous linear
constraints. A nonzero message remains. On those fifteen columns the mixer
fixes its input, so at most four output coordinates are nonzero; the last
column has at most sixteen. Hence the outer output has weight at most
`15*4+16=76`, below 35% of 256. Nonzero diagonal scalings do not change this
argument. The inner map would have to supply the required expansion.

The appropriate next proof gate is therefore a new RS[16,8] control, not a
transfer of an old receipt or a claim based only on the improved sparse law.

## First implementation and proof controls

The [native experiment](HD_TWO_SUM_GF16_NATIVE.md) passes its independent
scalar oracle and full-size correctness checks. Its initial producer is
slower: 4.044 ms versus 3.590 ms for the eight-row affine candidate in the
same serial three-round comparison. The unchanged Toeplitz control measures
3.772 ms. The algebraic operation saving has not produced a timing win.
The first implementation uses the existing Group16 producer, which lacks
the fused final RS stage used by Group8. Assembly inspection instead finds
an outlined store helper that materializes every mixed sixteen-vector column.
Removing that interface preserves all mathematics and removes 64 MiB of
intermediate stack traffic. In a fresh matched comparison, the revised
candidate takes 3.398 ms versus 3.616 ms for Group8 and 4.067 ms for the
original Group16 implementation. The revision does not change the RS schedule.

The first proof control fixes 17 active groups and covers all 81 source
matching classes and all 18 mixtures of one- and two-row groups. It does
not cover groups with three or more active rows. For each fixed source
set and fixed constituent row sets, it averages the exact one-hole law
over independent column diagonals. The regional matrices universally bound
each subset of the fixed source set; their chronological order is retained.

For a column containing `s` one-input groups and `d` two-input groups, set
`b=s+d`. The ordered positive recurrence is the existing hole recurrence
with sixteen regional factors and no-hole mass `239/255`. Its completed
action is bounded by one common potential. The column support polynomial
is raised to the sixteenth power, with support floor `9*(n1+2*n2)` and
normalization `255^(8*(n1+2*n2))`. The once-only multiplicity is
`K_17(a,b)*C(17,n1)*16^n1*120^n2`. The projective factor is applied once.

One inherited witness, without fitting, gives a FLOAT log2 union bound
of -511.43278211723936 over all 1,458 cells in 18.92 seconds. The worst
source class is `(0,0)`. This is a proposal for a scoped family, not directed
arithmetic or a full certificate. The improved sparse behavior does not
remove the dense four-row risk described above.

Six tests pass: the coset and all projective pair laws, all 128 binary
basis inputs, random circuit checks, the silent triples, ordered-region
envelopes, and exact noncommuting hole mixtures. The screen also checks
the new RS dimensions and once-only multiplicities.

Reproduce from `theory`:

```text
python -m unittest hd_two_sum_gf16_sparse_tests
python hd_two_sum_gf16_sparse.py --all-sources --output NEW_RESULT.json
```

The frozen helper SHA256 is
`e225f2fccb2f39d95dbfd7c6896c48957c3fd29898f4f574f406a090f6b95c17`.
The saved artifact `hd_two_sum_gf16_sparse_r17_all81_inherited_float.json`
has SHA256
`d609d3521631cad631d38678be40261c750b00948bc31effe77eebeb97d06380`.

## Reproduce the exact local checks

From `theory`, the existing AES arithmetic helper suffices:

```python
from collections import Counter
from itertools import combinations
from gf8_packet_maps import MUL, INV, power, xor_sum, matmul

K = [x for x in range(256) if power(x, 16) == x]
u = [0x06 ^ v for v in K]
H = tuple(tuple(int(i == j) ^ 1 ^ MUL[a][b]
                for j, b in enumerate(u)) for i, a in enumerate(u))
pairs = list(combinations(range(16), 2))
assert all(H[i][j] for i in range(16) for j in range(16))
assert not set(u).intersection(INV[a] for a in u)
assert matmul(H, H) == tuple(tuple(int(i == j) for j in range(16)) for i in range(16))
assert all(MUL[H[k][i]][H[l][j]] ^ MUL[H[k][j]][H[l][i]]
           for k, l in pairs for i, j in pairs)
expected = {65535: 239, **{65535 ^ (1 << i): 1 for i in range(16)}}
for i, j in pairs:
    counts = Counter(sum((H[k][i] ^ MUL[H[k][j]][ratio] != 0) << k
                         for k in range(16)) for ratio in range(1, 256))
    assert counts == expected
```
