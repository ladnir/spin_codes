# Affine two-sum outer: one dense group

The affine eight-row mixer has a simple sparse support law, but that law does not cover three or more active constituent rows.
This checkpoint uses its exact dense support kernel and keeps the dense row set fixed through all columns.
It concerns the same ideal g8 mode-4 setup as [the sparse-family replay](HD_TWO_SUM_AFFINE_SPARSE_CERTIFICATE.md).

The scope is exactly $17$ active groups in ordered source class $(8,0)$: sixteen groups have two nonzero RS$[32,16]$ constituents, and one has arity $3$ through $8$.
All $219$ dense row sets and all $28^{16}$ fixed sparse pair assignments are included.
For forbidden output weight at most $734003$, fresh directed replays give

$$
\log_2\Pr[\text{a forbidden word in this class}]
\le -\frac{184772727559433}{10^{12}}
=-184.772727559433.
$$

This does not cover other source classes, two dense groups, or the whole code.

## Fixed dense row sets

Write $q=256$, and let $H=I+\mathbf1\mathbf1^T+uu^T$ for $u=(8,\ldots,15)$ in the AES field.
For a fixed input row set $A\subseteq\{0,\ldots,7\}$, let $N_A(T)$ count the vectors in $\mathbb F_q^A$ whose image under $H$ has support exactly $T$.
Input zeros are allowed in this comparison count.
The implementation recomputes $N_A$ from the affine matrix by exact rank calculations and subset inversion.
It checks $N_A(\varnothing)=1$ and $\sum_T N_A(T)=q^{|A|}$.
It does not reuse the old XOR support kernel.

Let $B_a$ be the set of all row masks of arity $a$, for $3\le a\le8$.
The nonzero RS rows admit the comparison factor

$$
\beta_a=\frac{q^{32a}}{255^{16a}}.
$$

For a fixed sequence of output supports $T_1,\ldots,T_{32}$, this factor and the $q^{-a}$ column probabilities give

$$
\sum_{A\in B_a}\beta_a
\prod_{j=1}^{32}\frac{N_A(T_j)}{q^a}
=\sum_{A\in B_a}\frac{\prod_{j=1}^{32}N_A(T_j)}{255^{16a}}
\le\prod_{j=1}^{32}W_a(T_j),
$$

where Hölder's inequality gives the last step and

$$
W_a(T)^{32}:=\sum_{A\in B_a}\frac{N_A(T)^{32}}{255^{16a}}.
$$

The sum over $A$ occurs once, outside the $32$-column product, before applying Hölder.
The bound therefore does not resample the dense row set in each column.
The weights $W_a$ already contain both the dense comparison factors and the fixed-row-set union.
No additional dense factor $\beta_a$ or $|B_a|$ is used.

## Labelled dense masks and sparse holes

Fix the full source set $G$ and its dense group location.
Let $R_i(m)$ be the source-conditioned regional upper matrices from the sparse-family replay.
The source class fixes their parity order: even regions use full-pair count $8$, and odd regions use count $0$.
Let $v=(1,p)^T$.

For a column with $s$ sparse singleton inputs, $d$ sparse double inputs, and dense output mask $T$, set $t=s+d$.
The suffix recurrence keeps every bit of $T$:

$$
V_i^T[k]=\sum_{h=0}^{k}\binom{k}{h}
R_i\bigl(t-h+\mathbf1_{i\in T}\bigr)V_{i+1}^T[k-h],
$$

with $V_8^T[0]=v$ and $V_8^T[k]=0$ for $k>0$.
The sparse hole mixture is

$$
U_{s,d,T}=\frac1{255^d}
\sum_{k=0}^{d}\binom{d}{k}247^{d-k}V_0^T[k].
$$

The binary suffix implementation shares computation between masks with the same suffix.
It preserves their labelled order and does not commute regional matrices.
For each fixed mask, the recurrence is independent of $d$ until the final positive mixture.

For arity class $a$, define the empty and nonempty vector bounds

$$
E_{a,s,d}=W_a(\varnothing)U_{s,d,\varnothing},
\qquad
F_{a,s,d}=\sum_{T\ne\varnothing}W_a(T)U_{s,d,T}.
$$

Only after both positive mixtures does the code take the common-potential quotients

$$
c^0_{a,s,d}=\max_{j\in\{0,1\}}\frac{(E_{a,s,d})_j}{v_j},
\qquad
c^1_{a,s,d}=\max_{j\in\{0,1\}}\frac{(F_{a,s,d})_j}{v_j}.
$$

One global $z$ and $p$ apply to all masks, columns, and six arity classes.
The regional upper matrices dominate every allowed subset of the fixed source set.
Consequently, averaging the joint mask and hole assignments does not assume independent regional subsets.

## Support floors and probability bound

For $e\in\{0,1\}$, form

$$
P_a^e(X)=\sum_{s+d\le16}
\binom{16}{s}\binom{16-s}{d}2^s255^{s+2d}
c^e_{a,s,d}X^{s+2d}.
$$

The thirty-two sparse constituent words require total support at least $544$.
The dense group is nonzero in at least $17$ columns because it contains a nonzero RS word and $D_2HD_1$ is invertible.
The code retains both floors in the positive expansion of

$$
\bigl(P_a^0(X)+Y P_a^1(X)\bigr)^{32}.
$$

The resulting arity bound is

$$
\frac{K_{17}(8,0)\,17\,28^{16}}{255\,p\,255^{512}}z^{-734003}
\sum_{J\ge544}\sum_{D\ge17}
[X^JY^D]\bigl(P_a^0(X)+Y P_a^1(X)\bigr)^{32},
$$

where $K_{17}(8,0)=15708305994709246428578365440$.
The source count, dense-group location, sparse pair choices, and projective normalization each appear once.
Summing the six disjoint arity classes completes the stated scope.

## Directed replay

The frozen witness is

$$
z=\frac{8970988817957597}{9007199254740992},\qquad
\log_2p=-\frac{2266954815}{33554432}.
$$

The proposal is `theory/hd_two_sum_affine_dense_8_0_two_frozen.json`, SHA-256 `922d920286a6a513af0ae1a98ba41f51a5379c6fa26de2b9c249d5d64a393b9a`.
The replay reads the exact winning `additional` object and ignores all saved numerical estimates.
For each exact rational moment $N/D$, Arb proposes an upper dyadic endpoint $w$ for its $32$nd root.
An exact integer check verifies $w^{32}D\ge N$ before using that endpoint.

All subsequent operations use Arb, including the positive hole sums, dense mask sums, quotients, polynomial powers, and logarithms.
Maxima compare explicit upper endpoints.
The independent verifier sums the six rounded arity bounds at 512-bit precision.

Fresh 256- and 384-bit replays agree on every arity ceiling and the raw aggregate
$-92386363779717/500000000000$.
Their proof-computation times were $3.9162$ and $3.6199$ seconds.
The rounded-class union gives the slightly larger ceiling stated above.
Six focused tests check exact labelled hole enumeration, noncommuting order, affine moments, exact root endpoints, once-only normalization, and fail-closed scope.

| Artifact in `theory/` | SHA-256 |
| --- | --- |
| `hd_two_sum_affine_dense_certificate.py` | `eae905d1e713975e4e5c1fa6c638c670174fa8a84bf24cd4a861d7d89612fd68` |
| `hd_two_sum_affine_dense_h1_receipt256.json` | `41a4c2d606e230e9c7d976a67d4ff6420d1c2a3b67acd94deddb8964009234be` |
| `hd_two_sum_affine_dense_h1_receipt384.json` | `9acbc2fc144ace27a377a351ee723bba1f54cbdaa2e8160a92d443b8e2507950` |
| `hd_two_sum_affine_dense_h1_dual_verified.json` | `f202d10666a8e6d3b166c928d68a86aa7b5225fb0dbe9d22698be4ad6035e3a8` |

The next proof step is coverage of the remaining source classes and dense-group cohorts.
No encoder benchmark, native change, or promotion is part of this checkpoint.
