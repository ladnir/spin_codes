# Affine two-sum outer: directed sparse-family replay

This checkpoint concerns the new eight-row outer mixer with labels $u=(8,\ldots,15)$ in AES $\mathbb F_{256}$.
It does not transfer an old XOR or Toeplitz certificate.
The target family has exactly $17$ active outer groups, each containing either one or two nonzero constituent words.
The implementation replays all $81$ ordered source classes and all $18$ sparse cohort splits.
The family does not include a group with three or more nonzero constituent words, or any other active-group count.

The frozen proposal is `theory/hd_two_sum_affine_sparse_r17_all81_inherited_float.json`.
It uses only

$$
z=996248217/10^9,\qquad p=2^{-127/2}.
$$

Its FLOAT union is $-246.92130644901977$ bits.
The directed implementation reads its exact scope and witness, but ignores every saved FLOAT estimate.
The implementation and tests are `theory/hd_two_sum_affine_sparse_certificate.py` and `theory/hd_two_sum_affine_sparse_certificate_tests.py`.
Fresh 256- and 384-bit replays agree on every one of the $1458$ cell ceilings and all $18$ cohort ceilings.
The independent 512-bit sum of rounded cell bounds gives

$$
\log_2\Pr[\text{a forbidden word in this family}]
\le -\frac{61730326612253}{250000000000}
=-246.921306449012.
$$

The forbidden weight is at most $734003$ among $2^{21}$ output records.
The probability is over the ideal independent setup described below.

## Construction and conditional support law

Each group contains eight RS$[32,16]$ constituents over $\mathbb F_{256}$.
At each of the $32$ columns, the outer applies $D_2HD_1$, where

$$
H=I+\mathbf1\mathbf1^T+uu^T.
$$

Both diagonals have independent uniform nonzero entries, independently across groups and columns.
The existing scalar gauge normalization preserves each sampled matrix.
The route uses the alternating adjacent-pair matchings, independent uniform pair permutations, and independent fair swaps.
The inner is the existing $h=4$ DHD-refresh construction, not the experimental XOR refresh.

Fix the group locations, their constituent row sets, and all constituent support arrays before sampling setup.
The row sets remain fixed through all $32$ columns.
A column with one active input in a group has all eight outputs nonzero.
A column with two active inputs has full output support with probability $247/255$.
Each of the eight supports missing exactly one labelled output has probability $1/255$.
This law is identical for every fixed input pair.
Independent $D_2$ entries make the output amplitudes independent uniform nonzero elements, conditional on the full support assignment.

The exact local algebra is checked in `hd_two_sum_affine_outer.py` and its tests.
Those tests cover every fixed input pair and every nonzero projective ratio.
The uniform pair law removes the old constituent-pair multiplicity strata; it does not remove the dependence between regional subsets.

## Ordered column bound

Let $G$ be a fixed set of $17$ active group locations.
Write $(a,b)$ for its full-pair counts in the even and odd matchings.
For region $i$ and occupancy $m$, let $R_i(m)$ be the nonnegative two-state upper transfer matrix from `hd_xor_pair_route.py`.
The matrix averages the regional pair permutation, swaps, and nonzero amplitudes.
For source full-pair count $F\in\{a,b\}$, it takes an entrywise upper envelope over

$$
\max(0,m-17+F)\le f\le\min(F,\lfloor m/2\rfloor).
$$

Every size-$m$ subset of $G$ satisfies this interval.
Even regions use $F=a$ and odd regions use $F=b$.
The eight matrices remain in chronological order.

Suppose a column has $s$ groups with one active input and $d$ groups with two active inputs.
Set $t=s+d$ and $v=(1,p)^T$.
The exact positive mixture of regional upper matrices satisfies

$$
Q_{s,d}v
=\frac{d!}{255^d}[y^d]\,
e^{247y}
\prod_{i=0}^{7}
\left(\sum_{h=0}^{d}R_i(t-h)\frac{y^h}{h!}\right)v.
$$

The product is ordered; its factors need not commute.
This identity counts the labelled double-input groups assigned to each missing region and to the no-hole outcome.
Its scalar probability mass is $(247+8)^d/255^d=1$.

For each realized joint hole assignment, replace its regional matrix by $R_i(t-h_i)$.
Entrywise inequalities are preserved by products of nonnegative matrices.
Average only after this replacement, using the joint hole law above.
This argument does not resample the source geometry or assert independence of the regional subsets.

The implementation avoids rational factorials.
Starting with $V_8[0]=v$ and $V_8[k]=0$ for $k>0$, it computes

$$
V_i[k]=\sum_{h=0}^{k}\binom{k}{h}R_i(t-h)V_{i+1}[k-h],
\qquad
Q_{t-d,d}v=\frac1{255^d}
\sum_{k=0}^{d}\binom{d}{k}247^{d-k}V_0[k].
$$

For fixed $t$, the same vectors $V_i[k]$ serve every $d\le t$.
Only after the full mixture does the code define

$$
c_{s,d}=\max\{(Q_{s,d}v)_0,\ (Q_{s,d}v)_1/p\}.
$$

Thus $Q_{s,d}v\le c_{s,d}v$.
One common potential $v$ applies to every column and every category in the source class.
Starting in state zero and using $\mathbf1\le v/p$ introduces the endpoint factor $1/p$ once.

## Constituent support count and union

Fix a cohort with $n_1+n_2=17$ active groups of arities one and two.
Its total number of nonzero constituent words is $w=n_1+2n_2$.
For an RS$[32,16]$ word with a prescribed support of size $j$, the number of full-support words is at most $255^{j-16}$.
Every nonzero constituent has support size at least $17$.
Dropping individual support floors but retaining their sum yields the positive polynomial

$$
P(X)=\sum_{e=0}^{n_1}\sum_{s=0}^{n_2}\sum_{d=0}^{n_2-s}
\binom{n_1}{e}\binom{n_2}{s}\binom{n_2-s}{d}
2^s255^{e+s+2d}c_{e+s,d}X^{e+s+2d}.
$$

The two choices for each partial pair are counted by $2^s$.
The coefficient extraction from $P(X)^{32}$ counts labelled column support arrays.
The aggregate floor is $J\ge17w$ and the support-count denominator is $255^{16w}$.

Let $K_{17}(a,b)$ be the exact number of source sets in ordered class $(a,b)$.
The class/cohort failure probability is at most

$$
\frac{K_{17}(a,b)\binom{17}{n_1}8^{n_1}28^{n_2}}
{255\,p}\;z^{-734003}\;255^{-16w}
\sum_{J\ge17w}[X^J]P(X)^{32}.
$$

The source count, cohort coloring, constituent row-set choices, and projective normalization each appear once.
The code checks $\sum_{a,b}K_{17}(a,b)=\binom{8192}{17}$ before constructing a regional bank.
All constituent pairs obey the same local law, so this bound covers all $28^{n_2}$ fixed pair assignments without a multiplicity partition.

## Directed replay and scope checks

Arb encloses the rational mark, the positive potential, regional caps, hole recurrence, polynomial powers, and final logarithms.
Every maximum compares upper endpoints, not overlapping balls.
The regional bank is built once per precision.
Nine fixed-full-count envelopes supply all $81$ ordered source classes.
Each source class shares one category-cap bank across its $18$ cohorts.

The replay has a 600-second guard and fails without a complete receipt if that guard expires.
It writes a rational upper ceiling on a $10^{-12}$ grid for each of the $1458$ cells.
The verifier requires matching cell, cohort, and aggregate ceilings from fresh 256- and 384-bit replays.
It independently sums $2^{\text{cell ceiling}}$ at 512-bit precision.
Saved FLOAT estimates do not enter either replay or the rounded-cell union.

Seven focused tests check exact labelled hole enumeration, probability normalization, noncommuting order, source-conditioned envelopes, positive support floors, shared bank counts, and fail-closed metadata.
The result concerns only this sparse $r=17$ family in the ideal setup distribution.
It is not a whole-code certificate or a claim about a pseudorandom setup sampler.

The next proof step is the dense-arity family under the same affine mixer, with its own fixed-row-set count and support kernel.

## Replay records

The 256-bit run took $6.4797$ seconds and the 384-bit run took $6.9221$ seconds.
Both raw aggregate ceilings equal the independently recomputed rounded-cell ceiling above.
An independent 1024-bit recombination of all $1458$ rounded cell bounds reproduced the same rational ceiling.
These are proof-computation times, not encoder benchmarks.

| Artifact in `theory/` | SHA-256 |
| --- | --- |
| `hd_two_sum_affine_sparse_certificate.py` | `3fbf95494210ef7ee0aeaf82f7b9682fb386f5a52a3a034d7607adee3f362522` |
| `hd_two_sum_affine_sparse_r17_receipt256.json` | `363247b377a3fa728f2f35bbc4bc9c8b4b098a139406ebb2df5d4b5a53dfad09` |
| `hd_two_sum_affine_sparse_r17_receipt384.json` | `1f86a691d89291ff02981487b00729371ed71bfbfbaf413f1790051ca8647668` |
| `hd_two_sum_affine_sparse_r17_dual_verified.json` | `029b5b870f1c5466e8b8c968790b1fe5983c05292b88b876696937488bb11f74` |

The verifier checks the source pins and exact scope before comparing receipts.
The hashes identify the observed fresh computations; byte hashes alone do not establish their arithmetic correctness.
