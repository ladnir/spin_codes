# Directed small-cohort coverage for the g16 outer map

This adapter extends the selected-source calculation to every source support, while keeping its message-class scope explicit.
The construction is `two-sum-gf16-coset-g16-rs16x8-h4-swapped-pairs-k20-rate-half`.
It uses sixteen RS$[16,8]$ constituents per group, the $06+\operatorname{GF}(16)$ two-sum mixer, independent output diagonals, alternating fair-swapped pair routing, and the existing h4 DHD inner map.

Write $(n_1,n_2,h)$ for the numbers of active groups containing one, two, or at least three nonzero constituent words.
The adapter covers the following 322 disjoint cells:

- every $(n_1,n_2,h)$ with $2\leq r=n_1+n_2+h\leq17$ and $h\in\{0,1\}$: 320 cells;
- the two sparse cells at $r=1$: $(1,0,0)$ and $(0,1,0)$.

The cell $(0,0,1)$ at $r=1$ is excluded.
Cells with two or more dense groups and weights above 17 are also excluded.
This is not a whole-code certificate.

## Uniform source envelope and fixed row sets

Fix a source support $G$ of size $r$ and fix each group's constituent row set before averaging setup randomness.
For each outer column, the actual regional subsets may be correlated.
Let $R_{m,f}$ be the directed two-state regional matrix for a fixed regional subset of size $m$ with $f$ full pairs.
It averages the regional pair permutation, fair swaps, independent nonzero amplitudes, and inner setup.
Define the entrywise envelope

$$
\overline R_m=\max_{0\leq f\leq\lfloor m/2\rfloor}R_{m,f}.
$$

Every actual regional matrix is bounded by $\overline R_m$, regardless of $G$'s matching counts.
Nonnegative products preserve this order.
The comparison therefore applies to every fixed $G$, without an independence assumption on regional subsets.
The final source count is $\binom{8192}{r}$.

For a dense row set $A\subseteq[16]$, let $N_A(T)$ be its exact uniform-input output-support count, including zero input coordinates.
The closed kernel gives an exact numerator $S_a(w)=\sum_{|A|=a}N_A(T)^{16}$ for every fixed $T$ of weight $w$.
The merged weight is

$$
W_H(w)=\left(\sum_{a=3}^{16}\frac{S_a(w)}{255^{8a}}\right)^{1/16}.
$$

The implementation sums these rational moments over the common denominator $255^{128}$ before taking a checked 16th root.
Hölder's inequality bounds the sum over the 65,399 fixed dense row sets by a product of these column weights.
The same row set remains fixed through all 16 columns.
There is no additional dense row-set count or density factor outside $W_H$.
Merging arities is a relaxation; it is used only in the cells freshly replayed by this adapter.

## Positive column calculation

The recurrence and sparse-hole law are the ones proved in `HD_TWO_SUM_GF16_ARITY_CERTIFICATE.md`.
For each fixed witness, one ordered coefficient calculation supplies all categories up to 17 active sparse groups, and one dense-weight calculation supplies categories up to 16.
Using a larger precomputed category table does not invent additional active groups: each cell reads only its actual category entries.

Let $c_{s,d,e}$ be the positive common-potential cap after the exact sparse-hole mixture and dense-weight sum, with $e=0,1$ marking an empty or nonempty dense column.
The regional matrices act in chronological order.
The quotient against the common vector $v=(1,p)^T$ is taken only after both mixtures.
For a sparse cohort $(n_1,n_2)$, its column polynomial is

$$
P_e(X)=\sum_{u=0}^{n_1}\sum_{s=0}^{n_2}\sum_{d=0}^{n_2-s}
\binom{n_1}{u}\binom{n_2}{s}\binom{n_2-s}{d}
2^s255^{u+s+2d}c_{u+s,d,e}X^{u+s+2d}.
$$

Here $u$ selects active singleton rows, $s$ selects one active row of a pair, and $d$ selects both rows.
The factor $2^s$ counts which row is active in each partial pair.
For $h=0$, use the analogous polynomial without the dense marker.

Set $J=n_1+2n_2$.
The comparison retains $X$-degree at least $9J$ in the 16-column power and divides by $255^{8J}$.
When $h=1$, it also retains $Y$-degree at least 9 in $(P_0+YP_1)^{16}$.
These aggregate floors retain every actual nonzero RS row and dense group, although they may admit additional comparison configurations.

For each cell, the retained coefficient sum is multiplied by

$$
\binom{8192}{r}\binom rh\binom{r-h}{n_1}
16^{n_1}120^{n_2}\frac{z^{-734003}}{255p}.
$$

The source support, cohort placement, and sparse row sets are counted once.
The factors $255^{-1}$ and $p^{-1}$ are respectively the projective and terminal-vector normalizations.
The sum over the 322 cells is a union over disjoint message classes.

## Frozen choices and numerical checks

`theory/hd_two_sum_gf16_small_certificate.py` pins the proposal
`hd_two_sum_gf16_small_cover_r1_17_allsource_float.json`.
Its saved assignment chooses `dense-inherited` for five included cells and `h2-inherited` for 317 cells.
The unused `sparse-inherited` witness remains part of the pinned proposal, but its operator bank is not rebuilt.
Each chosen pair $(z,p)$ is used globally within its entire cell.
No witness depends on sampled setup, and the directed replay does not reselect witnesses using numerical bounds.

Before implementing the replay, all regional FLOAT entries were compared with 256-bit Arb entries at all three proposal marks.
Each mark has 90 $(m,f)$ matrices and 360 entries.
The four exact zeros at each mark matched the FLOAT negative infinities.
No positive-entry underflow, invalid sign, NaN, or zero mismatch occurred.

| Witness | Smallest positive Arb log₂ entry | FLOAT minus Arb log₂ range |
| --- | ---: | ---: |
| sparse-inherited | −158.15385219340328 | $[-1.279\cdot10^{-13},\ 2.843\cdot10^{-14}]$ |
| dense-inherited | −159.88981284318183 | $[-2.843\cdot10^{-14},\ 2.843\cdot10^{-14}]$ |
| h2-inherited | −178.89660711957174 | $[-8.527\cdot10^{-14},\ 2.843\cdot10^{-14}]$ |

This comparison is a diagnostic, not the certificate.
The directed adapter rebuilds the regional matrices, entrywise maxima, positive coefficient powers, and selected-cell bounds from exact rational literals.
Each root endpoint is checked by an exact integer power comparison.
No FLOAT estimate is used in its arithmetic.

Seven focused tests pass, covering the exact partition, rational merged moments, upper-endpoint envelopes, count and floor factors, shared bank use, scope rejection, and deadline guards.
Two independent source reviews cleared the stable adapter before the fresh serial replays.

| Replay | Arithmetic time | Log₂ upper bound |
| --- | ---: | ---: |
| 256-bit fresh replay | 4.0805 s | $-563328194368259/10^{12}$ |
| 384-bit fresh replay | 4.3457 s | $-563328194368259/10^{12}$ |
| 512-bit sum of rounded cell ceilings | — | $-281664097184129/(5\cdot10^{11})$ |

All 322 cell ceilings, all 17 active-weight ceilings, and the raw unions agree between precisions.
Across the 322 cells, the saved FLOAT bound minus the fresh directed ceiling lies between $-1.865\cdot10^{-11}$ and $4.093\cdot10^{-12}$ bits.
The directed replay therefore reproduces the screen's nonmonotonic behavior in $r$; regional underflow does not explain it.
The independent sum uses the rounded cell ceilings, not the saved raw union.
It gives failure probability at most $2^{-563.328194368258}$ for the selected classes and forbidden output weight at most 734,003.
This probability concerns the ideal construction's setup distribution; it is not an empirical distance measurement.

The artifacts are in `theory/`:

| Artifact | SHA-256 |
| --- | --- |
| `hd_two_sum_gf16_small_certificate.py` | `8102cd1431bb05a425a0fad5ad9f271f49e8cdca012d0665d1dfa8a984a1c546` |
| `hd_two_sum_gf16_small_certificate_tests.py` | `4249a345f49afcbd565835d267c05ff53e9cfd52a6c7e94e74a45080a32f5492` |
| `hd_two_sum_gf16_small_receipt256.json` | `e0d8e00016868233557c468a5771b594044c10e5770468185b21fd392111f63f` |
| `hd_two_sum_gf16_small_receipt384.json` | `6ac9b824b6bd8c5dd6bf550c64329484d7db2a3deed2c25fa0612eca6f1b8efa` |
| `hd_two_sum_gf16_small_dual_verified.json` | `0a7e0c131c52f68f912fd100595e54e5740f9065365acb695f135632b0c9d55e` |

The next step is to treat the excluded $r=1$ dense class with separate arities and full column matrices.
Any later combination must retain the separate construction and scope metadata.
