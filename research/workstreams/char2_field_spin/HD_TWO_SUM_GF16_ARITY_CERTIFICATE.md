# Directed g16 certificate for one dense group

This replay concerns the g16 two-sum outer map with labels
$C=06+\operatorname{GF}(16)$ inside AES $\operatorname{GF}(256)$.
It uses sixteen RS$[16,8]$ constituents, the alternating independently swapped pair route, and the existing h4 DHD inner map.
The construction identifier is
`two-sum-gf16-coset-g16-rs16x8-h4-swapped-pairs-k20-rate-half`.

The selected message class has 17 active source groups.
Its two ordered matching counts are $(8,0)$.
Sixteen groups have exactly two nonzero constituent words; one group has a fixed nonzero row set of arity 3 through 16.
The replay sums all 65,399 possible dense row sets and all $120^{16}$ sparse pair assignments.
It does not cover other source classes or two dense groups, and it is not a whole-code certificate.

The randomness is the construction's ideal setup distribution: independent nonzero diagonal entries, independent regional pair permutations and fair pair swaps, and the existing independent inner refresh setup.
Fix each source support and each constituent row set before averaging this randomness.
The row sets remain fixed through all 16 outer columns.

## Exact dense comparison

For a fixed dense row set $A\subseteq[16]$, let $N_A(T)$ count uniform-input assignments on $A$ whose mixed output support is $T$.
Zero input coordinates are allowed in this comparison.
The closed kernel in `theory/hd_two_sum_closed_kernel.py` computes these counts exactly.

For arity $a$, its power sum

$$
S_a(w)=\sum_{|A|=a}N_A(T)^{16},\qquad |T|=w,
$$

depends only on $w$.
Define the nonnegative weight

$$
W_a(w)=\left(\frac{S_a(w)}{255^{8a}}\right)^{1/16}.
$$

These weights are not probabilities.
The denominator includes the prescribed-support RS comparison for all $a$ nonzero constituent words.
For any fixed sequence of output supports $T_1,\ldots,T_{16}$, Hölder's inequality gives

$$
\sum_{|A|=a}\frac{\prod_{j=1}^{16}N_A(T_j)}{255^{8a}}
\leq \prod_{j=1}^{16}W_a(|T_j|).
$$

Thus the row-set multiplicity occurs once inside $S_a$.
There is no additional factor $\binom{16}{a}$.
The fixed row set is not resampled between columns.

## Ordered regional coefficient calculation

Let $R_i(m)$ be the directed two-state regional matrix for $m$ active positions in chronological region $i$.
Each matrix uses the source-conditioned maximum over feasible full-pair counts for matching count 8 or 0, according to the region's parity.
The regional averaging includes the uniform pair permutation and fair swaps.
Conditional nonzero amplitudes are independent because of the outer output diagonal.
The resulting envelope holds for every fixed regional subset; regional subsets need not be independent.

Fix a sparse column category $(s,d)$ and put $t=s+d$.
Every active singleton maps to all 16 regions.
Each active double has no hole with probability $239/255$, or a specified one-region hole with probability $1/255$.
These choices are independent across the $d$ labelled double groups.

The common state vector is $v=(1,p)^T$, where $0<p<1$.
For each $t$, compute positive suffix vectors indexed by dense weight $w$ and labelled hole count $k$:

$$
V_i(w,k)=\sum_{b=0}^1\sum_{h=0}^k
\binom{k}{h}R_i(t-h+b)V_{i+1}(w-b,k-h).
$$

The terminal value is $V_{16}(0,0)=v$; other terminal values vanish.
Negative weight indices also vanish.
The implementation constructs suffixes from region 15 to region 0, preserving the product $R_0\cdots R_{15}$.
Its dense bit choices enumerate each labelled output mask once.
Hence its coefficient of weight $w$ already includes the target-mask multiplicity; no factor $\binom{16}{w}$ is added.

For the category $(t-d,d)$, form

$$
U_{t-d,d}(w)=\frac{1}{255^d}\sum_{k=0}^d
\binom{d}{k}239^{d-k}V_0(w,k).
$$

This is the exact positive mixture over labelled sparse holes, after the regional envelope.
For arity $a$, multiply by $W_a(w)$ and sum separately over $w=0$ and $w>0$.
Only then take the largest coordinate ratio against $v$.
Denote the two resulting scalar caps by $c_{a,s,d,0}$ and $c_{a,s,d,1}$.
The same $z$ and $v$ are used for every region, category, and arity in this replay.

## Column floors and total bound

For $e\in\{0,1\}$ define the positive polynomial

$$
P_{a,e}(X)=\sum_{s+d\leq16}
\binom{16}{s}\binom{16-s}{d}2^s255^{s+2d}
c_{a,s,d,e}X^{s+2d}.
$$

The 32 nonzero sparse constituent words each have support at least 9.
The nonzero dense group has at least 9 nonempty columns.
The calculation therefore retains only terms of
$(P_{a,0}(X)+YP_{a,1}(X))^{16}$ with $X$-degree at least 288 and $Y$-degree at least 9.
Keeping only these aggregate floors enlarges the comparison set and preserves the upper bound.

Let $K=15\,708\,305\,994\,709\,246\,428\,578\,365\,440$ be the exact number of source supports in ordered class $(8,0)$.
For each arity, multiply the retained positive coefficient sum by

$$
\frac{K\cdot17\cdot120^{16}\,z^{-734003}}
{255\,p\,255^{256}}.
$$

The factors count the source support, the dense group's placement, and the sparse row pairs once.
The sparse RS comparison contributes $255^{-256}$.
Projectivization contributes $255^{-1}$, and the terminal-vector comparison contributes $p^{-1}$.
The weight marker contributes $z^{-734003}$ for the forbidden output-weight event.
The final bound sums the 14 arity bounds.

## Directed implementation and replay status

`theory/hd_two_sum_gf16_arity_certificate.py` fixes

$$
z=\frac{8970988817957597}{9007199254740992},\qquad
\log_2p=-\frac{2266954815}{33554432}.
$$

It pins the reviewed helpers and the frozen proposal
`hd_two_sum_gf16_arity_dense_8_0_all14_float.json`.
No saved FLOAT estimate enters the directed calculation.
All precision-sensitive operations execute inside one Arb precision context.
Each proposed dyadic root endpoint $u$ is checked with exact integer arithmetic:
$u^{16}255^{8a}\geq S_a(w)$.
The regional maxima and potential quotients use upper endpoints.
The coefficient calculation uses positive Arb polynomial products without numerical cutoffs or FFTs.

Seven focused tests pass. They check independent labelled-mask/hole enumeration with noncommuting matrices, exact multiplicity, root direction, normalization, fail-closed scope, and fresh rounded-union summation.
The adapter has a 600-second resource guard and refuses to overwrite receipts.
Two independent source reviews found no material discrepancy before the fresh serial replays.

| Replay | Arithmetic time | Log₂ upper bound |
| --- | ---: | ---: |
| 256-bit fresh replay | 1.7875 s | $-585542121955663/10^{12}$ |
| 384-bit fresh replay | 1.7154 s | $-585542121955663/10^{12}$ |
| 512-bit sum of rounded arity ceilings | — | $-585542121955663/10^{12}$ |

Every one of the 14 arity ceilings agrees between the two fresh replays.
The independent sum uses those rounded ceilings, not either saved aggregate.
An independent 1024-bit recombination of the 14 rounded ceilings gives the same rational upper bound.
The selected class therefore has failure probability at most
$2^{-585.542121955663}$ for output weight at most 734,003.
This is a bound for the stated ideal random construction and selected message class, not an empirical distance measurement or a whole-code result.

The artifacts are in `theory/`:

| Artifact | SHA-256 |
| --- | --- |
| `hd_two_sum_gf16_arity_certificate.py` | `c47e378f680a290c9b8c2c0fa5e852655aabc22acda47ae6ed9ac571c7eb466b` |
| `hd_two_sum_gf16_arity_certificate_tests.py` | `b8608134142a7d8365eba5f00977fd9f77dceb0d6084731c637fce5b8f292485` |
| `hd_two_sum_gf16_arity_h1_receipt256.json` | `7cd952add67f04e551a2ab4f7b53fee37f54ad7c650fa799dfb20f3825e5ec9e` |
| `hd_two_sum_gf16_arity_h1_receipt384.json` | `a6ab2d3822fb4a56b04f14195de59077decb666fe473451c50a137c781c97776` |
| `hd_two_sum_gf16_arity_h1_dual_verified.json` | `7e0cb6c90ed9f50819d478441c23e809b0b7ccb4c90b26879976d631ca1c7e2f` |

Hashes record provenance; they do not replace arithmetic verification.
The next proof step is a directed two-dense-group replay retaining separate arities: the merged-arity FLOAT control loses substantial margin.
Such work must not enlarge this certificate's scope without a fresh bound.
