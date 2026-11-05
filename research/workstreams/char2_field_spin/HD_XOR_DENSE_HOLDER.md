# Fixed-row-set Hölder compression for mixed XOR groups

The current mixed proof holds each dense constituent-row set fixed across all
32 outer columns. This note gives a positive measure envelope that avoids
enumerating all 219 row sets separately. The 12-type reduction now has a
**56.18-bit directed bound**, reproduced at 256 and 384 bits. It does not
cover two or more dense groups.

The construction is the byte-field XOR2 outer with the existing h4 DHD inner,
alternating pair permutations, and independent fair pair swaps. The control
has 17 active source groups: sixteen two-row groups and one dense group.
Its ordered source matching counts are $(8,0)$. It includes all $28^{16}$
fixed constituent-pair assignments. It does not cover other source classes
or other cohort counts.

## Dominate whole words, not independently resampled row sets

Put $q=256$. A dense group has a fixed active row set
$A\subseteq\{0,\ldots,7\}$ with $a=|A|\ge3$. Let
$H=I+\mathbf1\mathbf1^T$. For an output support $T$, define

$$
N_A(T)=\#\{x\in\mathbb F_q^8:
             \operatorname{supp}(x)\subseteq A,
             \operatorname{supp}(Hx)=T\},\qquad
\nu_A(T)=q^{-a}N_A(T).
$$

Inputs in this comparison include zeros. Conditional on the output support,
the independent nonzero output diagonal gives independent uniform nonzero
amplitudes. The existing RS[32,16] comparison charges

$$
\beta_a=\frac{q^{32a}}{(q-1)^{16a}}.
$$

Thus a 32-column support word $(T_1,\ldots,T_{32})$ receives comparison mass

$$
\beta_a\prod_{j=1}^{32}\nu_A(T_j)
=\frac{\prod_{j=1}^{32}N_A(T_j)}{(q-1)^{16a}}.
$$

Fix a deterministic class $B$ of row sets. Define a nonnegative mask weight

$$
W_B(T)=\left(
 \sum_{A\in B}\frac{N_A(T)^{32}}{(q-1)^{16|A|}}
\right)^{1/32}.
$$

Hölder's inequality gives, for every labelled support word,

$$
\sum_{A\in B}\beta_{|A|}\prod_{j=1}^{32}\nu_A(T_j)
\le \prod_{j=1}^{32}W_B(T_j).
\tag{1}
$$

The left side fixes $A$ through all columns. The right side is an upper
measure, not an equality or a claim that the construction resamples $A$.
The weights $W_B$ need not sum to one. Dense comparison factors and row-set
counts are already inside them; neither factor is charged again.

Inequality (1) is pointwise in the complete support word. It therefore
remains valid after multiplying by any nonnegative weight of that word and
summing. In particular, it preserves an ordered product of regional
matrices. It makes no commutativity assumption.

## Common-potential interface

For one fixed mark $z\in(0,1)$ and one fixed vector $v=(1,p)^T$, with
$0<p\le1$, the existing resource DP bounds each chronological eight-region
operator acting on $v$. Its regional tables retain the ordered source
matching counts. Sparse cancellation branches and dense masks are mixed
before taking the componentwise quotient by $v$.

For a class $B$, weight the empty dense mask by $W_B(\varnothing)$ and each
nonempty mask $T$ by $W_B(T)$. Apply the same potential to both branches.
After this step, the scalar column polynomial can be raised to the 32nd
power. The control retains both support floors:

$$
X\text{-degree}\ge17\cdot32=544,\qquad
Y\text{-degree}\ge17.
$$

The $X$ degree counts sparse nonzero constituent symbols. The $Y$ degree
counts nonempty columns of the dense group. Its floor follows because a
nonzero RS row has at least 17 nonzero columns and the per-column outer map
is invertible. The remaining factor is

$$
\frac{K_{17}(8,0)\,17\,28^{16}}{255\,p\,255^{512}}\,z^{-734003},
$$

where $K_{17}(8,0)=15708305994709246428578365440$ counts source supports.
The factors 17, $28^{16}$, $255^{-1}$, and $p^{-1}$ respectively count the
dense source position, sparse pair assignments, projective normalization,
and the terminal potential bound. Each occurs once.

Each deterministic class may choose its own complete witness $(z,p)$ before
its class bound is evaluated. This is a partition of the message union,
not selection based on sampled setup.

## Bounded controls and the obstruction they distinguish

The first control used all 219 row sets in one Hölder class. It tested five
existing witnesses, without fitting. The best FLOAT result was
$+237.741526$ bits. For each witness, Hölder's bound exceeded the exact
219-row-set sum at that same witness by less than two bits. Those same-witness
sums were themselves positive. This identified substantial loss from
forcing one witness on all row sets.

That observation did **not** establish that Hölder is cheap on every subset.
The next control partitioned row sets by their five selected witnesses.
It retained labelled masks and made no further occupancy relaxation.

| Frozen witness class | Row sets | Exact fixed-$A$ sum, FLOAT | One Hölder class, FLOAT |
|---|---:|---:|---:|
| `a3_mask7_fit` | 6 | −108.155534 | −105.543167 |
| `a4_mask240_fit` | 68 | −61.000916 | −54.471856 |
| `a8_fit` | 104 | −144.711249 | −14.988591 |
| `all_geometry_mask195_fit1` | 35 | −105.840905 | −95.512482 |
| `all_geometry_mask224_fit2` | 6 | −129.649790 | −129.486918 |

The five-class union was only $-14.988590825110746$ bits. The large loss in
the `a8_fit` class was hidden in the earlier all-row-set comparison by other,
larger fixed-$A$ terms. Keeping its 104 row sets exact and compressing the
other four classes gives $-54.471855956455784$ bits. This hybrid calculation
recombines saved values; it performs no new operator evaluation.

The final control split each witness class by arity, without changing any
witness. Exactly 12 nonempty intersections remain. Five operator banks were
rebuilt once; no new mark, potential, or fit was introduced.

| Compression | FLOAT log2 union |
|---|---:|
| Exact 219 fixed row sets, selected witnesses | −61.000916136180436 |
| Five witness classes | −14.988590825110746 |
| Exact `a8_fit` class; four Hölder classes | −54.471855956455784 |
| Twelve witness/arity intersections | **−56.17715193493034** |

The 12-type control loses 4.82376420125 bits relative to the exact sum.
Its dominant type contains the 36 arity-five sets assigned to
`a4_mask240_fit`. The single full-row-set type agrees with its exact result
to FLOAT rounding. The loss reported here includes the Hölder envelope and
the common post-mixture potential quotient; it is not an isolated estimate
of the analytic Hölder inequality's slack.

The separately verified exact-219 bound is
$-61.000916136179$ bits after summing rounded directed cells at 512-bit
precision. The new 12-type receipt verifies the compression; it does not
improve that stronger exact-219 bound.

## Directed replay of the 12 types

The isolated adapter `theory/hd_xor_dense_holder_certificate.py` reads the
frozen 12-type partition and its exact rational witnesses. It ignores all
saved FLOAT bounds. Both full replays rebuild the five regional banks,
ordered resource vectors, mask mixtures, and 12 positive polynomial tails.
The reviewed older helpers and artifacts remain unchanged.

Each moment is computed as an exact rational $N/D$. Arb proposes an upper
endpoint $u$ for its 32nd root. The adapter then checks $u\ge0$ and

$$
u^{32}D\ge N
$$

by exact integer arithmetic on the dyadic representation of $u$.
The function $x\mapsto x^{32}$ is increasing on nonnegative reals, so this
check proves $u\ge W_B(T)$. A failed check aborts the replay. The proof
therefore does not rely on an unverified approximate root. The zero moment
uses the exact endpoint zero.

All later matrix and polynomial expressions have nonnegative coefficients.
The adapter mixes the upper mask weights with both resource-vector
coordinates before taking the common-potential quotient. It includes
$28^{16}$ directly in the prefix, with no further dense comparison factor.
It retains the exact positive tails $X\ge544$ and $Y\ge17$.

The 256- and 384-bit runs took 77.64 and 80.46 seconds. All 12 rational
type ceilings agree exactly, as does the raw aggregate ceiling

$$
-\frac{5617715193493}{100000000000}.
$$

A separate 512-bit calculation sums $2^{b_i}$ for the 12 saved rational
type ceilings $b_i$. It does not reuse the cached aggregate. Its ceiling is

$$
\log_2\Pr[\text{a bad scalar codeword in the stated class}]
\le -\frac{56177151934929}{1000000000000}
=-56.177151934929.
$$

Here “bad” means nonzero output weight at most 734003 in length $2^{21}$.
Probability is over the ideal independent setup of the stated construction.
The class is precisely source $(8,0)$, cohort $(0,16,1)$, with all fixed
sparse constituent-pair assignments. This is not a whole-code certificate
or a new certificate for the proposed XOR state refresh.

## Interface for the next small-cohort experiment

Index the 12 deterministic types by $B_1,\ldots,B_{12}$. For $h$ dense
groups, retain a type histogram $(h_1,\ldots,h_{12})$ with
$\sum_i h_i=h$. There are $\binom{h+11}{11}$ such histograms: 78 when
$h=2$, and 364 when $h=3$. Charge the assignment factor
$h!/\prod_i h_i!$ once. A histogram bound must be uniform over assignments
of those types to the fixed dense source positions.

For each histogram, choose **one global mark and one global potential**.
The successful single-group witnesses cannot be assigned independently to
interleaved groups. All groups contribute to the same emitted word and
share the same chronological inner-state product. Group-specific marks or
potentials would require a separate proved interface.

For small $h$, an exact intermediate option is to convolve whole labelled
mask weights. At $h=2$, group the $256^2$ mask pairs by the ordered regional
occupancy vector and the nonempty-column indicators. This reduces repeated
matrix evaluations without permuting regions. Its cost and bound quality
remain unmeasured.

A cheaper positive envelope for larger $h$ uses, for each type,

$$
D_i(y,x,u)=\sum_{T\subseteq[8]}W_{B_i}(T)
 y^{\mathbf1[T\ne\varnothing]}x^{|T|}
 u^{\mathbf1[T=[8]]}.
$$

The product $\prod_iD_i^{h_i}$ records the number $e$ of nonempty dense
masks, their total occupancy $t$, and their number $f$ of full masks.
It has $O(h^3)$ possible degree triples. Generally $W_{B_i}(T)$ depends on
the labelled mask, not just its weight; the displayed sum must use those
actual weights.

If a sparse column category has $s$ singleton groups, $d$ double groups,
and $k$ cancellations, every combined regional occupancy satisfies

$$
d-k+f\le m_j\le s+d+e,\qquad
\sum_{j=0}^{7}m_j=7s+8d-6k+t.
$$

The ordered box/sum resource DP can bound all such vectors. This adds a
separate relaxation of dense-mask placement; it was not used in the
12-type control. Mix its positive feature weights and sparse cancellation
probabilities before applying the common-potential quotient.

For general sparse cohort sizes $(n_1,n_2)$, let $c_{s,d,e}$ denote the
resulting scalar caps after the dense-feature mixture. A positive column
polynomial is

$$
P(X,Y)=\sum_{a,s,d,e}
 \binom{n_1}{a}\binom{n_2}{s}\binom{n_2-s}{d}2^s
 255^{a+s+2d}c_{a+s,d,e}X^{a+s+2d}Y^e.
$$

Here $a$ counts active one-row groups in the column; $s$ and $d$ count
two-row groups contributing one and two nonzero symbols, respectively.
In $P^{32}$ retain $X\ge17(n_1+2n_2)$ and $Y\ge17h$, then divide by
$255^{16(n_1+2n_2)}$. The second floor relaxes the separate 17-column
floor for each dense group. It is a necessary condition, not an equivalent
one. Source supports, source colors, sparse row choices, and sparse pair
assignments remain outside and must be charged exactly once. No extra
dense row-set factor is added.

## Exact/coarse hybrid and stopping point

With the 12-type $h=1$ replay complete, the next bounded step is to test
selected $h=2$ histograms using a small, fixed global witness bank per
histogram. Retain labelled occupancy vectors: the companion control in
`hd_xor_dense_two_feature.py` found the coarse feature envelope already
loses about 467 bits on the same dominant $h=1$ type. Its first $h=2$
control also fails. These failures concern the coarse proof envelope, not
the encoder. Do not infer useful larger-$h$ bounds merely from the
polynomial number of types.

At the opposite end, the route-independent inverse-image argument in
[HD_XOR_COVERAGE_PLAN.md](HD_XOR_COVERAGE_PLAN.md) covers total active
constituent count $A_{\rm tot}\ge61191$. Its directed log2 ceiling is
$-20721908276069/125000000000=-165.775266208552$.
The threshold counts constituent words, not active source groups. That
tail avoids any row-set enumeration. It leaves a genuine middle range
between small exact/type-compressed cohorts and the terminal bound.

No native code, encoder ensemble, or benchmark was changed by this work.

## Reproduction and pins

All files below are in `theory/` beside this note. The helpers are
`hd_xor_dense_holder.py`, `hd_xor_dense_holder_classes.py`, and
`hd_xor_dense_holder_arity.py`, each with an isolated `_tests.py` file.
The tests cover exact GF4 mask counts, every three-column Hölder inequality
in a tiny model, non-symmetric classes, zero weights, class partitions,
and mixture-before-quotient ordering. No production replay is hidden in
those tests.

The frozen proposal is
`hd_xor_mixed_all_geometry_8_0_refine64.json`, SHA256
`44fd35d8bf2e78391572a14d071cd283dc547f7bc204b3802c180132eaaab4e6`.
It stores the exact five selected rational witnesses and all 219 class
assignments. The controls reject a changed proposal.

| Diagnostic artifact | SHA256 |
|---|---|
| `hd_xor_dense_holder_8_0_h1_control5.json` | `dd340a4fe73529e49ac0f34a4a2ec394d9b0c4c78d5376e1b55d37603e0433df` |
| `hd_xor_dense_holder_8_0_h1_classes5.json` | `46365397077b1ed038f3b258c4cd438cc5a455ede2771d7a4c4688974d8e8484` |
| `hd_xor_dense_holder_8_0_h1_arity12.json` | `4577166d2aefbe1912ac6bc096c8f08a3604439c03f72ff0a2dd7ad4fef74a0e` |

The exact-219 verified artifact is
`hd_xor_mixed_all_geometry_8_0_dual_verified.json`, SHA256
`af53e419971b68619c23489a71090c09a1d706b13288d833c1c1bb315178a095`.
The directed adapter was independently reviewed, including fixed-A witness
selection, both positive floors, and the once-only $28^{15}$ transfer from
the common-pair count. Its five focused tests pass. Byte hashes identify
the observed artifacts; they are not substitutes for the directed arithmetic.

The new adapter was independently reviewed before either full replay.
Its SHA256 is
`0934613d5690c280f81d94ecf41f56f25b1a8bd020f532d42b264458a701428e`.
Eight focused adapter tests cover exact root direction, rational mixtures,
both support floors, fail-closed scope checks, five-bank reuse, once-only
prefix factors, and the observed dual receipts. Together with the 11
earlier Hölder tests, all 19 tests pass.

| Directed artifact | SHA256 |
|---|---|
| `hd_xor_dense_holder_8_0_h1_receipt256.json` | `33ad1a1b1550b640497a97b1d8e68101fe9bd618c8868bb2579e9bbb76f82ced` |
| `hd_xor_dense_holder_8_0_h1_receipt384.json` | `b264538ea3e0ba2df0712b7d21e65dc4ed7af86bb960c0cb3648ce746ec4813c` |
| `hd_xor_dense_holder_8_0_h1_dual_verified.json` | `a823fddd3d6e6e92cadce141f43fc06e8ae30cb44f403db1533ae845cca7f54c` |
