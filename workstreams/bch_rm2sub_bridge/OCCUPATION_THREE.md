# Occupation three without a cubic coefficient table

Occupation two left open messages with three or more nonzero BCH rows.
For three rows, auxiliary Bernoulli supports give a compact bound after
paying exact conditioning costs. The construction itself is unchanged.

## Model and result

Use the fixed BCH [256,128] code, t128_s15 map, and setup in `README.md`.
There are L = 8192 rows and n = 256 regions. The message length is 2^20,
the output length is 2^21, and the bad-weight cutoff is H = 209716.
The state starts at zero, output precedes the update, and each epoch uses
an independent nonzero field multiplier. All setup randomness is shared
across messages, as before.

Let Z_3 count messages with exactly three nonzero rows and output weight at
most H. The checked rational upper bound U_3 satisfies

\[
\Pr[Z_3>0]\le\mathbb E[Z_3]\le U_3<2^{-126}.
\]

Its diagnostic margin is 126.5064403982 bits. Together with the retained
occupation-one and occupation-two bounds,

\[
\Pr[Z_1+Z_2+Z_3>0]\le U_1+U_2+U_3<2^{-49}.
\]

No bound for occupations four through 8192 is claimed here. The probability
space remains the ideal fresh-multiplier setup, not a seeded implementation.

## Three input bits in an epoch

The exact map audit gives minimum kernel weight four for B. Consequently,
every nonzero input of weight at most three has nonzero syndrome.
Retain the state classes Z,D,L, M = 2^15 - 1, kappa = M/(M-1), and d_A = 48.
For j in {0,1,2,3}, define

\[
m_j(z):=\frac1{M\binom tj}
\sum_w a_w\sum_v\binom wv\binom{t-w}{j-v}z^{w+j-2v},
\]

where a_w is the nonzero A-spectrum and impossible binomial terms vanish.
The weight-three transfer is

\[
T_3(z)=\begin{pmatrix}
0&z^3&0\\
z^{d_A-3}/M&0&z^{d_A-3}\\
\kappa m_3(z)/M&0&\kappa m_3(z)
\end{pmatrix}.
\]

The same weighted-measure argument as for T_1 and T_2 applies. Zero activates
an arbitrary nonzero state. A D state emits at least d_A - 3 bits. An L state
has moment at most kappa m_3. Conditioned on the entering state and input,
fresh multiplication gives termination probability 1/M, independent of the
output weight. Conditional on survival, the next law belongs to L.

For E = L/t, the region envelopes are

\[
R_j(z):=\frac1{\binom Lj}[u^j]
\left(\sum_{a=0}^3\binom ta u^aT_a(z)\right)^E,
\qquad 0\le j\le3.
\]

These coefficients count all placements of up to three distinct positions
in a region. They include all epoch collisions without allowing two input
bits to occupy the same position. Independent region permutations give this
uniform subset law for each fixed triple of outer supports.

## Auxiliary supports and exact conditioning

Fix three nonzero outer words of weights w_1,w_2,w_3. Their independently
permuted supports are independent uniform subsets of those sizes.
We now define a separate experiment solely to bound its moment.

In the auxiliary experiment, row i has independent Bernoulli(p_i) entries
across the n regions. The three rows are independent. Then apply the same
random region permutations and fresh-multiplier inner. Let W_i count its
row-i ones and V be the resulting output weight.

For 0 < p_i < 1, write

\[
\beta_p(w):=\binom nw p^w(1-p)^{n-w}.
\]

Conditioning on W_i = w_i for all three rows recovers exactly the fixed-weight
support law. Each conditioning probability is positive. Therefore

\[
\begin{aligned}
\mathbb E[z^V\mid W_i=w_i\ (i=1,2,3)]
&=\frac{\mathbb E[z^V\mathbf1_{\{W_i=w_i\ \forall i\}}]}
{\prod_i\beta_{p_i}(w_i)}\\
&\le\frac{\mathbb E[z^V]}{\prod_i\beta_{p_i}(w_i)}.
\end{aligned}
\]

The inequality uses nonnegativity, not an approximation to the conditioning
event. The special choice p_i = 1 is used only for the singleton weight n;
that row is deterministic and its conditioning probability is one.

Let e_Z = (1,0,0), and let 1 be the three-entry column of ones. Independent
auxiliary coordinates allow a single region matrix

\[
M(p_1,p_2,p_3;z):=
\sum_{b\in\{0,1\}^3}
\left(\prod_i p_i^{b_i}(1-p_i)^{1-b_i}\right)R_{b_1+b_2+b_3}(z).
\]

The state persists across regions, so

\[
\mathbb E[z^V]\le e_ZM(p_1,p_2,p_3;z)^n\mathbf1.
\]

Combining this inequality with Chernoff's inequality bounds the bad-output
probability for each fixed outer-weight triple. No input-one deletion or
expected-weight monotonicity is required.

## Weight groups and counting

Let U_w be the retained deterministic BCH shell caps. Partition all supported
nonzero weights into the following 13 disjoint groups:

```
38..42, 44..48, 50..54, 56..64, 66..78, 80..100,
102..128, 130..154, 156..176, 178..198, 200..210,
212..218, {256}.
```

Every range advances by two. For a group G, define its conditioning cost

\[
C_G(p):=\sum_{w\in G}\frac{U_w}{\beta_p(w)}.
\]

For each ordered group triple (G_1,G_2,G_3), any admissible p_1,p_2,p_3 and
z in (0,1) give the contribution bound

\[
\binom L3 z^{-H}
e_ZM(p_1,p_2,p_3;z)^n\mathbf1
\prod_{i=1}^3 C_{G_i}(p_i).
\]

Here the three chosen row positions are ordered by index. Their local message
values may coincide. Summing over the group partition covers every nonzero
message with occupation three exactly once before applying the upper bounds.

The expression is invariant under simultaneous permutations of the groups
and their probabilities. We evaluate the 455 unordered group triples and
multiply by 1, 3, or 6, according to their number of distinct orderings.
Different triples may use different auxiliary distributions and tilts;
each is a separate valid upper bound for the original setup.

## Search, arithmetic, and independent checks

Choosing each p from the group's central weight produced a vacuous bound.
That diagnostic is retained in `t128_s15_q3_conditioning_outward.json`.
It is not a successful distance certificate.

The successful search shifts the auxiliary log-odds away from that initial
choice, separately for each group triple and tilt. The tilt grid is
z = exp(-exp(j/10)), with j = -85,-80,...,-55. Search is binary64 only.
Its output probabilities are subsequently treated as exact dyadic rationals.
No global optimality claim is needed: any admissible choices prove a bound.

The compact certificate recomputes all 455 cases with 256-bit Arb. It computes
conditioning costs as exact rationals, then rounds them upward to dyadics
with 192 bits of relative precision. This avoids growing rational denominators
without weakening the inequality's direction. The final aggregation is exact.

An independent 512-bit verifier checks every case. It uses binary polynomial
powering for regions, explicit enumeration of the eight auxiliary support
patterns, and 256 successive state-vector updates. The producer instead uses
sequential region coefficients and binary powering of the mixed region matrix.
Every independently recomputed moment lies below its stored upper bound.
The verifier also checks every conditioning cost, group multiplicity, source
hash, contribution, and the exact partial sum U_1 + U_2 + U_3.

Additional exact tests cover 2,016 GF(16) weighted-prefix comparisons with
epoch weights zero through three. Small-region enumeration checks every
placement of up to three ones. Enumeration of all support triples on three
regions checks the auxiliary-mixture identity and its conditional inequalities.

The current receipt is `generated/t128_s15_q3_compact_outward.json`.
It contains scalar bounds and witnesses, not a three-dimensional experiment
table. Generated artifacts remain ignored by Git. From this worktree root:

```powershell
python -B workstreams/bch_rm2sub_bridge/test_occupation_three.py
python -B workstreams/bch_rm2sub_bridge/verify_q3_compact.py
```

For fresh generation, use `optimized_conditioning.py screen`, followed by
`certify_q3_compact.py`. Both refuse to overwrite their output files. The
unrounded experimental aggregation in `optimized_conditioning.py certify`
is not the compact production path and is unnecessary for replay.

## Next step

Occupation four introduces genuine zero-syndrome inputs: the selected kernel
has weight-four words. The next transfer must account for those inputs in
both activation and live-state transitions. The auxiliary-conditioning method
can then be tested as a compact route to additional occupations. Its success
for occupation three does not establish a bound for the remaining tail.

All changes remain in our integration directory. The other worktree and the
completed Q1/Q2 source files and certificates were not modified.
