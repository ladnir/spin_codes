# Goal 08: C38 information-set preflight

## Result

The 38-node target passes the information-set preflight. Every zero-anchor
kernel has at least 35 pairwise disjoint information sets. The result holds
for all 38 anchor positions.

The initial, centered, and terminal kernels were also optimized exactly. Their
maximum packing numbers are 36, 35, and 36. The exact successor tests rule
out one additional information set in each case.

This result gives a uniform finite reduction for the proposed bound

\[
D_{64}(38)\ge228.
\tag{A38}
\]

If a state violates (A38), then one of 1,330 systematic charts represents
that state by two 64-bit values of weights at most five and six. The first
value is one output block. The second value is the restriction to an
information set.

The reduction is materially stronger than the 24-node reduction. In
particular, the restriction radius is six at every anchor. The centered and
terminal 24-node shortenings required radius seven.

The preflight does not prove (A38). A naive enumeration of the two sparse
64-bit values remains too large. The next proof step needs a joint decoder or
an exact low-distance search in the resulting systematic charts.

## Zero-anchor kernels

Let

\[
Q_{38}(s):=(q_0(s),\ldots,q_{37}(s))
\in\mathbb F_2^{2432}
\]

be the 38-node observation word of a lifted state
$s\in\mathbb F_2^{128}$. For an anchor position
$j\in\{0,\ldots,37\}$, define

\[
K_j:=\{s\in\mathbb F_2^{128}:q_j(s)=0\}.
\]

The map $s\mapsto q_j(s)$ has rank 64. Therefore $K_j$ has dimension 64.
Delete the anchor block from $Q_{38}(K_j)$. The resulting shortened code
has parameters

\[
[2368,64].
\]

An information set is a set of 64 output coordinates whose restriction map
is injective on this shortened code. Since both spaces have dimension 64,
the restriction map is then bijective.

An exact linear-matroid union algorithm packed 35 disjoint information sets
in every $K_j$. An independent audit reconstructed every kernel and checked
the rank of all 1,330 information sets.

## Exact representative packing numbers

The preflight first considered three representative anchors.

| Anchor position | Maximum information sets | Union size for one more copy | Required size |
|---:|---:|---:|---:|
| 0 | 36 | 2,347 | 2,368 |
| 19 | 35 | 2,290 | 2,304 |
| 37 | 36 | 2,312 | 2,368 |

For example, the initial shortening contains 36 disjoint bases. In the union
of 37 copies of its represented matroid, the maximum independent union has
size 2,347. Since 37 bases would require size 2,368, no such packing exists.

The uniform sweep then packed 35 sets at every remaining anchor. It proves
feasibility, not optimality, away from the three representative positions.

## Uniform sparse-chart reduction

Suppose a nonzero state $s$ violates (A38). Then

\[
\operatorname{wt}(Q_{38}(s))\le227.
\]

Some anchor position $j$ satisfies

\[
\operatorname{wt}(q_j(s))
\le
\left\lfloor\frac{227}{38}\right\rfloor
=5.
\]

Fix 35 disjoint information sets
$I_{j,1},\ldots,I_{j,35}$ for $K_j$. These coordinate sets exclude the
anchor block. Their restrictions are disjoint subsets of $Q_{38}(s)$.
Consequently,

\[
\sum_{i=1}^{35}
\operatorname{wt}\bigl(Q_{38}(s)|_{I_{j,i}}\bigr)
\le
227-\operatorname{wt}(q_j(s)).
\]

For every anchor weight $a\in\{0,\ldots,5\}$,

\[
\left\lfloor\frac{227-a}{35}\right\rfloor=6.
\]

Therefore some $i\in\{1,\ldots,35\}$ satisfies

\[
\operatorname{wt}\bigl(Q_{38}(s)|_{I_{j,i}}\bigr)\le6.
\]

For fixed $j$ and $i$, define

\[
\Phi_{j,i}(s)
:=
\left(q_j(s),Q_{38}(s)|_{I_{j,i}}\right)
\in\mathbb F_2^{128}.
\]

The map $\Phi_{j,i}$ is bijective. Its first component determines the affine
fiber of the anchor map. Within that fiber, the second component determines
the unique element of the kernel $K_j$.

It follows that every counterexample lies in the exact cover

\[
\bigcup_{j=0}^{37}
\bigcup_{i=1}^{35}
\left\{
s:
\operatorname{wt}(q_j(s))\le5,
\quad
\operatorname{wt}(Q_{38}(s)|_{I_{j,i}})\le6
\right\}.
\tag{SC38}
\]

There are $38\cdot35=1{,}330$ systematic charts in (SC38).

## Comparison with 24 nodes

The systematic zero-anchor certificate uses the disjoint sets as follows.
If a word has weight at most $B$, let

\[
p:=\left\lfloor\frac{B}{m}\right\rfloor,
\qquad
r:=B-mp,
\]

where $m$ is the number of disjoint sets. Enumerate restriction weights
zero through $p-1$ in all $m$ sets. Then enumerate weight $p$ in
$r+1$ fixed sets.

The representative comparisons are:

| Position | C24 sets | C24 radius | C24 candidates | C38 sets | C38 radius | C38 candidates |
|:---|---:|---:|---:|---:|---:|---:|
| Initial | 22 | 6 | 1,082,372,342 | 36 | 6 | 1,198,623,204 |
| Center | 19 | 7 | 8,415,660,131 | 35 | 6 | 1,640,165,779 |
| Terminal | 20 | 7 | 4,150,424,788 | 36 | 6 | 1,198,623,204 |

Across these three positions, the candidate count falls from
13,648,457,261 to 4,037,412,187. The C38 count is 29.6 percent of the C24
count.

The initial position becomes 10.7 percent more expensive. The centered and
terminal positions become 5.13 and 3.46 times cheaper. More importantly, all
38 C38 anchors have radius six.

The number of 64-bit vectors of weight at most six is

\[
\sum_{w=0}^{6}\binom{64}{w}=83{,}278{,}001.
\]

The corresponding radius-seven list has size 704,494,193. Reducing the
radius from seven to six shrinks this list by a factor of 8.46.

## Remaining joint problem

The cover (SC38) is finite and uniform, but direct enumeration is still
impractical. One chart contains at most

\[
\left(\sum_{a=0}^{5}\binom{64}{a}\right)
\left(\sum_{w=0}^{6}\binom{64}{w}\right)
=691{,}509{,}957{,}277{,}633
\]

candidate systematic vectors. Across all 1,330 charts, the naive count is
approximately $9.20\cdot10^{17}$.

The two sparse parts cannot be processed independently. For one chart, write
the systematic generator as

\[
G_{j,i}(u,v)=G^{(1)}_{j,i}(u)+G^{(2)}_{j,i}(v),
\]

where $\operatorname{wt}(u)\le5$ and
$\operatorname{wt}(v)\le6$. The remaining task is to rule out

\[
\operatorname{wt}(G_{j,i}(u,v))\le227.
\]

This is a structured bichromatic closest-pair problem. A useful next method
must combine the two lists before their Cartesian product forms. Candidate
methods include Stern-style information-set decoding, a meet-in-the-middle
search with exact bucket bounds, or a solver encoding that uses the 128
systematic coordinates exposed by $\Phi_{j,i}$.

## Refutation search

A deterministic diagnostic search used 10,000 pseudorandom restarts, all 128
basis states, the known 18-node transient, coordinate descent, and a final
pair descent. It found no counterexample to (A38).

The best word has weight 720. It is the 38-node extension of the known
18-node transient. This result is positive diagnostic evidence only. Local
search supplies no lower bound and cannot exclude a weight-227 word.

## Assessment

The 38-node route is materially better than the 24-node route. The preflight
removes the varying restriction radius that made the centered and terminal
C24 slices expensive. It also converts every possible counterexample into
the same two-part sparse form.

The proof is still algorithmically open. The remaining obstruction is now a
specific 128-coordinate joint decoding problem rather than an unstructured
collection of affine cosets.

The next goal should implement one joint-chart decoder on the width-16 BCH
instance. The complete width-16 census from Goal 07 provides exact minima for
every anchor position and anchor weight. A correct decoder must reproduce
those minima before it is trusted on the width-64 C38 charts.

## Evidence and scope

The exact maximum-packing source is
`scripts/prepare_riffle_packetmul_wrapmul_2lap_g4_goal08_c38_information_sets.py`.
The all-anchor feasibility source is
`scripts/prepare_riffle_packetmul_wrapmul_2lap_g4_goal08_c38_uniform35.py`.
The independent audit is
`scripts/audit_riffle_packetmul_wrapmul_2lap_g4_goal08_c38_preflight.py`.

The authenticated summary is
`constructions/riffle_packetmul_wrapmul_2lap_g4/receipts/goal08_c38_preflight_audit.json`.
The 38 all-anchor receipts and three exact maximum receipts are in the same
receipt directory.

The audit verifies the rank and disjointness of every information set. The
maximum-packing claims additionally rely on the exact augmenting-path
implementation of linear-matroid union. No artifact in this goal proves
(A38).
