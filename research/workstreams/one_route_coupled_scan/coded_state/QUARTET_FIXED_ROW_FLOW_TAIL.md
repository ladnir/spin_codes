# A concrete moment target for the fixed-row flow

2026-11-01. The exact flow interface controls every state mask for a fixed
eligible row set. Its remaining task is a probability bound on flow
deficits. This note states one sufficient window moment and explains why
generic concentration from a small observed deficit is inadequate.

Use the network and window construction in
[the fixed-row flow note](QUARTET_FIXED_ROW_FLOW.md). Fix `H` with
`|H|=4096` before sampling the routing and quartet partitions. Scale the
network capacities by two. A 64-update window has `32768` source units;
let `delta` be its integer maximum-flow deficit. The physical deficit is
`tau=delta/2`.

The 64 disjoint windows have independent, identically distributed
deficits `delta_1,...,delta_64`, with exact regional half-row quotas.
Their flows compose, so the global integer deficit is at most their sum.
A structural obstruction with this eligible-row count requires global
integer deficit at least `104860`.

## A sufficient moment, including the entire row-set union

It is enough to prove

```
E[2^(delta/8)] <= 2^75.                              (1)
```

Indeed, Markov's inequality, window independence, and one row-set union
give

```
Pr[some structural obstruction with eligible count4096]
 <= binom(8192,4096) * 2^(-104860/8)
                     * (E[2^(delta/8)])^64
 < 2^(8192-13107.5+64*75)
 = 2^(-115.5).                                      (2)
```

This bound requires no union over state masks or profile coordinates.
It concerns the whole eligible-row class, but other eligible-row counts
still require their own bounds.

One sufficient decomposition of (1) is

```
E delta <= 512,
log2 E[2^((delta-E delta)/8)] <= 11.                  (3)
```

These are targets, not established estimates. A variance estimate alone
would not prove the second inequality. More generally, a constructive
feasible flow with integer loss `G>=delta` can replace `delta` in (1).
This avoids analyzing the exact maximum flow if `G` has a tractable
moment.

The unavoidable terminal-source component is already part of `delta`.
It must not be charged again in (1) or (2). The exact 64-update diagnostic
returned `delta=467`; this is one realization, not a mean or a moment
bound. The existing backward greedy construction returned loss `3021`
on that instance. Neither sample establishes an exponential moment.
The next branch seeks a sharper constructive flow or a direct
deficit-tail argument meeting (1).

## Why generic bounded-range concentration is insufficient

Sending every non-`H` source to its own input gives a feasible flow of
`16384` units in each balanced window. Hence `0<=delta<=16384`.
The generic bounded-range inequality at the marker in (1) gives only

```
log2 E[2^((delta-E delta)/8)]
 <= (ln2)*16384^2/512
 = (ln2)*2^19.
```

This is far above the target eleven bits in (3). The needed improvement
must use the actual connectivity, regional balance, and repair of local
matching defects. Dividing into independent windows by itself does not
control each window's exponential moment.

## Exact balance does not imply variance proportional to deficiency

Consider a separate family of bipartite matching graphs, indexed by an
integer `b>=1`. There are `4b` sources of demand one and `4b` candidate
inputs. Exactly `2b` inputs are chosen uniformly to have capacity two;
the others have capacity zero. This has exact total capacity balance.

Split sources and candidate inputs into corresponding halves of size
`2b`. Each source connects to every input in its own half and to none
in the other half. Let `X` count the positive-capacity inputs in the
first half. Then

```
X ~ Hypergeometric(4b,2b,2b),
matching deficiency Delta = 2*|X-b|.                 (4)
```

The formula follows by taking the minimum of demand and capacity in
each of the two complete components. Exact calculation gives

```
E Delta = b * binom(2b,b)^2 / binom(4b,2b),
E Delta^2 = 4b^2/(4b-1).                             (5)
```

For the first identity, telescope

```
(k-b)*binom(2b,k)^2
 = b*[binom(2b-1,k-1)^2-binom(2b-1,k)^2]
```

over `k>b`, then use symmetry. The second identity is four times the
hypergeometric variance. Stirling's formula yields

```
E Delta ~ sqrt(2b/pi),
Var Delta ~ (1-2/pi)*b.                              (6)
```

Thus the variance divided by the mean is unbounded, despite exact
balance and a small mean deficiency compared with the graph size.
This example is not the encoder's nine-neighbor graph. It rules out an
argument based only on balanced quotas, small mean deficiency, and the
generic matching interpretation. The actual random-quartet connectivity
must supply an additional property.

## Minimal Hall witnesses have an alternating representation

In the restricted matching graph, let `N(A)` contain the non-`H` input
neighbors of a source set `A`. Take an inclusion-minimal set satisfying
`|A|>2|N(A)|`. Removing any source leaves no deficient subset, so Hall's
theorem and integrality imply

```
|A|=2|N(A)|+1.                                      (7)
```

For any root in `A`, match the remaining sources. This matching fills
both slots of every input in `N(A)`. Explore alternately from the root:
from a source visit its input neighbors, and from an input visit its two
matched sources. The exploration reaches all of `A`. Otherwise its
reachable sources form a proper set with one more source than available
input slots, contradicting minimality.

Thus minimal defects have connected alternating witnesses, with two
matched source children per newly discovered input. This is a more
specific encoding than an arbitrary source subset. It is not a tail
bound: witnesses for several deficit units may overlap, and large
deficient components remain possible.

## One normalized Hall comparison still fails at `y=sqrt(2)`

The asymmetric normalization in
[the dual Hall counting note](QUARTET_SHORT_MATCHING_HALL_COUNT.md)
has an exact random-cover interpretation. Conditional on the setup,
select each nonempty quartet at updates two through 64 independently,
with probability

```
p_j=(1-y^(-2j))/(1-y^(-(8+j))),
```

where `j` is its number of non-`H` slots. Set `C_1=J_1` and let later
`C_t` consist of the non-`H` slots in selected quartets. If
`Z_cov=32768-|N^-(C)|`, that normalized Hall majorant equals
`E_coins[y^Z_cov]`. The flip argument bounds `y^delta` by this moment;
individual values of `Z_cov` need not dominate `delta`.

At `y=sqrt(2)`, write `g_j=1-p_j`, with `g_0=1`. Direct inequalities give

```
g_1>=15/32,  g_2>=7/32,  g_3>=1/10,  g_4>=3/64.
```

For independent half labels, the mean quartet absence probability is
therefore at least

```
g_bar >= (1+4*(15/32)+6*(7/32)+4/10+3/64)/16
      = 1483/5120 > 37/128.
```

For a source at times one through 62, its own label being `H` and both
future quartet coins being absent imply that the source is uncovered.
Fix the partitions first. These events involve nine distinct input
positions. Under the exact regional half quotas, each prescribed
nine-label pattern has probability at least `(511/512)^9` times its
independent-half probability: at each exposure the required label still
occupies at least `4096-8` of the original 8192 regional positions.
This factor is at least `503/512>31/32`. Consequently every such source
is uncovered with probability greater than

```
(31/32)*(1/2)*(37/128)^2
 =42439/1048576 > 1/25.
```

The last layer contributes an additional expected 256 uncovered `H`
sources. Ignoring the penultimate layer and applying Jensen gives

```
E Z_cov >= 62*512/25+256 = 38144/25,
log2 E[y^Z_cov] >= 19072/25 = 762.88.
```

At this marker, a sufficient window exponent for forty failure bits is
`(104860/2-8192-40)/64=690.59375`. Thus the asymmetric normalized Hall
comparison cannot establish that target at `y=sqrt(2)`. This is a lower
bound on the comparison, not on the true deficit moment. It does not
exclude markers `y=2` or `4`. The derivation received an independent
mathematical audit; no numerical job was used.

## Next proof objective

The unit-flow compression theorem in the fixed-row flow note proves
that the nine-neighbor matching has exactly the full flow's value: an
integral maximum flow can be decomposed into paths using at most one
temporal edge. Its Hall deficiency is therefore exactly `delta`. This
identity holds for every routing and partition realization at the
half-row capacities. A useful proof would charge that deficiency to
explicit deficient components or bounded alternating-path failures, then
control their joint exponential moment under exact regional quotas.

The moment must account for larger deficient components as well as
isolated sources. The existing interior isolation expectation does not
bound the total deficiency. Establishing (1), directly or through a
constructive loss, is a concrete stopping criterion for this branch.
No deficit-tail bound or whole-code certificate is asserted in this note.

## The exact capped-row-rank network already includes eligibility

For arbitrary nonnegative real `a,b`, the objective

```
min_S [R(V(S))-a*L(S)+b*D(S)]                         (8)
```

is an ordinary minimum cut. Give every state position a source arc of
capacity `a` and every forward coordinate transition an arc of capacity
`b`. Use infinite-capacity implications from a selected state to its
current footprint position and from a selected old state through its
fresh quartet to all four next footprint positions. Each footprint
position has a capacity-one arc to its outer-row auxiliary, whose arc
to the sink has capacity 128. There is no final death arc.

For a fixed selected state set, minimizing the footprint and row
auxiliaries contributes exactly
`sum_rows min(128,degree_V)=R(V)`. The source and temporal arcs contribute
`a*(N-L)+b*D`. Thus the minimum-cut value minus `a*N` equals (8).
An optimizing row auxiliary chooses the 128-cost option exactly when
its footprint degree exceeds 128; either choice has the same cost at
degree 128. Actual row eligibility is already present in this network.

This is an optimization identity, not a tail bound. In particular, the
full state set has `R=K`, `L=N`, and `D=0`, so it contributes `K-a*N`.
The structural obstruction also requires `R+D<=W-1`; that condition
cannot simply be dropped. A hard degree threshold, or a convex penalty
for exceeding a rank budget, is not generally another graph-cut term.
For example, the two-position set function `1[|S|>=2]` violates
submodularity on the two singleton sets. Auxiliary minimization of a
cut function preserves submodularity, so it cannot represent this
example. This does not rule out a constrained-cut method or a random
geometry theorem; it identifies the missing step in the plain network.

## Eligibility is almost automatic at the primary exact profile

Consider the exact profile

```
h=4096, R=819200, L=838861, f=65536, F=1081344.
```

Its numbers of covered `H` and outside slots are

```
C_H=F-R+128*h=786432=192*h,
C_out=R-128*h=294912=72*(8192-h).                    (9)
```

Fix the partitions, the footprint geometry, the row set `H`, and every
regional mask specifying which slots carry labels from `H`. Conditional
on this information, identities within `H` and within its complement
are uniform permutations, independently across regions. If region `j`
has `c_j` covered `H` slots, an individual `H` row has a Bernoulli
coverage indicator with probability `c_j/h`. Its total degree `Z` is
therefore a sum of independent regional indicators, with mean 192.
The same argument gives mean 72 for an outside row. Independence
between different rows is neither asserted nor needed.

The product moment bound and arithmetic-geometric mean inequality give

```
Pr[an H row has Z<=128]
 <= 3^128 * product_j(1-2*c_j/(3*h))
 <= (3/4)^128 < 2^-48,

Pr[an outside row has Z>=129]
 <= Pr[Z>=128]
 <= (41/32)^256 / 2^128 < 2^-32.                    (10)
```

The last strict comparisons follow from `3^8<2^13` and `41^8<2^43`.
A union over the two groups of 4096 rows bounds classification failure
by `2^-36+2^-20<2^-19`, uniformly over the fixed information above.

Consequently, imposing the exact eligible-row classification on any
nonnegative geometry weight that depends only on this fixed information
can reduce its weighted count by at most a factor `1/(1-2^-19)`.
Its binary logarithm is less than `2^-18`. This limitation applies to
such geometry counts, not to weights that retain within-group row
identities or RS amplitude equations. It rules out a large improvement
from row qualification alone at this exact profile; it is not a lower
bound on code failure probability.

## An exact obstruction to the `y=2` pair-cover comparison

The 129-state pair-cover operator in
[the Hall counting note](QUARTET_SHORT_MATCHING_HALL_COUNT.md) drops
the own-input absence condition. Its full exact-quota comparison is
already too large at `y=2`, independently of any Cauchy approximation
or choice of potential. The following single positive path proves this.

Keep exactly 96 hole quartets at all 4095 cover updates. The first
overlap, from the full-hole start, contributes `2^384`. Each of the
remaining 4094 overlaps is at least the single hypergeometric term

```
H_lo=2^300 * binom(384,300)*binom(128,84)/binom(512,384).
```

At this marker the exact hole probabilities are

```
g=(1,127/511,21/341,31/2047,1/273),   p_j=1-g_j,
G_j=binom(4,j)*g_j,                 P_j=binom(4,j)*p_j.
```

For a histogram `v`, write `Mult(v)=(sum v)!/product_j v_j!`.
One positive term in the 16-update regional quota coefficient is

```
h=(178,508,544,260,46),    k=(0,28,106,216,162),
c_lo=binom(128,96)^16 * Mult(h)*Mult(k)
     * product_j G_j^h_j P_j^k_j.
```

The two histograms contain respectively 1536 and 512 quartets, and
2560 and 1536 non-`H` labels. They therefore meet the regional quota
exactly. The analogous term for the first region, which has fifteen
cover updates and a free input-one label layer, is

```
h_end=(167,476,510,244,43),
k_end=(0,26,100,202,152),
e_lo=binom(128,96)^15 * Mult(h_end)*Mult(k_end)
     * product_j G_j^h_end_j P_j^k_end_j * binom(512,256).
```

Here the cover updates contain 2400+1440 non-`H` labels; input one adds
256. A zero exponent contributes one, including the unused `P_0=0`.

Let `C0=binom(8192,4096)`. The full pair-cover tail comparison is at
least

```
2^(512+384-104860) * H_lo^4094 * c_lo^255 * e_lo / C0^255.
```

The parent verified by exact rational arithmetic that

```
floor(log2 H_lo)=291,   floor(log2 c_lo)=4070,
floor(log2 e_lo)=4322,  ceil(log2 C0)=8186.
```

All histogram and quota identities passed in the same serial check.
An independent mathematical audit passed the path, coefficients,
endpoint powers, and quota normalizers.
These directed integers give a logarithmic lower bound

```
896-104860+4094*291+255*4070+4322-255*8186 = 42132.
```

Thus even exact regional coefficient extraction cannot make this
particular `y=2` pair-cover comparison certify the desired tail.
The lower bound concerns the comparison expression, not the actual
deficit moment, a realizable state support, or code failure. Retaining
own coverage or changing the proved comparison remains necessary.
