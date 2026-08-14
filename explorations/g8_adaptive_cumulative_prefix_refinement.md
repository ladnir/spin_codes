# Adaptive cumulative-prefix refinement after the `h=2` diagnostic

The level-two cumulative diagnostic has 165 cells and 6,435 vertex
incidences.  No cell passes its binary64 target.  The worst diagnostic
contribution exceeds the target by about 1.45 million bits.  A broad level-three
grid is therefore premature.  The next wave should test whether local
prefix refinement changes the witness landscape materially.

## Refinement invariant

For the full-support profile, define

```text
x_i=a_i-1,                  T=M-9,
y_k=sum_(i=0)^k x_i         for 0<=k<8.
```

Every node stores integer intervals

```text
L_k <= y_k <= U_k,
```

together with `0<=y_0<=...<=y_7<=T`.  A split chooses one coordinate `k`
and one integer `t` with `L_k<=t<U_k`:

```text
left:  y_k<=t,
right: y_k>=t+1.
```

In profile coordinates, the canonical left split is

```text
sum_(i=0)^k a_i <= t+k+1.
```

The left child replaces `U_k` by `t`.  The right child replaces `L_k` by
`t+1`.  Thus every descendant remains a cumulative-interval cell.  No
arbitrary affine direction enters the tree.

The producer must recompute monotone closure after each split:

```text
L_k <- max_(i<=k) L_i,
U_k <- min_(i>=k) U_i.
```

Reject a child if some `L_k>U_k` or its exact integer count is zero.  This
test prevents a formally nonconstant split from creating an empty child.

Do not install a selected threshold as a global grid boundary.  A local
threshold belongs only to the chosen node and its descendants.  Global
installation can multiply unrelated cells without improving their witness
fit.

## Selecting a prefix coordinate

Fix one failed cell and its current vertex set `V`.  Use the complete frozen
discovery bank.  For every `v in V`, record a leading witness `w(v)`, with
ties resolved by the witness reference.

For each coordinate `k`, group vertices by the other seven cumulative
coordinates.  Within one group, sort the vertices by `y_k`.  An adjacent pair
is a transition if its leading witnesses differ.  Each transition defines a
candidate one-dimensional edge segment.

For a transition from `u` to `v`, let the two leading affine witness scores
be `F_p` and `F_q`.  Their common normalization cancels.  Restrict
`F_p-F_q` to

```text
z(lambda)=(1-lambda)u+lambda v.
```

If the binary64 affine difference changes sign, compute its crossing
coordinate `c`.  Propose

```text
t in {floor(c), ceil(c)-1} intersect [y_k(u),y_k(v)-1].
```

If numerical cancellation makes the crossing unusable, propose the integer
midpoint.  The crossing is discovery data only.  The serialized integer
prefix split determines ownership.

Evaluate every proposed `(k,t)` by constructing both exact child hulls and
rerunning the cell minimax discovery.  Rank candidates lexicographically by:

1. the smaller maximum child contribution `log2(n_child)+U_child`;
2. the smaller maximum child witness bound `U_child`;
3. the larger minimum exact child count;
4. the smaller total child vertex count;
5. smaller `(k,t)`.

Require a positive diagnostic improvement over the unsplit parent.  Otherwise,
leave the cell unresolved.  This rule prevents subdivision that changes only
the artifact size.

## Exact ownership and counts

The root-to-leaf integer comparisons form a disjoint partition.  Every
integer profile takes exactly one branch at each node.  Local refinements in
different parents cannot overlap in ownership.

The h2 run-product formula applies only while every coordinate uses one
shared base-bin interval.  Independent local splits can give different
intervals within an old run.  Use the general monotone-chain dynamic program
instead.

Let `D_k(z)` count prefixes ending at `y_k=z`.  Initialize

```text
D_0(z)=1 for L_0<=z<=U_0.
```

For `k>0`, compute

```text
D_k(z)=sum_(u=L_(k-1))^min(z,U_(k-1)) D_(k-1)(u)
        for L_k<=z<=U_k.
```

One prefix-sum pass computes each layer.  The exact node count is

```text
n=sum_z D_7(z).
```

The cost is `O(8T)` big-integer additions and `O(T)` memory per node.  The
verifier must recompute this count before expanded aggregation.  Collapsed
aggregation may continue to use the exact root census.

Every accepted split has positive child counts, and the recurrence must
satisfy

```text
n_left+n_right=n_parent.
```

## Vertex guarantees

A cumulative-interval cell is an order polytope on a chain with individual
lower and upper bounds.  Refinement changes endpoints but introduces no new
facet direction.  Its vertices are exact integers because the chain
constraint matrix has the consecutive-ones property and is totally
unimodular.

The generic exact verifier can enumerate the hull from the path constraints.
A depth-independent conservative bound in dimension eight is

```text
V_0=1,
V_d=sum_(r=1)^d 2r V_(d-r),
V_8=21728.
```

This bound counts anchored constant blocks and applies after arbitrary local
prefix refinements.

The first wave has a sharper bound.  An h2 parent is a product of run
polytopes with at most 81 vertices.  One prefix cut divides one run interval
at one coordinate.  Direct anchored-block counting gives at most 320 vertices
per child.  Hence eight first-wave splits create at most 5,120 child vertex
incidences.  The earlier value 2,560 was an arithmetic error.  The producer
must record the actual total and reject any child above the per-child cap.

## Bounded first wave

Rank the 165 failed cells by the diagnostic value

```text
log2(exact_count)+mixture_upper.
```

The present top cells begin with indices

```text
4, 162, 10, 161, 20, 35, 160, 163.
```

Process exactly these eight cells.  For each cell:

1. generate transition candidates for all eight prefixes;
2. retain at most two thresholds per transition;
3. deduplicate exact `(k,t)` pairs;
4. accept at most one split;
5. recompute both child counts and exact hulls;
6. rerun rational-mixture discovery on each child.

The wave adds at most eight split nodes and sixteen children.  It performs
no outward hardening and makes no completion claim.  Keep all other h2 cells
unchanged.

Stop after this wave.  Continue only if the maximum child contribution drops
materially on several parents.  A useful planning gate is at least 16 bits
of improvement on four of the eight parents.  This threshold is diagnostic,
not a theorem condition.  Given the million-bit h2 gaps, a weak improvement
means that cumulative refinement is not the controlling remedy.

## Artifact requirements

The adaptive diagnostic must bind the h2 artifact digest and current
manifest digest.  Each split record contains `k`, `t`, the canonical profile
coefficients, exact child counts, and exact vertex digests.  Store transition
witnesses and crossing values under `diagnostic` only.

A later certificate shard discards all stored vertex claims.  The independent
verifier reconstructs the path hull, verifies the integer count when used,
resolves one manifested fixed selector per leaf, and performs every outward
vertex inequality.

## First-wave result

The bounded wave accepted one split for each prescribed parent.  The minimum
improvement was `22533.496307770256` bits.  The maximum improvement was
`66713.07030917914` bits.  Thus all eight parents exceed the 16-bit planning
threshold, and the stop gate passes with count eight.

Replacing the eight parent contributions by their sixteen children reduces
the global log2 union diagnostic from `1454125.2595775214` to
`1415166.1517534186`.  The new worst leaf is the right child of parent 162.
Its contribution is `1415166.1517534186`.

The sixteen children have 540 vertex incidences in total.  The largest child
has 48 vertices.  Both values are below the corrected conservative caps.
Every child selector uses at least one base witness and one supplementary
witness.

The mean parent improvement is about `46583.17` bits.  A purely linear
extrapolation needs about 31 comparable waves to move the current aggregate
below `-40`.  Using the observed minimum improvement gives 63 waves.  These
figures are planning diagnostics only.  Later refinements can change the
number of active children and the attainable improvement.

The wave artifact is
`out/g8_full_support_adaptive_cumulative_wave1.json`.  Its SHA256 digest is
`a73ea19d44a59662983c3c7c42ad6069a5800fe02fb3480c9a796038fd1f7e9a`.
