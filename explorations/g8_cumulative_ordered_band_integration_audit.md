# Integration audit for cumulative boxes and ordered-band slabs

Both prototypes define sound proof cells for the full-support `g=8` root.
Neither needs relaxed proof geometry.  Both can use the current exact split
rule

```text
left: A(a)<=t,             right: A(a)>=t+1.
```

The current artifacts are not certificate shards.  The ordered-band result
also shows that its present witness bank does not certify any leaf.

## The theorem interface

Fix one full-support integer profile set

```text
Lambda={a in Z^9 : a_j>=1, sum_j a_j=M}.
```

The physical-weight constraint is redundant because `sum_j j>=36`.  A leaf
needs two objects.

- Its owned set is the subset of `Lambda` selected by the root-to-leaf
  integer comparisons.
- Its proof polytope is the root simplex intersected with the corresponding
  closed real halfspaces.

The exact-gap rule assigns every integer profile to one child.  It can leave
an open real gap between the child polytopes.  That gap is harmless because
the theorem sums only over integer profiles.

For each leaf, the verifier fixes one support-eligible witness or one fixed
rational mixture.  The witness exponent is convex on the proof polytope.
A nonnegative fixed mixture preserves convexity.  Therefore, outward bounds
at every exact vertex bound every owned integer profile in that leaf.

This argument requires no stored vertex list.  The verifier must reconstruct
the vertex hull from the root and path constraints.

## Cumulative boxes

For active classes `s_0<...<s_8`, put

```text
x_i=a_(s_i)-1,              T=M-9,
y_k=sum_(i=0)^k x_i         for 0<=k<8.
```

The integer map `a -> y` is unimodular, and its image is

```text
0<=y_0<=...<=y_7<=T.
```

Balanced integer bins partition `{0,...,T}`.  A cell fixes one bin for every
`y_k`.  Its bin-index tuple is nondecreasing.  Thus every profile owns one
cell, and no profile owns two cells.

Equal consecutive bin indices form runs.  A run of length `r` in an integer
interval `[L,U]` has the closed proof polytope

```text
L<=z_1<=...<=z_r<=U.
```

Its vertices are `(L,...,L,U,...,U)`, with one transition position.  Distinct
runs use disjoint ordered intervals.  Hence the complete cell is a product
of run polytopes.  Its stored product vertex list is exact, but the verifier
must derive the same list independently.

If a run interval has `w=U-L+1` integers, that run owns

```text
C(w+r-1,r)
```

integer sequences.  The cell count is the product over runs.  For full
support, no exclusion correction is needed.  The counts over all
nondecreasing bin tuples sum to `C(M-1,8)`.

Each bin boundary uses a prefix split

```text
sum_(i=0)^k a_(s_i) <= U+k+1.
```

Because `k<8`, the full-support pivot class is absent.  This split already
uses the manifest's support-relative canonical form.

The prototype's level-two artifact has 165 cells, at most 81 vertices per
cell, and 6,435 vertex incidences.  Its exact cell counts sum to the frozen
full-support census.

## Ordered-band slabs

Fix a proper class band `B`, with `q=|B|` and `r=9-q`.  A leaf fixes

```text
L <= m_B(a) <= U,           m_B(a)=sum_(j in B) a_j.
```

Integer thresholds partition the feasible masses `q,...,M-r` into disjoint
intervals.  Therefore, every full-support profile owns exactly one slab.

At an internal boundary mass, a vertex has one nonunit band coordinate and
one nonunit complement coordinate.  Endpoint slabs also retain the relevant
simplex corners.  These points are the complete vertex set of the real slab,
not merely integer samples.  A slab has at most

```text
9+2q(9-q) <= 49
```

vertices.  The current diagnostic observed at most 40.

The exact slab count is

```text
sum_(m=L)^U C(m-1,q-1) C(M-m-1,r-1).
```

Adjacent leaf intervals are disjoint, and their counts sum to the root
census.  Canonicalizing a band that contains the pivot exchanges its two
serialized children.  The producer already records this exchange.

The ordered-band artifact has seven leaves and exact count conservation.
However, all seven leaves are `UNRESOLVED`.  Its binary64 discovery pass
found zero candidate passes.  Exact geometry and counts do not replace the
missing outward witness inequalities.

## Compatibility with the current verifier

The geometry is compatible with `support-relative-primitive-affine-gap-v2`.
No new split semantics or relaxed proof polytope is necessary.  The current
generic vertex enumerator can reconstruct both families from canonical path
constraints.

The cumulative prototype is a cell list rather than a binary tree.  A shard
producer must serialize a reachable tree whose leaves are precisely the
nonempty bin tuples.  It must avoid impossible tuples because the verifier
rejects an empty or redundant split child.

The ordered-band artifact stores a nested diagnostic tree.  The shard schema
requires a flat node array with stable `node_id`, `left`, and `right` fields.
Diagnostic `cell` and stored vertex fields cannot be proof inputs.

Neither prototype is a complete shard.  A complete artifact needs:

1. schema `packet-group-g8-support-shard-v1`;
2. the current manifest digest and full-support root census;
3. one flat, reachable canonical split tree;
4. no `UNRESOLVED` node;
5. one manifested fixed selector at each `CERTIFIED_LEAF`;
6. state `COMPLETE` and a declared aggregation mode.

Only `atlas-3934dae` is a witness source in the current manifest.
Supplementary rows may guide discovery, but a leaf cannot reference one.
Using a supplementary row as a selector requires a new manifested witness
source and therefore a new manifest digest.  Unmanifested diagnostic inputs
should be omitted from the immutable proof shard.

## Counts and aggregation

Collapsed aggregation already suffices for either geometry.  It uses

```text
C(M-1,8) * 2^(max_leaf U_leaf).
```

The verifier needs only the manifest root count, complete tree coverage, and
the outward leaf bounds.  A producer may omit nonroot count records.  This
is the smallest schema-compatible route to a complete shard.

Expanded aggregation needs two verifier changes.

1. Add independently recomputed exact methods, for example
   `cumulative-run-product-v1` and
   `two-band-positive-composition-scan-v1`.
2. Feed those verified terminal counts into aggregation.

The second change is essential.  The current verifier constructs expanded
terms only from bounded profile enumeration.  Full support exceeds that
limit, so adding count parsers alone would still leave expanded aggregation
unavailable.

For a cumulative count record, bind the support order, complete bin
partition, bin-index tuple, `T`, and exact product.  For an ordered-band
record, bind the canonical band, `L`, `U`, `M`, and the exact binomial sum.
The verifier must reconstruct these parameters from the path.  It must not
trust a producer's redundant `cell` description.

Exact counts must be positive at every serialized child and conserve at each
split.  Upper counts may support collapsed aggregation but must not appear in
an exact conservation identity.

## Recommendation

Use cumulative boxes when the available fixed mixtures certify most of the
165 level-two cells.  They give materially finer localization with a modest
vertex total.  Use the ordered-band tree only if a later witness bank closes
its large transverse slabs; the current diagnostic provides no such evidence.

First build a collapsed cumulative-box shard with manifested selectors.
That path needs no count-method extension.  Add expanded aggregation only
after an independent exact-count replay and aggregation data path exists.
