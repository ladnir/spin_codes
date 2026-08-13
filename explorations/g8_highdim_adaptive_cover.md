# High-dimensional `g=8` adaptive coverage producer

The high-dimensional lane covers exact supports of dimensions 8, 7, and 6.
These 46 supports contain nearly all feasible profiles.  Every such support
also has minimum physical weight at least 21.  Therefore, its root count is
the positive-composition count without a low-weight correction.

The producer is
`scripts/produce_packet_group_g8_highdim_adaptive_shards.py`.  It binds to the
frozen manifest and atlas by their repository-stable SHA-256 digests.  It also
recomputes every canonical atlas-row digest from the manifest before work.

## Ownership geometry

For one exact support, the root has integer coordinate bounds

```text
a_j >= 1 for active j,       a_j = 0 otherwise,       sum_j a_j = M.
```

The producer uses coordinate slabs only.  A split on class `j` has the exact
integer ownership rule

```text
left: a_j <= t,              right: a_j >= t+1.
```

The threshold approximately balances the exact child counts.  The producer
computes each coordinate-box count by inclusion--exclusion.  It records the
lower bounds, upper bounds, total mass, active classes, and zero weight-cut
exclusion with every count.  Both child counts must sum to the parent count
before serialization.  The serialized verifier-facing record uses
`kind=upper` because the verifier permits `coordinate-box-ie-v1` as a
recomputed upper count.  A diagnostic flag records that the count is exact
for the producer's coordinate-only path.

The vertices of a box intersected with the sum hyperplane are integral.  The
producer enumerates them by fixing all but one active coordinate at a bound.
This costs at most `s*2^(s-1)` candidates for support size `s`.  Thus a
full-support node has at most 2,304 vertex candidates before deduplication.

## Discovery and certification boundary

All 510 frozen affine atlas witnesses have positive fugacities.  The producer
may therefore compare them on every high-dimensional support.  It chooses the
witness with the smallest binary64 maximum across a node's exact vertices.
The resulting score is diagnostic because its normalization and affine data
use binary64 arithmetic.

No producer terminal is a `CERTIFIED_LEAF`.  A diagnostic pass becomes an
`UNRESOLVED` node with reason `AWAITING_INDEPENDENT_OUTWARD_REPLAY`.  A node
that reaches a depth or node budget also remains `UNRESOLVED`.  Consequently,
every output shard declares `state=INCOMPLETE`.  A later hardener must resolve
the frozen witness reference, reconstruct the path domain, replay all exact
vertices with outward arithmetic, and replace a successful terminal with a
canonical certified leaf.

Arbitrary affine leaves do not claim exact owned counts.  This producer uses
coordinate boxes, so it records exact leaf counts.  Nevertheless, each shard
requests collapsed aggregation.  This remains compatible with later affine
refinement, where only the exact ancestor count may be available.

## Scaling and first batch

One node performs at most about `510 * 2304` affine scores and at most
`9 * log2(M) * 2^9` small exact-integer count terms for split selection.  A
63-node shard therefore has a conservative ceiling near 74 million affine
scores.  Actual vertex counts shrink after splitting.  The output size is
linear in the node budget.  A split adds exactly two nodes, so the node budget
must be odd.

The first root-run batch should measure the full-support shard alone with a
31-node discovery budget:

```powershell
python scripts\produce_packet_group_g8_highdim_adaptive_shards.py `
  --output-dir out\g8_highdim_batch000 `
  --dimensions 8,7,6 `
  --shard-start 0 `
  --max-shards 1 `
  --max-nodes-per-shard 31 `
  --max-depth 10 `
  --hardening-reserve-bits 16
```

The command is a discovery batch, not a benchmark.  Run it alone.  The
conservative ceiling is about 37 million affine scores, and actual vertex
counts should be smaller.  Inspect `batch_index.json` before scheduling the
dimension-7 roots.  Rank residual nodes by their stored collapsed candidate
contribution, then decide whether the next stage should increase the
coordinate-slab budget or add an affine cut policy.
