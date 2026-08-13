# Machine-readable contract for `g=8` support shards

This contract separates immutable run inputs, one proof tree per support, and
independent replay results.  Discovery fields may be present, but the verifier
must ignore them.

## Canonical encodings

Artifacts use UTF-8 JSON and a declared canonical-JSON profile.  SHA-256 is
computed over the canonical bytes.  Integers are decimal strings when they can
exceed the interoperable JSON integer range.  A rational is a reduced `p/q`
string with `q>0`.  A vector has exactly nine entries in class order `0,...,8`.

The immutable manifest has this logical shape:

```json
{
  "schema": "packet-group-g8-support-manifest-v1",
  "run_id": "...",
  "parameters": {
    "group_bits": 8,
    "N": "2097152",
    "K": "1048576",
    "M": "262144",
    "minimum_profile_weight": "21",
    "probability_target_log2": "-40/1"
  },
  "arithmetic": {
    "canonical_json": "RFC8785",
    "integer_split_rule": "primitive-affine-gap-v1",
    "outward_evaluator": "...",
    "outward_evaluator_sha256": "..."
  },
  "census": {
    "rule": "exact-support-positive-compositions-weight-cut-v1",
    "feasible_masks": [{"mask": "0x002", "count": "1"}],
    "infeasible_masks": ["0x000", "0x001"],
    "feasible_mask_count": 510,
    "feasible_profile_count": "553169839211945865258921061892182603726"
  },
  "witness_sources": [
    {"source_id": "atlas-0001", "sha256": "...", "schema": "..."}
  ],
  "fixed_sources": [{"source_id": "spectrum-01", "sha256": "..."}]
}
```

The manifest is frozen before work begins.  A shard records the SHA-256 of
the exact manifest bytes.  It refers to sources by `source_id` and digest, not
by an environment-dependent path.

## Support shard

One file owns one feasible mask:

```json
{
  "schema": "packet-group-g8-support-shard-v1",
  "manifest_sha256": "...",
  "support_mask": "0x1ff",
  "active_classes": [0,1,2,3,4,5,6,7,8],
  "root_count": "...",
  "root_node": "r",
  "nodes": [],
  "aggregation_requested": "expanded|collapsed|minimum",
  "state": "COMPLETE|INCOMPLETE"
}
```

The mask, active classes, and root count must equal the manifest census row.
The root domain is reconstructed, not trusted: active coordinates are at
least one, inactive coordinates are zero, their sum is `M`, and
`sum(j*a_j)>=21`.  Exact support selects the shard before any tree decision.

Each node has exactly one of these states.

- `EMPTY` contains `node_id` and an emptiness method.  The verifier proves
  integer emptiness; it does not trust a stored zero count.
- `UNRESOLVED` contains only diagnostic references.  Its presence makes the
  shard incomplete.
- `SPLIT` contains `node_id`, `left`, `right`, and
  `split={"coefficients":[...],"threshold":"t"}`.  The left child owns
  `A(a)<=t`; the right child owns `A(a)>=t+1`.
- `CERTIFIED_LEAF` contains `node_id`, one fixed witness selector, and optional
  count data.  The selector is either one witness or one fixed rational
  mixture.

For every split, the nine coefficients and threshold are integers.  The
coefficient vector is nonzero and primitive.  The verifier rejects rational
or floating cuts.  To canonicalize an integer form, divide its coefficients
by their gcd and replace the threshold by its floor after the same division.
If the first nonzero coefficient is negative, negate the form, replace `t` by
`-t-1`, and exchange the children.  These operations preserve integer
ownership.  The verifier reconstructs both child domains from the root and path,
so serialized vertices and incidence claims are diagnostic only.  Every node
must be reachable exactly once from `r`; cycles, duplicate parents, missing
children, and unreachable nodes are errors.

A witness reference has

```json
{"source_id":"atlas-0001", "source_sha256":"...", "row":17,
 "row_sha256":"..."}
```

The row digest binds the canonical witness row and prevents row-reordering
ambiguity.  A mixture stores two or more such references with nonnegative
reduced rational weights whose exact sum is one.  The verifier checks that
every zero fugacity occurs outside the shard's support.  It reconstructs each
leaf's exact rational vertices and outward-evaluates the fixed selector at
all vertices.  Stored vertex hashes, predicted scores, and claimed upper
bounds are not proof inputs.

## Counts and aggregation

Each count record declares one of:

```text
{"kind":"exact", "method":"coordinate-box-ie-v1", "value":"..."}
{"kind":"exact", "method":"laminar-convolution-v1", "value":"..."}
{"kind":"upper", "method":"verified-min-v1", "value":"..."}
```

The verifier recomputes every value.  It subtracts the manifest's exact table
of weight-at-most-20 exclusions.  Whenever both children claim exact counts,
their counts must sum to the exact parent count.  An upper count must not be
used as an exact conservation value.

Expanded aggregation is permitted only when every contributing terminal has
an exact owned count.  It evaluates

```text
sum_leaf n_leaf * 2^U_leaf.
```

Collapsed aggregation needs an exact ancestor count and complete terminal
coverage.  It evaluates

```text
n_ancestor * 2^(max_leaf U_leaf).
```

If both modes are available, the verifier may retain their smaller outward
upper bound.  A witness mixture and a tree branch are proof devices, not new
events, and receive no additional union factor.

## Replay receipt and gates

The independent verifier emits a receipt rather than modifying a shard:

```json
{
  "schema": "packet-group-g8-support-replay-v1",
  "manifest_sha256": "...",
  "shard_sha256": "...",
  "support_mask": "0x1ff",
  "status": "VERIFIED_COMPLETE|REJECTED",
  "root_count": "...",
  "node_counts": {"empty":0,"leaf":1,"split":0,"unresolved":0},
  "aggregation_used": "collapsed",
  "contribution_log2_interval": ["lo","hi"]
}
```

The verifier applies these gates in order.

1. **Manifest gate.** Validate parameters, source digests, arithmetic policy,
   all 510 feasible masks, the rejection of `{0}`, each exact support count,
   and their sum to the stated feasible profile count.
2. **Ownership gate.** Require one shard for the selected mask.  Replay every
   integer split using `<=t` and `>=t+1`; reject cross-shard node references.
3. **Tree gate.** Reconstruct all domains and exact vertices.  Require one
   reachable root, valid node arities, and no unresolved node.
4. **Witness gate.** Resolve and hash every witness row, check support
   eligibility and mixture weights, then replay each outward vertex bound.
5. **Count gate.** Recompute counts and conservation identities.  Check that
   the requested aggregation mode has the required exact counts.
6. **Shard-completion gate.** Emit `VERIFIED_COMPLETE` only after all preceding
   gates pass.  A shard's own `COMPLETE` claim is insufficient.
7. **Integration gate.** Require exactly one successful receipt for each of
   the 510 feasible masks, no duplicate mask, one manifest digest, root counts
   summing to the exact census, and an outward log-sum-exp upper endpoint at
   most `-40`.

Reconnaissance may stop with incomplete shards and diagnostic binary64
scores.  Certification integration consumes only immutable shards and their
independent `VERIFIED_COMPLETE` receipts.
