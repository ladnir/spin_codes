# Exact cover strategy for `g=8` support dimensions zero through five

This lane owns the 464 feasible supports of dimensions zero through five.
The frozen witness source is `out/g8_support_seed_atlas.json` at SHA-256
`c48644fa62ad74619a4db30b07950b14aaf1419d040760f7d659f43288fca45e`.
Its 510 rows are binary64 discovery witnesses.  They are not certificates.
Every plan and shard binds `G8_SUPPORT_MANIFEST.json` at SHA-256
`dafff5c978d51515355960293832223a8a4736ca4b2c0405496c89d395ac3890`.

## Direct root certificates

Fix an exact support `S` of size `s`, and put `M=262144`.  Its proof domain is

```text
P_S={a in R^9 : a_j>=1 for j in S, a_j=0 for j not in S,
                    sum_j a_j=M, sum_j j*a_j>=21}.
```

A frozen combined witness `w` has the profile exponent

```text
F_w(a)=C_w-q_w*a-log2(M!/prod_j a_j! * prod_j binom(8,j)^a_j).
```

The map `F_w` is convex on `P_S`.  Indeed, `log Gamma(x+1)` is convex, and
all remaining terms are affine.  Therefore,

```text
sup_{a in P_S} F_w(a)=max_{v in vertices(P_S)} F_w(v).
```

This identity makes a one-leaf root the cheapest sound certificate.  The
independent verifier reconstructs `P_S`, outward-evaluates one fixed selector
at every exact vertex, and retains their largest upper endpoint `U_S`.

The root vertices have a closed form.  Write `a_j=1+x_j` on `S`, and put
`T=M-s`.  Without an active class zero, the vertices allocate all `T` units
to one active class.  The physical-weight cut is then redundant.

Suppose class zero is active, and put

```text
h=21-sum_{j in S, j>0} j.
```

If `h<=0`, the same pure residual vertices apply.  If `h>0`, the class-zero
vertex is removed.  The remaining vertices are:

```text
x_p=T                         for each p in S minus {0};
x_p=h/p, x_0=T-h/p           for each p in S minus {0}.
```

All values are exact rationals.  A dimension-five root has at most ten
vertices, so direct outward replay is small.

The exact root count is

```text
n_S=C(M-1,s-1)                         if 0 is not in S;
n_S=C(M-1,s-1)-E(S minus {0},h)        otherwise,
```

where `E` is the finite weighted-partition count defined in
`g8_domain_decomposition.md`.  The final root contribution is

```text
n_S * 2^U_S.
```

The low-dimensional lane is complete when all 464 roots have verified
one-leaf shards, or verified trees for the residual roots, and their receipts
pass the integration gate.  Binary64 planning results do not change this
condition.

## Selector order and fallback cells

For each root, discovery tries selectors in this order.

1. Try the atlas row tuned on `S`.
2. Try every support-eligible single atlas row.
3. Try a fixed rational mixture of a small nondominated witness set.
4. Split only a root that still fails its assigned ledger budget.

All frozen atlas fugacities are positive.  Thus every atlas row is eligible
on every support.  The final verifier must still perform the eligibility
check from the manifest source.

A mixture is useful before geometry because it preserves one root leaf.  The
planner solves the finite vertex game by dual column generation.  Every round
prices all eligible atlas rows, so a per-vertex candidate cap cannot hide an
LP-useful witness.  The final binary64 solution is rounded by stable largest
remainder to a declared common denominator.  Zero components are removed.
The emitted weights are nonnegative reduced rationals with exact sum one.
The planner then re-evaluates the rational mixture and records any loss from
rounding.  The verifier uses only these rational weights and outward vertex
evaluations.  A one-component result is encoded as a witness, because the
canonical mixture selector requires at least two components.

The first geometric fallback uses an active coordinate `j` and integer `t`:

```text
left: a_j<=t,                 right: a_j>=t+1.
```

The split is already primitive and follows `primitive-affine-gap-v1`.
Choose `j` from the largest difference between the active selectors' affine
charges.  Choose `t` near their diagnostic equality point, clamp it to the
node's integer range, and reject an empty child.  This rule uses a floating
score only to propose a split; exact integer ownership does not depend on the
score.

If coordinate splits stagnate, use a laminar band mass

```text
A_B(a)=sum_{j in B} a_j
```

for `B=S intersect {0,...,r}` or `B=S intersect {r+1,...,8}`.  These splits
retain exact `laminar-convolution-v1` counts.  Arbitrary witness-difference
cuts are a last resort because they complicate exact counts and canonical
integer scaling.

Every terminal polytope is reconstructed from the root and its integer path.
A fixed witness or fixed rational mixture must pass every exact terminal
vertex.  A centroid, sampled profile, or binary64 incidence calculation is
never a coverage argument.

## Artifact alignment

A successful direct root maps to `packet-group-g8-support-shard-v1` as follows.

```json
{
  "support_mask": "0x...",
  "active_classes": [],
  "root_count": "...",
  "root_node": "r",
  "nodes": [{
    "node_id": "r",
    "state": "CERTIFIED_LEAF",
    "selector": {
      "kind": "witness",
      "witness": {
        "source_id": "atlas-3934dae",
        "source_sha256": "c486...a45e",
        "row": 0,
        "row_sha256": "..."
      }
    },
    "count": {"kind":"exact","method":"root-census-v1","value":"..."}
  }],
  "aggregation_requested": "collapsed",
  "state": "COMPLETE"
}
```

The planning script copies each row digest from the frozen manifest and emits
`proof_state: REQUIRES_OUTWARD_REPLAY`.  The manifest uses
`sorted-compact-json-v1`.  The planner cannot produce `COMPLETE` shards.

## First bounded run

Run dimensions zero through two first.  This checks 128 supports and uses
only affine binary64 evaluation of the frozen rows.  It does not invoke the
inner kernel or optimize a witness.

```powershell
python scripts/plan_packet_group_g8_lowdim_root_cover.py `
  --manifest G8_SUPPORT_MANIFEST.json `
  --min-dimension 0 --max-dimension 2 --candidate-scope all `
  --selector-mode minimax-mixture --mixture-seed-per-vertex 1 `
  --mixture-denominator 1073741824 `
  --output out/g8_lowdim_root_mixture_d0_2.json
```

If nearly all roots pass the uniform target, run dimensions three through
five with the same planner.  Outward-harden only the selected rows.  If many
roots fail, inspect the largest exact contributions before introducing any
split tree.
