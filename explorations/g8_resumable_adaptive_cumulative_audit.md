# Audit of the adaptive cumulative discovery engine

The engine preserves exact integer ownership and the cumulative-cell family.
Its current scoring and persistence rules are not yet suitable for a
resumable discovery ledger.  This note separates those defects from later
certificate conversion.

## Required frontier state

At wave `w`, let `F_w` be the current leaves.  Each leaf `ell` has an exact
path, exact count `n_ell`, exact vertex hull `V_ell`, one fixed rational
discovery selector, and binary64 diagnostic bound `U_ell`.

The producer must derive `n_ell` and `V_ell` from the path.  Stored parent
vertices and counts are inputs to check, not facts to inherit.  Every split
must satisfy

```text
n_left>0, n_right>0, n_left+n_right=n_parent.
```

The root-to-leaf integer comparisons then prove disjoint ownership.  The
later certificate verifier reconstructs the paths and ignores producer vertex
claims.

## Ranking and contribution accounting

For one leaf, define its expanded diagnostic term

```text
E_ell=log2(n_ell)+U_ell.
```

Rank leaves by decreasing `E_ell`, then by canonical path.  Ranking by
`U_ell` alone can prioritize a negligible leaf.

For a candidate split, the two valid aggregate diagnostics are

```text
expanded:
E_split=log2(n_left*2^U_left+n_right*2^U_right),

collapsed:
C_split=log2(n_parent)+max(U_left,U_right).
```

The current engine instead uses the larger individual child term.  That
value omits one child, can be one bit below `E_split`, and is not `C_split`.
It can mis-rank candidates and overstate improvement.

Record both valid aggregates.  Use expanded improvement as the primary key
because exact child counts are available.  Use collapsed improvement as the
secondary key because the current proof verifier supports that mode.

The parent and child selectors must use the same bank, denominator,
tolerances, and tie rules.  The engine reads the parent's old contribution
while allowing another child denominator.  It must recompute the parent under
the current configuration or reject the mismatch.

## Transition discovery

Binary64 affine witness differences may propose a threshold.  The emitted
integer prefix split alone determines ownership.

The current routine examines adjacent vertices that agree in seven cumulative
coordinates.  It finds axis-aligned hull transitions, but not edges whose
endpoints change several coordinates.  Therefore, an emitted proposal is
sound, but an empty proposal set does not establish that refinement is useless.

For each coordinate without a transition, also evaluate one exact count
median.  Choose the smallest threshold whose left child owns at least half
the parent profiles.  This deterministic fallback prevents a false no-go
caused by incomplete edge sampling.

Reject nonfinite scores and crossings.  Deduplicate `(k,t)` pairs before
evaluation.  Rank candidates by

```text
(expanded aggregate, collapsed aggregate, -minimum child count,
 total child vertices, k, t).
```

Serialize binary64 ranking values by hexadecimal form or raw 64-bit encoding.

## Fixed mixtures

Every child receives one fixed rational selector.  A later proof cannot
choose a mixture component independently at each vertex.

The bank includes 32 supplementary rows that are not manifest sources.  They
may guide discovery but cannot enter a shard.  Before proof conversion,
either regenerate the manifest or recompute terminal selectors using only
manifested rows.

The independent verifier then outward-evaluates each fixed rational selector
at every reconstructed vertex.  Binary64 LP optimality and crossings have no
proof role.

## Deterministic checkpoints

The present script writes one final wave artifact.  It has no resumable
checkpoint state.

A checkpoint must bind the manifest, parent artifact, witness sources, code
closure, numerical configuration, exact frontier, completed-result digests,
and ordered pending queue.  Define a work item by `(wave, leaf_path)`.

Sort vertices, witnesses, transition groups, candidates, and mixture
components canonically.  Process the smallest pending key.  Write canonical
LF JSON atomically after each complete item.

Do not hash wall time, host paths, process identifiers, or timestamps into
the checkpoint.  Put telemetry in a separate receipt.  On resume, verify all
bindings and completed-result digests.  Never merge checkpoints from different
numerical configurations.

## Safe stops and finite-wave decision

Resource stops are safe when every unfinished leaf remains `UNRESOLVED`.
Stop only at work-item boundaries.  Caps may limit accepted splits, frontier
size, depth, vertex incidences, candidate evaluations, or elapsed time.

A numerical failure also leaves its leaf unresolved.  It must never create an
empty or certified node.

The 16-bit gate is too weak against diagnostic gaps near 1.45 million bits.
Use two validation waves before broad work.  Each wave selects eight frontier
leaves by `E_ell` and accepts at most one split per leaf.  Declare `GO` only if:

1. at least six leaves split in each wave;
2. the exact expanded aggregate over selected leaves decreases in each wave;
3. the lower-quartile positive improvement `Delta_25` satisfies

   ```text
   ceil(G_max/Delta_25)<=8,
   ```

   where `G_max=max_ell(E_ell+40)` after wave two;
4. expanding each dominant branch to that projected depth yields at most
   2,048 frontier leaves;
5. every accepted split satisfies count, vertex, and replay caps.

If `Delta_25<=0`, declare `NO_GO`.  Record `INCONCLUSIVE` instead when missing
transition proposals caused failure and count-median fallbacks were absent.
A no-go redirects effort; it does not prove that cumulative certification is
impossible.

The bounded first wave has eight splits, sixteen children, and a correct
worst-case cap of `16*320=5120` child vertex incidences.

## Proof conversion

A checkpoint becomes a shard candidate only with no unresolved leaf.
Conversion must flatten the tree, discard producer vertex claims, use only
manifested selectors, preserve every split, and request collapsed aggregation
unless verified exact counts feed aggregation.

Independent outward replay at every reconstructed leaf vertex remains the
certificate gate.  Until then, the checkpoint is a deterministic search
ledger, not a probability certificate.
