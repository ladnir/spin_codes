# Exact weight-18 shell and a limitation of the coupled model

Updated: 2026-09-05.

The new code-specific identity is A18(HQ-perp)=0. Adding it to the previous
weight-14 and weight-16 identities gives the deterministic cap
A38(C) <= 226,659,579,388,193. The full original RandomStepConv-M22 bound has
margin 37.57061695409 bits. The 40-bit goal remains open.

## Completeness from a known supercode

Use the codes defined in BCH_HULL_DUAL_LOW_SHELLS.md: A is the dimension-69
hull of the published code L, D=A-perp, and HQ-perp is a subcode of D.
The known spectrum of A gives A18(D)=6,332,160 by exact MacWilliams transform.
Both D and HQ-perp are invariant under the affine coordinate maps x -> ax+b
over F256, with a nonzero. This invariance is checked on full binary bases.

The candidate search forms symmetric differences of three affine
three-flats. It uses exact syndrome matching and a necessary intersection
size filter. These choices are search heuristics, not assumptions about
the form of every low-weight word. An initial linear-only search supplied
49 affine orbits, which did not exhaust the shell. Translating the component
flats supplied 25,492 distinct candidate words.

The verifier checks that each candidate has weight 18 and lies in D.
For each as-yet uncovered candidate it explicitly generates all 65,280
affine images and deduplicates them. Removing all covered candidates
ensures that the resulting orbits are disjoint: intersecting group orbits
would be identical and their later seed would already have been removed.
The resulting 97 orbits contain exactly 6,332,160 distinct words. Equality
with A18(D) certifies completeness, independently of the candidate search.

No orbit representative lies in HQ-perp. Affine invariance then excludes
every word in all 97 orbits, proving A18(HQ-perp)=0. The published exact
spectrum remains an explicit mathematical input; no randomized counting
guarantee is used.

## Effect on the bound

With A(HQ)=r+255s and the length-256 Krawtchouk transform K, the new identity
is (K(r+255s))18=0. The resulting small model has 196 variables and 1154
normalized rows. Rational primal and dual checks certify its optimum and
the cap stated above. The full failure bound retains the previously
certified contributions from all other weights and occupations.

## The full intersection constraints are still insufficient

Two separate feasibility probes combine all the intersection constraints
with, respectively, the exact weight-14/16 identities and the exact
weight-14/16/18 identities. Each probe adds h38 >= 5e12 only to seek a
hypothetical spectrum. This auxiliary inequality is not a BCH property and
must never be used as an upper-bound constraint.

For each model, an exact rational feasible point satisfies every constraint,
including all 1475 raw original intersection rows after restoring the
eliminated fixed coordinates. The final augmented model has 1483 raw rows.
At each point, 31 h38 times the certified lower coefficient for the true
weight-38 inner tail exceeds 2^-40. Thus these relaxations cannot establish
the target even with exact inner-tail evaluation. This is a limitation of
the constraints and first-moment route, not an actual BCH spectrum or a
lower bound on the actual code's failure probability.

## Replay

Run sequentially:

    python -B code/certify_bch_hull_affine_shell.py 18 hull_weight18_affine_triples bch256_hull_dual_weight18_complete --verify
    python -B code/audit_bch_extended_hull_cap.py prepare_bch_hull_dual18_probe oa21_hull_dual16_count_probe --verify
    python -B code/bch_full_hull_feasibility_probe.py verify prepare_bch_hull_dual16_count_probe hull_full_exact16_feasibility
    python -B code/bch_full_hull_feasibility_probe.py verify prepare_bch_hull_dual18_probe hull_full_exact18_feasibility

The complete shell receipt is generated/bch256_hull_dual_weight18_complete.json;
the earlier receipt without the `_complete` suffix is only partial. The
new cap certificate is generated/oa21_hull_dual18_probe/audit.json. All
hashed scripts, inputs, and receipts are frozen.

Next: obtain additional code-specific information, testing its effect on
the bound before undertaking expensive enumeration. Do not keep iterating
cuts from a relaxation now known to be insufficient.
