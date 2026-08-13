# Independent audit and brainstorming request: grouped-permutation distance proof

You are being asked for an independent second opinion on a computer-assisted
coding-theory proof under active development.  Please audit the existing
argument, identify any hidden assumptions or unsound reductions, and propose
better proof strategies.  Do not assume that the team's current preferred
direction is correct.

## Objective

The construction is intended for OT generation.  It encodes a vector of
`GF(2^128)` elements through the transpose of a binary linear code.  The proof
problem itself is binary: certify that the frozen randomized binary code has
relative distance at least `0.09` at `N=2^21`, with total first-moment bad-event
probability at most `2^-40`.

The performance-motivated construction uses packet/group permutations.  The
deployment target is group size `g=8`; `g=2` and `g=4` are proof-development
ladder rungs.  BCH inner and outer blocks remain `[128,64,22]`; shrinking the
code blocks is not part of the proposal.

## What is currently proved

### `g=2`

There is an end-to-end outward-certified integer-profile proof.  Its final
cell-local union interval is

```text
[-62.71954806419136023328, -62.71947138524163556114].
```

Thus it gives 62.7194713852 bits of total security, or 22.7194713852 bits
beyond the requested 40.  The older global-maximum allocation gave only
40.0941387824 bits.  The improvement used the same construction, witnesses,
and exact hybrid triangle/singleton cover; it sharpened the union accounting.

Please independently assess whether the hybrid convex-triangle/singleton
argument, exact profile ownership/completeness, fixed-mixture convexity, and
cell-local union accounting are sound.

### Exact `g=4` infrastructure

The five-class profile is `a=(a0,...,a4)` with total atom count `M=524288`
and physical-weight constraint `sum(j*a_j)>=21`.  Exact positive supports
partition the integer profiles.  There are 30 feasible supports; `{0}` is
empty.

For each support, active coordinates are shifted by `a_j=1+b_j`.  The code
constructs the exact integer hull, including the nontrivial low-weight kink
facets for zero-containing supports.  For full support the independently
reconstructed integer hull has 15 vertices and 9 maximal facets.

The producer builds an exact rational regular triangulation.  Floating Qhull
is used only to propose lifted facets.  Every retained facet is reconstructed
over exact rationals and checked against all anchors.  The verifier separately
checks hull vertices/facets, extremality, determinant signs, paired facets,
boundary facets, orientations, and exact total volume.

Fixed mixtures use exact rational weights in the emitted ledger.  The outward
verifier independently reloads SHA-bound witness sources, hardens each named
witness, enforces whole-cell support eligibility, evaluates every cell vertex,
and recomputes the cell-local union.  It does not trust the floating LP or its
candidate-selection procedure.

## What is not proved

There is no end-to-end `g=4` theorem yet, and therefore no `g=8` theorem.

The current best full-support `g=4` diagnostic has:

```text
full-support anchors                    1653
exact cells                            37750
cell-local log2 bad-event upper   +67252.4553
largest cell branch               +67198.7278
```

The target is at most `-40`, so this is still far from a certificate.

However, every deliberately targeted pointwise profile has closed.  Several
rounds of legal witness-transition profiles remained safe; the hardest
observed transition still had about 1944 bits of pointwise margin.  This is
evidence, not a theorem.

The original global mesh was badly unsound as a sharp bound because cells
crossed exact-support faces: zero-fugacity witnesses become infinite off their
support.  Exact support stratification removed that defect.  Within full
support, unconstrained regular simplices still cross regions where different
fixed witnesses or mixtures are useful, making convex interpolation extremely
loose.

Structured refinement results were:

```text
method                                      worst branch (approximately)
initial full-support exact mesh                         +922259
128 worst-cell centroids                                +882949
grid on recurring pure (1,2,3) face                    +497145
coarse feasible full-simplex grid                       +301859
first dominance-transition round                       +148586
subsequent transition rounds                            +67199
```

A proposed local subdivision of only the complete star of a pure triangle was
rejected correctly: edge midpoints/grid points produced hanging facets in
neighboring cells containing an edge but not the whole triangle.  The exact
mesh audit caught the first unpaired facet.  Current refinements therefore add
anchors globally and rebuild the exact regular triangulation.

The finite-atlas minimax mixture LP has also been audited.  Per-vertex top-k
selection is not complete: a balanced witness can be useful without winning
at any vertex.  Dual column generation now scans every eligible atlas column.
It matches literal all-column LPs on deterministic tests, uses at most five
active witnesses per `g=4` cell, and reports a final pricing gap below
`8.5e-5` bits in the current run.  Hence ordinary finite-atlas candidate
truncation is no longer the apparent owner of the `+67k` residual.

## Current proposed direction (please challenge it)

All fixed witness bounds share a common profile-normalization term:

```text
F_i(a) = C_i - <q_i,a> - log2 Q(a).
```

Therefore differences `F_i-F_j` are affine.  The current proposal is to replace
unconstrained point refinement by a globally conforming polyhedral complex cut
along selected affine equality/dominance hyperplanes of important fixed
witnesses or fixed rational mixtures.  Each final polytope would then be
triangulated exactly and independently certified at its vertices by one fixed
witness/mixture.  Hyperplane dominance would be discovery geometry only, not a
trusted proof assertion.

We need to know whether this is actually the best route, whether it can scale
to `g=8` (nine profile classes and 510 feasible exact supports), and whether a
more analytic proof can avoid an enormous chamber complex.

## Questions for your audit

1. Is the `g=2` end-to-end proof genuinely sound?  Identify any exact place
   where expectation, independence, convexity, support eligibility, ownership,
   or union accounting might be misused.
2. Is exact-support stratification for `g=4` mathematically sufficient and
   correctly formulated?  Audit the shifted integer-hull argument, especially
   zero-containing supports and the weight-21 clipping/kink facets.
3. Is the regular-triangulation producer/verifier separation adequate for a
   computer-assisted proof?  Look for ways floating proposal generation could
   still contaminate the exact conclusion.
4. Does the fixed rational mixture convexity argument really certify every
   profile in a simplex from its vertices?  Check the Gamma extension and
   common-domain requirements.
5. Does the cell-local weighted union ledger safely handle shared cell
   boundaries and support strata?  Is any event being omitted or counted with
   the wrong multiplicity?
6. Given that targeted pointwise profiles close but convex cells remain very
   loose, what is the best next proof representation?
7. Critique the proposed affine dominance-hyperplane complex.  Estimate its
   likely combinatorial growth at `g=4` and `g=8`, and suggest a constrained
   subdivision/cutting algorithm whose global conformity can be checked.
8. Look for alternatives that could dominate the geometry approach, such as:
   stronger analytic curvature/Hessian bounds within cells; Fenchel/entropy
   duality; branch-and-bound with interval Taylor models; support-specific
   global inequalities; integrating/summing witness bounds without a profile
   union; exact chamber arrangements; or a different enumerator/composition
   theorem.
9. What additional experiment would most efficiently distinguish “the theorem
   is true but the cover is loose” from “there is a genuine difficult profile
   family”?
10. Give a prioritized path to `g=8`, including explicit stop/go criteria.

Please label every conclusion as one of:

- **proved from the supplied artifacts/code**;
- **plausible but requiring an additional check**;
- **diagnostic evidence only**; or
- **incorrect/unsound**.

We value counterexamples and objections more than agreement.  If an existing
claim is not independently checkable from the repository, state exactly what
artifact or invariant is missing.

## Recommended reading order

1. `PROOF_STATUS.md` — current proof ledger and precise numerical status.
2. `explorations/g2_triangular_cover_math.md` — `g=2` convex/hybrid theorem.
3. `explorations/g4_delaunay_cover_math.md` — `g=4` geometry, counting, and
   support-hull audit notes.
4. `scripts/certify_packet_group_triangle_ledger.py` — outward `g=2` hardener.
5. `scripts/certify_packet_group_g4_anchor_mesh.py` — independent exact `g=4`
   geometry and outward verifier.
6. `scripts/build_packet_group_g4_support_regular_mesh.py` — support-stratified
   exact mesh producer.
7. `scripts/reattach_packet_group_g4_support_mixtures.py` and
   `scripts/probe_packet_group_g4_stellar_cover.py` — mixture optimization and
   finite-atlas column generation.
8. `scripts/generate_packet_group_g4_transition_anchors.py` — dominance
   transition discovery; this is diagnostic, not certificate logic.
9. The relevant manuscript files: `localCodeCertificate.tex`,
   `localCodeInner.tex`, `localCodeProjections.tex`, and
   `localCodeStructured.tex`.

Avoid relying on obsolete Singer/PAP-era probes or `mtx.tex`; they are not part
of the current claimed proof chain.  Also do not treat diagnostic binary64 LP
scores, sampled profiles, or Qhull adjacency as certificate evidence.
