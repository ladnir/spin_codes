# Independent second opinion: grouped-permutation distance proof

Auditor: Claude (Fable 5), 2026-08-12.  Independent of the team's prior
sessions.  Everything below was re-derived from the repository files listed at
the end; I did not execute any verifier, so "proved from the supplied
artifacts/code" means "the committed code implements a sound argument as
written," not "I reproduced the run."

Labels used, per the request:

- **PROVED** — proved from the supplied artifacts/code;
- **PLAUSIBLE** — plausible but requiring an additional check;
- **DIAGNOSTIC** — diagnostic evidence only;
- **UNSOUND** — incorrect/unsound as stated.

---

## 0. Summary of headline verdicts

1. The `g=2` end-to-end integer-profile certificate architecture is sound.  I
   independently re-derived the convexity lemma, the fixed-mixture rule, the
   hybrid ownership theorem, the Pick lattice counts, the 2-chain mesh
   conformity audit, the cell-local weighted union, and the Collatz
   composition, and checked the verifier code against each.  **PROVED**, with
   three caveats below (F3, F6, F7), of which one (F3, float underflow in the
   "outward" DP) is a strict-soundness defect that is almost certainly
   immaterial but should be fixed before the result is called fully outward.
2. The single most load-bearing *semantic* assumption — that conditional on a
   packet-weight profile the interleaved word is uniform on the whole
   weight-profile class, which is what justifies `Q(a)` — is justified by the
   lane-bijection factorization stated in
   `explorations/riffle_group_chain_proof.md:25-32,116-123`.  It is currently
   documented only in an exploration note.  **PLAUSIBLE** until the frozen
   construction spec (and the C++ implementation) is checked to actually
   sample independent within-packet lane orders.  If the deployed interleaver
   permutes packets *without* randomizing internal packet order, every
   profile-normalized certificate in this program is bounding the wrong
   ensemble.  This deserves a one-page lemma in the manuscript and an
   implementation cross-check, not a note.
3. Exact-support stratification for `g=4` is correctly formulated and the
   verifier genuinely re-derives it.  I independently proved the
   integer-hull candidate lemma the verifier relies on (the "corner band"
   midpoint argument) — it is correct.  **PROVED**.
4. The producer/verifier separation is adequate: floating Qhull, floating LPs,
   and power iterations only *propose* objects (facets, mixture supports,
   Collatz vectors); every retained object is re-proved exactly or outward.
   Failure modes of the discovery layer make the verifier reject, never
   accept.  **PROVED** (architecture), with the F3 caveat on the "outward"
   float DP layer.
5. The `+67k` full-support `g=4` residual is a *minimax-game artifact*, not a
   witness-atlas or pointwise deficiency.  My analysis (§7) says the proposed
   dominance-chamber direction is essentially right, and moreover that it has
   a clean termination property the team has not stated explicitly: **if the
   finite atlas's pointwise lower envelope closes on the whole stratum, then
   the exact dominance-chamber complex with one owner witness per chamber
   closes with no further refinement.**  The open question is therefore not
   "how much refinement" but "does the envelope close" — a question the
   chamber complex itself certifies at its vertices.  Estimated chamber
   complexity is manageable at `g=4` and prohibitive at `g=8` without witness
   compression (§7.3).
6. **Recommended plan modification (§7.5): do not build a globally
   conforming chamber complex at all.**  The cell-local union ledger only
   needs a valid upper bound per *existing* cell, so the dominance geometry
   can be realized as an independent interior BSP decomposition inside each
   failing cell of the frozen, already-audited mesh.  This keeps every
   existing conformity proof closed forever, removes all pairing/orientation
   obligations for the new cuts, is embarrassingly parallel and fully
   adaptive, and has the same termination property as the global complex.
   The outward log-sum-exp of 37,750 terms can exceed its largest term by at
   most `log2 37750 ≈ 15.2` bits, so the residual is necessarily owned by
   the top of the cell queue; a shallow BSP on the worst few hundred cells
   should collapse the residual by orders of magnitude before any
   systematic pass is needed.  The decisive first experiment (§9) is one
   day of work on the single worst cell.

---

## 1. Q1 — soundness of the `g=2` end-to-end proof

### 1.1 What I verified line-by-line (all PROVED)

**Convexity and mixtures** (`explorations/g2_triangular_cover_math.md`,
`certify_packet_group_triangle_ledger.py`):

- Lemma 1 (trigamma) is correct; `F_W` is convex on the mass plane, so the
  max over a closed cell is attained at a vertex.
- The mixture rule (min ≤ fixed convex average) is correctly *not* Jensen and
  correctly requires nonnegative weights summing to exactly one; the verifier
  enforces exact rational weights with `sum == 1`
  (`_mixture_from_row`, `parse_mixture`) and hard-rejects any vertex-varying
  mixture schema (`forbidden` keys, `len(set(resolved)) != 1`).
- Mixtures may mix `subtract_normalization=True` inner witnesses with
  outer-only affine witnesses: the combined function is
  `affine + μ·(−log2 Q)` with `μ ∈ [0,1]`, still convex.  Sound.
- Support eligibility: `evaluate_vertex` raises when a component has zero
  fugacity on a class that is *positive at any vertex*.  Because each
  coordinate `a_j` is affine, `a_j = 0` at all three vertices forces
  `a_j ≡ 0` on the closed triangle, so the vertex check is exactly the
  "identically zero on the cell" condition required by the math note.  The
  `g=4` verifier additionally enforces stratum-level eligibility
  (`require_support_eligible`).  Sound.

**Exact geometry** (`_mesh_covers_polygon`):

- The audit is a rigorous 2-chain/degree argument: all triangles CCW inside
  the hull, exact primitive-line edge splitting handles conforming
  T-junctions, every interior edge piece has multiplicity 2 with opposite
  orientations, every boundary piece has multiplicity 1 and exactly tiles a
  hull edge (rational interval sweep), and doubled areas sum to the hull
  area.  Crossing any interior edge changes the coverage count by 0 and the
  hull boundary by 1, so coverage ≡ 1 a.e.  Combined with vertex membership
  and nondegeneracy this proves an exact partition up to shared boundaries.
  Sound.
- Both accepted domains cover `Λ`: the pentagon is the integer hull of `Λ`
  (I checked all five vertices, including the `(1,10)` kink), and every
  lattice point of the pentagon satisfies the hull-valid inequality
  `a1+2a2 ≥ 21`, hence is in `Λ` — so pentagon lattice counts contain no
  sub-threshold slack at all.  The full-support quadrilateral + three edge
  segments alternative also covers `Λ` exactly (two-zero-coordinate profiles
  are segment endpoints).

**Lattice counting and terminal cells**:

- `_closed_triangle_lattice_count` implements Pick correctly:
  `I + B = (2A + B + 2)/2`.  `a0 = M − a1 − a2 ≥ 0` holds inside by
  convexity of the vertex constraints.
- `_triangle_lattice_profiles` scanline: I checked the half-plane feasibility
  logic (including the vertical-edge case and floor/ceil of `Fraction`s) —
  exact.
- Terminal cover: every terminal-cell lattice point must appear in the
  deduplicated singleton assignments (`missing = expected - assigned`
  rejected), digests bind counts, and each singleton is outward-evaluated
  individually.  Duplicate profile/mixture pairs would only overcount —
  safe.

**Union accounting**:

- Cell-local ledger: `n_r 2^{U_r}` per closed triangle plus deduplicated
  singletons, outward log-sum-exp (`_log2_sum_exp` over guarded Decimal
  exp/ln with relative pads).  Closed-boundary duplicates retained as a
  nonnegative overcount — valid per equation (12) of the math note.
  Profiles, not leaves/witnesses/components, are the union events.  Sound.
- The final claim `E[Z_d] ≤ 2^{-62.71947…}` at `d = 188743` implies
  `d_min ≥ 188744`, and `188744/2^21 > 0.09`.  Note the inner map of this
  construction is *always* a bijection for fixed state permutations
  (`riffle_group_chain_proof.md:110-114`), so no separate noninjectivity
  event is needed — unlike the RM/EBCH line, where the manuscript correctly
  bounds the joint event.

**Witness hardening (Collatz composition)**:

- `witness()` power iteration only proposes a positive vector; soundness is
  the Collatz–Wielandt inequality: `image ≥ A v` outward, so
  `λ̂ = max(image/v)` satisfies `A v ≤ λ̂ v` and
  `1ᵀAᴮe₀-type sums ≤ (max 1/v_q) λ̂ᴮ v₀`.  I re-derived this chain and it
  is correct, including monotonicity of the nonlinear upper operator `F`.
- `outward_transport_image`: the two-pointer greedy is the exact optimum of
  the rank-one transportation *relaxation* (row sums ≤ enlarged supplies,
  column sums ≤ split caps, nonnegative profits), and the true mass vector is
  feasible for it, so the greedy value dominates.  The
  "supplies exceed sink capacity" case is handled correctly (leftover mass
  contributes zero in the ≤-row relaxation while the true equality-row mass
  always fits).  Zero caps with positive mass hard-exit.  Exact `Fraction`
  arithmetic throughout.  Sound.
- The drive-stratified lemma itself
  (`explorations/packet8_drive_stratified_transfer.md`): I verified the
  shared-column inequality (2) (disjointness over drive weights under the
  accumulator bijection), the rank-one rearrangement bound, and the
  Chernoff direction of the `pole^{-D}` term.  Sound.

**Interval layer** (`outward_log2.py`):

- Directed Decimal rounding, min/max recomputed under each rounding mode,
  Robbins `ln n!` bounds with Archimedean rational `π` bounds, and a
  `self_check`.  The `1e-110` pad dominates the 130-digit rounding unit for
  every magnitude that occurs (largest `ln` arguments here are ≲ 10^5·digits,
  giving ulps ≲ 10^{-125}).  Sound.
- `outward_normalization` implements `Q(a) = M!/∏a_j! · ∏ C(g,j)^{a_j}`,
  i.e. includes the `2^{a1}` (g=2) and `4^{a1+a3}6^{a2}` (g=4) class factors
  exactly once.  The fractional-vertex variant uses the chord upper bound on
  `ln Γ` (valid by log-convexity) on the side that makes `F` larger, and a
  crude but valid `> −1` floor on the unused side.  Sound.

### 1.2 Findings (caveats)

**F1 (PLAUSIBLE — the most important open check).  Ensemble uniformity on
the profile class.**  All profile certificates divide by `Q(a)`, the size of
the *weight-profile class*, which presumes that conditional on the profile
the interleaved word is uniform over that class.  A uniform permutation of
packets alone does *not* give this: its orbit is the concrete-value multiset
class, which is smaller.  The justification is the factorization of the
outer lane bijection into (packet assignment) × (independent within-packet
orders) (`riffle_group_chain_proof.md:25-32`), which makes internal packet
patterns independent and uniform within each weight class
(lines 116-123) and hence the conditional law uniform on the profile class.
This argument is correct — but it lives in an exploration note, is nowhere
stated as a numbered lemma consumed by the `g=2` theorem, and I could not
check that the deployed C++ construction actually samples those within-packet
orders.  Action: promote to a manuscript lemma ("interleaver acts
transitively on each weight-profile class") and add an implementation
conformance note.  Everything downstream is conditioned on this.

**F2 (PLAUSIBLE).  Outer-witness semantics.**  The mixture theorem requires
all components to bound the *same* `Z(a)`.  This is a semantic obligation the
verifier cannot check mechanically.  I verified the composition shape and two
branches end-to-end: the `full_bijection` witness
(`Z(a) ≤ 2^K·Vol(N,d)/Q(a)` with `Vol(N,d) ≤ 11^N/10^{N−D}` — exact), and
the three-band linear-BL polytope lemma (the algebra
`Σ p_b dim P_b(V) ≥ dim V` whenever the two smallest coefficients sum to at
least 1, given the certified rank-64 direct-sum decomposition — correct).
The finite-field linear Brascamp–Lieb inequality itself is standard.  What I
did not re-derive: the full `packet_group_outer_profile.py` moment algebra
(`symmetric_s2`, band norms, `256·42/256·86` exponents) and the exact
graph-puncture Maclaurin branch.  These are shared, hashed modules; they
deserve one dedicated audit pass of their own.

**F3 (UNSOUND as stated; almost surely immaterial; cheap fix).  The
"outward" float DP has no underflow guard.**  `_up_mul` returns `raw`
unchanged when `raw == 0.0`, and `_up_add` drops `right == 0.0`
contributions.  IEEE-754 products of ≥ ~1075 magnitude-bits underflow to
exactly `0.0`.  With sharp support-local witnesses (fugacities `≲ 2^{-34}`)
a 32-atom path product can underflow, silently zeroing a structurally
positive histogram entry in what is claimed as an entrywise *upper* bound;
the enclosure self-check (`histograms_upper >= diagnostic_histograms`)
cannot catch it because the round-to-nearest diagnostic underflows
identically, and `outward_transport_image` then skips the stratum
(`if histogram == 0.0: continue`).  The induced error is bounded by the
dropped mass itself (≲ 2^{-1000} against margins of ≥ 0.059 bits), so no
realistic artifact is wrong — but the certificate's "one-ulp upward" claim
is strictly false in the underflow regime, and one certified local margin is
as thin as 0.0596 bits.  Fix: run a boolean reachability DP over the exact
structural-zero pattern of `compatible_pair_polynomials` (whose zeros I
verified are genuinely structural: a nonnegative float sum is 0 iff all
terms are 0, and `multiplicity × fugacity` cannot underflow since
`multiplicity ≥ 1`) and assert the float DP is strictly positive wherever
the boolean DP is reachable; or scale the DP per atom step.

**F4 (PROVED, positive finding).  Fatal-shortcut enforcement.**  Each fatal
shortcut listed in the two math notes has an explicit code rejection: per-
vertex witness selection, vertex-varying mixtures, binary64 weights, zero-
fugacity leakage, degenerate cells, mismatched reported counts/digests, and
nonconforming local refinements (the hanging-facet rejection is real:
`verify_mesh` facet multiplicity/orientation audit).

**F5 (DIAGNOSTIC).  Artifacts not reproduced.**  I did not re-run the
outward parallel run or check the SHA-256 values in `PROOF_STATUS.md`
against `out/`.  The chain "code is sound ⇒ the recorded interval
`[-62.7195480…, -62.7194713…]` is certified" additionally requires the
recorded artifacts to be the ones this code produced.

**F6 (PLAUSIBLE).  Frozen constants.**  `K`, `D`, `INNER_BLOCKS`,
`GRAPH_REPLACEMENTS`, and the 42/86 band sizes are imported from shared
modules I did not re-derive against the frozen construction.  A one-page
"constant provenance" table would close this.

**F7 (PLAUSIBLE).  EBCH spectrum provenance.**  Acknowledged in the
manuscript for the RM/EBCH rows; the same declared-input status applies to
every packet-group certificate through `split_cap_table` /
`EBCH128_64.wd` / the systematic split slices.  Cheap hardening: verify the
committed table is a MacWilliams fixed point / transform-consistent with its
dual and matches the published literature value; better, regenerate the
split-slice table from the committed generator with an independent method
and compare digests.

---

## 2. Q2 — `g=4` exact-support stratification

**PROVED**, with the following independently re-derived pieces:

- The support partition is exact and the shifted formulation is correct:
  for `0 ∉ S` the weight-21 cut is redundant (`Σ j a_j ≥ Σ_{j∈S} a_j = M`),
  `|Λ_S| = C(M−1,|S|−1)` is stars-and-bars, and for `0 ∈ S` the continuous
  hull (17) has *exactly* `Λ_S` as lattice points (integral change of
  variables; `b_0 ≥ 0 ⇔ Σy ≤ R`).  No profile can be omitted.
- The verifier's integer-hull certification is genuinely independent of the
  producer: it (i) enumerates its own candidate set
  (`integer_hull_candidates`), (ii) rebuilds every supporting facet of the
  claimed vertices by exact determinant sign tests, (iii) checks every
  candidate lies inside those facets, and (iv) checks every claimed vertex is
  a candidate and is facet-extreme.  I proved the candidate lemma the code
  relies on: an integer point with `b_0 > 0` and physical shift
  `≥ threshold + max(S)` is the midpoint of `y ± e_i` for any occupied
  positive class `i` (both neighbors feasible), and hull vertices on the
  `b_0 = 0` face must be vertices of that face's plain simplex, i.e. pure.
  Hence all hull vertices are pure or lie in the enumerated corner band.
  Correct, including the `(2,5)` kink for `P={1,3}` which I spot-checked.
- The kink-facet table in `explorations/g4_delaunay_cover_math.md` is
  therefore *not* in the trusted base (the verifier reconstructs facets);
  label the table itself PLAUSIBLE, consequence-free.
- `verify_support_stratum_mesh` / `verify_mesh` implement the 4-D analogue
  of the 2-D conformity audit: exact determinants, orientation-signed facet
  pairing, unpaired facets forced onto exact hull-boundary hyperplanes, and
  exact total volume against an independently constructed pulling
  triangulation of the certified hull (or the closed-form
  `M^4 − 21^4/4!` for the flat clipped simplex, which I verified).  Since
  all anchors are checked inside the hull, the degree argument gives an
  exact cover a.e.; hanging facets are structurally rejected (and were, per
  the `(0,1,4,15)` episode).  Sound.
- Fractional hull vertices use rigorous Gamma chords (F above); integer
  anchors use exact factorials.  The per-cell lattice bound (best
  coordinate box over the five omitted-coordinate charts, capped by the
  global count) is a valid upper bound; zero-width boxes correctly prove
  empty cells.

One formulation caveat: the strata meshes triangulate the *integer hull* of
each stratum, which is correct for covering `Λ_S`, but cells of a
zero-containing stratum may contain lattice points *below* the stratum's own
continuous cut only if those points are in the integer hull — they are not,
since the covering facets are valid inequalities.  No slack, no omission.
Verified.

---

## 3. Q3 — producer/verifier separation

**PROVED (adequate)**.  The trust boundary is clean:

- Qhull (`Qx Qt`/`Qc`, and the joggled diagnostic meshes) only proposes
  candidate lower facets.  `exact_lower_facet` re-solves the interpolating
  hyperplane over `Fraction` and requires *strict* below-ness of every other
  lifted anchor; the deterministic rational lifting
  (`‖x‖² + prime-residue perturbation`) breaks ties exactly.  A missing
  facet cannot be silently accepted: the downstream conformity audit
  (facet pairing + volume) fails.  A wrong facet cannot be accepted:
  `exact_lower_facet` rejects.  Failure directions all point to rejection.
- Binary64 witness parameters (`fugacities`, `pole`, log variables) are
  converted by `Fraction.from_float` and thereafter treated as *the* exact
  parameters; tuning quality affects only margins, never soundness.
- LP mixture weights are re-parsed as exact rationals, re-summed to exactly
  one, and every vertex inequality is re-evaluated outward.  The floating LP
  selects supports only.
- The residual floating channel is the "outward" float DP itself (F3).
  That is the only place where a floating computation is *load-bearing* for
  an inequality rather than a proposal, and it is the place to spend the
  next hardening effort.
- Minor: `certify_packet_group_g4_anchor_mesh.py` mutates
  `g2.GROUP_BITS`/`g2.BLOCK_ATOMS` at runtime (`configure_generic_hardener`)
  and clears the normalization cache.  This works (Python late-binds module
  globals, and the cache clear prevents cross-`g` contamination), but it is
  fragile under future refactors; pass `g` explicitly.

---

## 4. Q4 — fixed-mixture convexity and the Gamma extension

**PROVED** as implemented, with the standing semantic obligation F2:

- The Gamma extension is only ever *evaluated* at cell vertices; convexity
  does the interior work.  Fractional vertices use outward chords on the
  correct side.  Integer vertices are exact.
- Mixture components are evaluated unconditionally at every vertex; a
  support-ineligible component poisons the whole cell (raise), never a
  silent drop.  Correct per the math notes.
- Mixed normalization flags keep convexity (coefficient of the convex term
  in `[0,1]`).
- Caratheodory sparsity (≤3 / ≤5 components) is an efficiency fact, not a
  soundness requirement; the verifier accepts any exact mixture.
- The "common domain / same `Z(a)`" requirement is enforced by convention
  and code structure (all four branch families — linear-BL outer × drive
  inner, total-weight outer × drive inner, exact graph-puncture outer ×
  drive inner, full bijection — bound `E[#bad codewords with profile a]`),
  but no machine check exists nor can one.  Keep the family list closed and
  documented.

---

## 5. Q5 — union ledger, shared boundaries, multiplicities

**PROVED**.  Specifically:

- Events are profiles; leaves, witnesses, mixture components, and duplicate
  closed-cell memberships are proof devices.  The implementation counts every
  closed-cell boundary profile at least once and possibly more — a
  nonnegative overcount, never an omission.
- Omission is excluded by: exact geometric cover of a set containing `Λ`
  (pentagon = integer hull at `g=2`; 31 disjoint stratum hulls at `g=4`),
  plus exhaustive terminal-cell lattice matching at `g=2`.
- Cross-stratum double counting is impossible (exact supports partition
  `Λ`); within-stratum duplicates are safe.
- The `g=4` `union_accounting_mode` gate forces cell-local accounting for
  support-stratified ledgers, preventing an accidental global-max union with
  a stale count.  Good defensive design.
- One wording nit: `PROOF_STATUS.md` says shared boundaries are
  "deliberately counted more than once, which is a safe nonnegative
  overcount" — true, but note the *terminal* singletons are globally
  deduplicated while also potentially lying inside adjacent convex cells;
  that, too, is only an overcount.  No error.

---

## 6. Q6 — why cells are loose while points close (diagnosis)

This is the strategic crux, so I state it precisely.

Write `F_i(a) = A_i(a) − log2 Q(a)` with `A_i` affine.  The pointwise bound
is the envelope `E(a) = min_i A_i(a) − log2 Q(a)` — a *concave* piecewise
affine part plus a convex part.  A convex cell with fixed mixture `λ` is
charged `max_k Σ_i λ_i F_i(v_k)`; optimizing `λ` gives the matrix-game value

```
G(cell) = min_λ max_k Σ_i λ_i F_i(v_k) = max_μ min_i Σ_k μ_k F_i(v_k),
```

by LP duality, where `μ` ranges over probability vectors on the vertices.
By convexity `Σ_k μ_k F_i(v_k) ≥ F_i(a_μ)` at the barycenter `a_μ`, so

```
G(cell) ≥ E(a_μ*) ,  and moreover
G(cell) ≤ max_μ [ min_i A_i(a_μ) ] − min over cell of log2 Q + Jensen terms.
```

Two distinct looseness sources follow:

1. **Envelope-vs-mixture gap**: if different witnesses own different parts
   of the cell, `max_μ min_i Σ μ_k A_i(v_k)` sits far above
   `min_i A_i` anywhere — this grows linearly with (cell diameter) ×
   (spread of active charge vectors), i.e. thousands of bits for long cells.
2. **Jensen gap of `−log2 Q`**: bounded by
   `½ Σ_j diam_j² · ψ₁(l_j+1)/ln 2 ≤ ½ Σ_j diam_j²/(l_j ln 2)` on a cell
   with coordinate lower bounds `l_j` — hundreds of thousands of bits for
   cells spanning pure corners, matching the observed `+922k → +67k`
   trajectory as cells shrank.

Both terms vanish for a cell owned by a single witness whose dominance
region contains it: then the optimal `λ` is pure, `G(cell) = max_k F_{i*}(v_k)`,
which is the *exact* maximum of `F_{i*}` over the cell — no interpolation
slack at all.  The column-generation result (pricing gap `< 8.5e-5` bits)
already told you truncation is not the owner; this analysis says the owner is
the vertex-averaging in the game value, which no amount of witness addition
fixes without geometry aligned to dominance.  **DIAGNOSTIC** (analysis), but
the inequalities above are elementary and checkable.

**Key consequence (the termination property).**  Suppose the atlas envelope
closes pointwise on a stratum: `E(a) ≤ τ(a)` for all real `a` in the hull.
Build the exact chamber complex of the minimization diagram of `{A_i}`
(restricted to the hull), assign each chamber its owner `i*`, and triangulate
chambers if desired.  Every chamber vertex `v` satisfies
`F_{i*}(v) = E(v) ≤ τ`, and convexity extends this over each chamber.  So the
cover closes *immediately*, with zero refinement waves, if and only if the
envelope closes at the (finitely many) chamber vertices.  The chamber complex
is simultaneously the certificate and the exhaustive search for pointwise
counterexamples: any chamber vertex where the envelope fails is a concrete
profile where the atlas is genuinely insufficient (new witness or a real
hard region).  This reframes Q9 (see §9).

---

## 7. Q7 — critique of the affine dominance-hyperplane proposal

**Verdict: the proposal is sound and the dominance geometry is the right
lever at `g=4` — but its best realization is *not* a globally conforming
complex.  §7.1–7.4 assess the proposal as written; §7.5 gives the
recommended per-cell variant, which achieves the same certified numbers
without ever reopening the global geometry proof.**

### 7.1 Soundness

Dominance hyperplanes as discovery only, with every final cell certified by
one fixed witness/mixture at its vertices, adds no new trust obligations —
the existing vertex-inequality + conformity machinery already covers it.
Correct as proposed.

### 7.2 Two modifications

1. **Cut globally, not locally.**  The rejected pure-triangle star episode
   generalizes: local subdivision against a background complex creates
   hanging facets.  Cutting the *entire stratum complex* with one affine
   hyperplane at a time preserves conformity by construction (every cell
   crossed is split by the same plane, facets stay matched), and the
   existing facet-pairing audit re-verifies it.  Iterate: pick the dominance
   hyperplane separating the two best witnesses in the worst residual cell,
   cut the whole stratum, re-certify.  This avoids both the global
   retriangulation cost and the nonconformity trap.
2. **Assign owners, drop mixtures where possible.**  Inside a dominance
   chamber a pure witness is optimal and its cell bound is exact (no game
   gap).  Mixtures should survive only on cells that deliberately straddle
   chambers (e.g., to avoid tiny slivers).  This also removes most of the
   LP machinery from the certificate path.

### 7.3 Combinatorial growth estimates

The minimization diagram of `W` affine functions on `R^d` is the projection
of the upper envelope of a polyhedron in `R^{d+1}` with `W` facets: by the
Upper Bound Theorem its total face count is `O(W^{⌊(d+1)/2⌋})`.

- `g=4`, full support, `d=4`: `O(W²)` chambers.  With `W` = number of
  *envelope-active* witnesses (those minimal somewhere — prune by LP; likely
  a small fraction of the 1653-anchor atlas, plausibly `10²`), expect
  `10⁴–10⁵` chambers and a similar order of vertices — comparable to the
  current 37,750-cell mesh.  Feasible.
- `g=8`, full support, `d=8`: `O(W⁴)`.  Even `W = 300` gives `~10¹⁰`
  worst-case faces; with 510 supports (most low-dimensional, but the
  handful of high-dimensional strata dominate) this is **not** feasible by
  direct chamber enumeration.  `g=8` therefore needs witness compression
  (few parametric families with analytically controlled envelopes) or a
  stronger local inequality *before* geometry — see §10.

Incremental global cuts (7.2.1) typically realize far fewer cells than the
worst case because only crossed cells split; measure the growth rate on
`g=4` and extrapolate before committing `g=8`.

### 7.4 What could still fail

- The envelope may genuinely not close somewhere (`g=4` history: the
  `(375821,…)` profile needed a *new branch*, not refinement).  The chamber
  complex will surface such points as failing chamber vertices — that is a
  feature, but budget for a discovery loop, not a single pass.
- Chamber vertices are intersections of up to `d` hyperplanes and hull
  facets: exact rational but with large denominators; Gamma-chord
  evaluations stay cheap, but expect heavier rationals than mesh anchors.
- The union ledger needs per-chamber lattice-count uppers; the existing box
  bound works on any polytope via its vertex coordinate ranges.

### 7.5 Recommended variant: per-cell dominance BSP inside the frozen mesh

This is the concrete plan I would execute instead of the global chamber
complex.  It uses the same dominance geometry but relocates it *below* the
cell level, where the union ledger imposes no conformity obligations.

**7.5.1 The observation.**  The cell-local ledger (equation (13)/(14) of the
`g=4` note) needs, for each closed cell `σ` of the verified mesh, only a
number `U_σ` with

```
max over Λ ∩ σ of log2 Z(a)  ≤  U_σ.
```

Nothing in the ledger requires `U_σ` to be produced by one fixed witness or
mixture evaluated at `σ`'s own five vertices.  Any exact interior
decomposition of `σ` may be used, because it is invisible to every other
cell: the mesh conformity audit, the volume audit, and the stratum
disjointness proof all remain exactly as already certified.  The hanging
facet problem that forced global retriangulation (the rejected `(0,1,4,15)`
episode) cannot occur, since no facet of an interior decomposition is ever
matched against a neighboring cell.

**7.5.2 The per-cell certificate (BSP form).**  For a failing cell `σ` with
locally active witness set `I(σ)`:

1. Build a binary space partition: split `σ` by an exact rational affine
   hyperplane (normally a dominance hyperplane `A_i = A_j`, but soundness
   does not require this), split each child recursively, stop at leaves.
   Each split partitions its parent exactly *by construction* — the two
   children are the parent intersected with the two closed half-spaces —
   so completeness of the leaf cover needs no audit beyond checking each
   recorded split.
2. At each leaf `L` (a convex rational polytope: `σ`'s five facets plus the
   accumulated cut half-spaces), choose one owner witness `i_L ∈ I(σ)`
   (discovery may pick the argmin at the leaf barycenter or the minimizer
   of the vertex max; the choice affects only sharpness).
3. Enumerate the leaf's vertices exactly and set

   ```
   U_L = max over vertices v of L of outward-upper F_{i_L}(v),
   U_σ = max over leaves of U_L.
   ```

**7.5.3 Soundness (one-paragraph lemma).**  For every real `a ∈ σ` there is
a leaf `L ∋ a` (BSP leaves cover `σ` by construction).  On `L`,
`log2 Z(a) ≤ min_i F_i(a) ≤ F_{i_L}(a)`, and `F_{i_L}` is convex (Lemma 1 of
the `g2` note, unchanged), so `F_{i_L}(a) ≤ max_{v ∈ vert(L)} F_{i_L}(v) ≤
U_L ≤ U_σ`.  Support eligibility is inherited from the stratum-level rule
already enforced (`require_support_eligible`): every witness usable anywhere
in the stratum is finite on all of it, so leaf owners need no new
eligibility machinery.  Lattice counting stays at the cell level with the
existing box bound (or, for sharper accounting, per-leaf box bounds summed —
also valid, since leaf duplicates on shared cut planes are a nonnegative
overcount exactly as for closed cells).

**7.5.4 Verifier obligations (all local, all existing machinery).**

- each recorded cut hyperplane has exact rational coefficients;
- each leaf's half-space list equals its ancestors' cuts with the recorded
  signs (tree well-formedness — a purely syntactic check);
- exact vertex enumeration of each leaf (small 4-D polytopes: ≤ 5 simplex
  facets + depth-many cuts; rational LP-free vertex enumeration over
  `Fraction` is adequate at these sizes);
- outward evaluation of the owner at each leaf vertex — the existing
  `evaluate_profile` path, including the fractional-Gamma chord for
  non-integer vertices;
- `U_σ` recomputed as the max, and the final ledger log-sum-exp unchanged.

No facet pairing, no orientation audit, no global volume identity, and no
change to any already-passed check.

**7.5.5 Why this dominates the two alternatives on the table.**

- Versus *global conforming hyperplane cuts* (§7.2.1): same certified
  numbers in the limit, but the global variant reopens the stratum geometry
  proof after every cut, serializes the pipeline on retriangulation, and
  spreads cuts into thousands of cells that were already closing.  The BSP
  variant touches only failing cells and leaves every closed cell's
  certificate byte-identical.
- Versus *more atlas growth / LP retuning at current geometry*: the
  column-generation result (final pricing gap `< 8.5e-5` bits) plus the
  game-value analysis (§6) prove this cannot help — the residual is vertex
  averaging, which is a property of the geometry, not of the witness pool.
- Adaptivity matches where the mass is: the outward log-sum-exp of the
  37,750 cell terms exceeds its largest term by at most
  `log2 37750 ≈ 15.2` bits, so the residual total is always dominated by
  the heaviest few cells and drops with them.  (Consistency note on the two
  reported diagnostics: the total `+67252.4553` sits `~53.7` bits above the
  largest *branch* `+67198.7278`, which exceeds the 15.2-bit LSE spread —
  so the ledger's largest *term* must carry roughly 38–54 bits of per-cell
  lattice-count factor on top of its branch.  The worst-cell queue should
  therefore be ordered by `U_σ + log2 n_σ`, not by `U_σ` alone.)  Expect
  the first shallow pass (depth 3–5 on the worst few hundred cells) to
  collapse the total by orders of magnitude.

**7.5.6 Termination and failure semantics.**  As leaves shrink toward
dominance chambers, the owner's vertex max converges to the pointwise
envelope, which the transition diagnostics place ≥ ~1900 bits below target.
Hence, *if* the envelope closes on the stratum, a bounded-depth BSP closes
every cell (same termination property as the global complex, §6).  If some
fully refined leaf still fails, its failing vertex is a concrete profile
where the envelope itself fails — a proof-content gap requiring a new
branch (as with the graph-puncture profile), which no subdivision scheme of
any kind would have fixed.  Both outcomes are actionable; neither is wasted
work.

**7.5.7 Cautions.**

- Leaf vertices are intersections of up to four cut planes with simplex
  facets: exact rationals with materially larger denominators than mesh
  anchors.  The Gamma-chord evaluation handles them, but budget for heavier
  `Fraction` arithmetic; interning/hashing repeated vertices across sibling
  leaves is worth it.
- Cap BSP depth (say 8–10).  In the rare cell where the local decomposition
  does not converge under the cap, fall back to inserting one global anchor
  and re-running the standard retriangulation once per *batch* of such
  cells, not per cell.
- Hyperplane selection is discovery: prefer dominance planes between the two
  best witnesses at the current worst leaf vertex, weighted by residual
  contribution.  A bad choice costs depth, never soundness.
- Keep the per-cell BSP records in the ledger schema (cut list + leaf owner
  list + leaf vertex digests) so the outward verifier replays them without
  rerunning discovery.

---

## 8. Q8 — alternatives assessment

- **Analytic Hessian/curvature bounds within cells**: useful as a *certified
  Jensen-gap bound* — `ψ₁(l+1) ≤ 1/l` gives a cheap per-cell upper bound on
  the interpolation slack, letting the refinement loop prove "this cell can
  never close with any mixture; cut it" instead of discovering it by LP
  failure.  Complementary, not a replacement.  **PLAUSIBLE / cheap win.**
- **Branch-and-bound with interval Taylor models over boxes**: strictly
  worse than vertex evaluation here.  Splitting
  `max(affine) + max(−log2 Q)` over a box pays `min(var affine, var convex)`
  ≈ tens of bits per lattice unit of box width, forcing `≫10¹⁴` boxes.  The
  simplex-vertex evaluation captures both terms jointly and exactly.
  **Not recommended.**
- **Fenchel/entropy duality**: the profile sum is a large-deviations object
  and the envelope-vs-`log2 Q` saddle could in principle be bounded
  analytically, but certifying a continuous sup reintroduces the same
  geometry (or interval B&B).  Worth one exploratory note only if `g=8`
  witness compression succeeds, since a 2-parameter analytic family +
  1-D interval certificates per support would sidestep chamber explosion.
  **PLAUSIBLE, second priority.**
- **Integrating/summing witness bounds without a profile union**: the sum
  `Σ_a 2^{F_i(a)}` for a single fixed witness is max-term dominated (the
  `1/Q(a)` weight concentrates it on the worst corner), which is exactly why
  a single global witness fails by ~10⁶ bits.  Per-chamber Brion/Barvinok
  exponential sums would sharpen `count × max` to a true sum, but that gap
  is only `log2(count) ≲ 70` bits — negligible against the current residual.
  **Not the bottleneck.**
- **Support-specific global inequalities / stronger local inequality
  retaining drive-mass distribution**: this is the right lever for the
  `g≥8` middle-class point-cap obstruction documented in `PROOF_STATUS.md`
  (the `2^{-60}` representative-fugacity tradeoff).  The point-cap
  abstraction throws away the concrete-drive mass *distribution*; replacing
  `c[q,d]` by a small histogram of caps (a second-moment or quantile cap)
  would directly attack the class-32 obstruction.  **Recommended research
  direction for `g=8`.**
- **Exact chamber arrangements**: same as Q7 — endorse.
- **Different enumerator/composition theorem**: no easy lift between `g`
  values exists; the `g=4` group neither contains nor is contained in the
  `g=8` ensemble in a way that transfers first moments.  The ladder is
  method-development only, as the team already treats it.

---

## 9. Q9 — the discriminating experiment

**Run the §7.5 BSP construction on the single worst cell first** — worst by
union *term* `U_σ + log2 n_σ` (see the consistency note in §7.5.5), which
may differ from the `+67198.7278` worst-*branch* cell.  This is roughly one
day of work, needs no new
certificate machinery (the outward evaluation path already exists), and is
decisive in both directions:

- if the fully refined leaves close — equivalently, if
  `max over the cell of E(a) = min_i F_i(a)` is below target — the `+67k`
  there is *provably* cover looseness, the same computation emits the exact
  cuts and owners that certify the cell, and the systematic pass (§10,
  step 2) is de-risked before any pipeline investment;
- if some fully refined leaf fails at a point `a*`: `a*` is a concrete
  profile where the finite atlas envelope cannot close — feed it to the
  existing pointwise tuning; if tuning also fails, `a*` is the first
  candidate for a genuinely difficult profile family, and only then is the
  expensive realizability machinery (product-kernel search, thermodynamic
  integration) worth running *at `a*` specifically*.  No subdivision
  strategy of any kind would have closed that cell, so the negative result
  saves the entire optimization effort.

The cell is small, the active witness set is ≤ tens, and the restricted
minimization diagram is a tiny arrangement, so the exact computation is
cheap.  This is strictly more informative than more refinement waves or more
sampled tuning.  Complement: the per-cell certified Jensen-gap bound (§8,
item 1) distinguishes "no mixture can ever close this cell" from "atlas too
small" without any LP run, and can pre-sort the failing-cell queue.

---

## 10. Q10 — prioritized path to `g=8` with stop/go criteria

1. **(Now) Close the strict-soundness items on `g=2`**: underflow guard
   (F3), the ensemble-uniformity lemma in the manuscript (F1), constant
   provenance (F6), spectrum MacWilliams check (F7).  Cheap; converts the
   `g=2` result into an unqualified theorem.  *Stop/go: all four are
   mechanical; no research risk.*
2. **(`g=4`, weeks) Per-cell dominance-BSP certificate on the full-support
   stratum** (§7.5): first the worst-cell experiment (§9), then a parallel
   sweep over the failing-cell queue ordered by residual contribution, with
   pure owner witnesses at leaf vertices, cell-level (or leaf-level) box
   lattice counts, and the existing outward hardening and ledger.  The
   frozen mesh and all closed cells stay byte-identical.  *Go if* the
   worst-cell experiment closes and the residual total drops by orders of
   magnitude on the first shallow (depth 3–5) pass; *fall back* to batched
   global anchor insertion + one retriangulation only for cells that
   exhaust the depth cap; *stop and pivot to new branches if* any fully
   refined leaf vertex fails pointwise after dedicated tuning — that is a
   proof-content gap, not a geometry gap.  The globally conforming
   chamber-cut variant (§7.2) remains a sound fallback but should not be
   the default, for the reasons in §7.5.5.
3. **(`g=4`) Assemble the 30-stratum ledger** with the cell-local union and
   the `g=2`-style aggregation.  Target: the analogue of the 62.7-bit
   cell-local result.  This is the template deliverable for `g=8`.
4. **(gate before `g=8` geometry)** (a) verify the `g=8`
   ensemble-uniformity lemma and implementation conformance; (b) run
   pointwise envelope transition scans on the dominant high-dimensional
   supports; (c) measure chamber growth on `g=4` and extrapolate `W⁴`
   against budget.  *Go to direct chambers only if* extrapolated cells
   `≲ 10⁷`; otherwise:
5. **(`g=8` research track) Witness compression / stronger local
   inequality**: replace the point-cap by a quantile/histogram cap (§8) to
   kill the middle-class obstruction, and search for a ≤2-parameter witness
   family per support whose envelope is analytically monotone in
   identifiable directions, certified by low-dimensional interval methods.
   *Stop/go: if neither reduces the effective `W` for the full-support
   `d=8` stratum below ~100, a direct `g=8` first-moment chamber theorem is
   out of reach with this architecture; consider proving `g=8` with a
   coarser but analytically summable outer (e.g., the total-spectrum
   envelope) at a reduced security target, or accept `g=4` for deployment.*
6. **(parallel, cheap) Second-moment sanity at the deployment parameters**:
   the first-moment margins at `g=2` (22.7 bits beyond target) suggest room;
   a small-`N` exact IOWE surrogate with packet-orbit denominators (already
   planned) guards against construction-level surprises that the first
   moment cannot see.

---

## 11. Missing artifacts / invariants (per the request to name them)

1. A numbered manuscript lemma for interleaver transitivity on weight-profile
   classes, plus an implementation conformance statement (F1).
2. An underflow-freedom invariant for every `*_outward` float DP (F3): a
   boolean reachability DP compared against float positivity, recorded in
   the witness report.
3. Constant-provenance table binding `K, D, N, INNER_BLOCKS,
   GRAPH_REPLACEMENTS`, band sizes 42/86, and `BLOCK_BITS` to the frozen
   construction document (F6).
4. MacWilliams/self-consistency check artifact for `EBCH128_64.wd` and the
   systematic split-slice table, with digests (F7).
5. A dedicated audit note for `packet_group_outer_profile.py` (the only
   sizable certificate-path module I did not fully re-derive) (F2).
6. For the recorded `g=2` result: a re-verification log tying the quoted
   SHA-256 values in `PROOF_STATUS.md` to a fresh verifier run (F5).

---

## 12. Register of labeled conclusions

| # | Conclusion | Label |
|---|---|---|
| 1 | `g=2` convexity/mixture/hybrid/union architecture and code | PROVED |
| 2 | `g=2` mesh conformity audit (2-chain argument, T-junctions) | PROVED |
| 3 | Pick counts, terminal enumeration, dedup, ownership | PROVED |
| 4 | Outward interval layer, Robbins bounds, normalization `Q(a)` | PROVED |
| 5 | Drive-stratified lemma + rank-one transport + Collatz composition | PROVED |
| 6 | `outward_transport_image` greedy dominates true optimum | PROVED |
| 7 | Ensemble uniformity on profile classes (justifies `Q(a)`) | PLAUSIBLE — needs manuscript lemma + implementation check |
| 8 | All mixture components bound the same `Z(a)` | PLAUSIBLE — semantic convention, needs closed documented family |
| 9 | Outward float DP underflow handling | UNSOUND as stated (immaterial in practice; cheap fix) |
| 10 | `g=4` support partition, shifted hulls, candidate lemma | PROVED |
| 11 | `g=4` verifier hull/facet/volume/conformity audits | PROVED |
| 12 | Kink-facet table in the g4 note | PLAUSIBLE (not in trusted base) |
| 13 | Producer/verifier separation (Qhull, LP, power iteration as discovery) | PROVED |
| 14 | Three-band BL polytope lemma algebra + rank-64 decomposition use | PROVED (modulo standard finite-field BL + hashed rank certificate) |
| 15 | Full outer moment algebra (`symmetric_s2`, band norms, graph terms) | PLAUSIBLE — module not fully re-derived |
| 16 | `+67k` residual owned by minimax-game/Jensen gap, not atlas truncation | DIAGNOSTIC (analysis; elementary and checkable) |
| 17 | Chamber-complex termination property (envelope closes ⇒ zero refinement) | PROVED (as a theorem statement; envelope closure itself is open) |
| 18 | Chamber feasibility at `g=4` (`O(W²)`), infeasibility of naive `g=8` (`O(W⁴)`) | PLAUSIBLE (worst-case bounds; measure actual growth) |
| 19 | Interval-Taylor box B&B inferior to vertex evaluation here | PROVED (split-max slack argument) |
| 20 | Point-cap histogram/quantile strengthening as the `g≥8` lever | PLAUSIBLE — research direction |
| 21 | Recorded `g=2` numerical interval and hashes | DIAGNOSTIC until independently re-run |
| 22 | EBCH spectrum tables | declared input; PLAUSIBLE pending MacWilliams check |
| 23 | Per-cell BSP bound `U_σ` validity (§7.5.3: leaf cover by construction, convex owner vertex max, inherited stratum eligibility, overcount-safe leaf duplicates) | PROVED (as a theorem statement) |
| 24 | Per-cell BSP dominates global conforming cuts and atlas/LP retuning (§7.5.5) | DIAGNOSTIC (engineering judgment on top of proved inequalities) |
| 25 | Residual mass concentration (LSE total exceeds its largest term by ≤ `log2 37750 ≈ 15.2` bits ⇒ the heaviest cells own the residual; the 53.7-bit total-vs-branch gap implies the worst terms carry ~38–54 bits of lattice-count factor, so order the queue by `U_σ + log2 n_σ`) | PROVED (elementary LSE bound + the reported `+67252.4553` / `+67198.7278` figures) |
| 26 | Worst-cell BSP experiment is decisive in both directions (§9) | PROVED (dichotomy follows from #17 and #23) |
