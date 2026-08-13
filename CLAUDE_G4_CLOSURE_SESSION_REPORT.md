# Second independent session: `g=4` closure audit

Auditor: Claude (Fable 5), 2026-08-12.  Follows `CLAUDE_SECOND_OPINION_REPORT.md`.
This session inspected the BSP producer/verifier code line-by-line, replayed the
stored `s1fr032620` artifacts numerically, and ran three deliberately bounded
local experiments (each seconds-to-minutes; no Peach runs, nothing concurrent).
Labels: **PROVED** / **DIAGNOSTIC** / **PLAUSIBLE** / **SPECULATIVE**, with the
promotion requirement stated per item.

---

## 1. Headline verdicts

1. **No fatal flaw.  The BSP architecture is sound as implemented** (§2).  The
   single-cell verifier independently replays geometry, vertices, owners, and
   outward values; producer discovery cannot contaminate it.  Three verifier
   invariants are missing before a global certificate can be assembled (§2.3).
2. **The stubborn profile closes.  The "300.25-bit miss" at `s1fr032620` is
   an evaluator artifact, not a proof-content gap.**  At the already-stored
   witness `(pole, fugacities)`, converging the Collatz power iteration
   (160+ iterations instead of 40) improves the inner bound by ~318 bits and
   the combined value crosses the target with **+22.48 bits of margin**
   (measured, binary64; §3).  Mechanism: the diagnostic λ̂ was 0.010284
   bits/block at 40 iterations vs 0.000124 converged, and `×32768` blocks
   amplifies that to ~330 bits.  **DIAGNOSTIC**, promotable to PROVED by one
   single-cell outward rerun with a converged (or frozen) Collatz vector.
3. **The failure census is systematically inflated by the same artifact.**
   Measured on the 24 smallest-gap census rows: **19 of 24 close by
   convergence alone** — re-evaluating their stored best witness at 320
   iterations, touching nothing else (gains 0–450 bits per witness; §3.4).
   The tuner's internal objective used `witness_iterations=16`, so every
   Powell run in the pipeline optimized an objective whose noise (hundreds of
   bits) exceeds the margins being chased.  This also explains the observed
   "regression" to 346.62 bits — incomparable evaluations, not a landscape
   fact (§4).
4. **The combined bound is exactly separable**: at a fixed profile,
   `combined = outer(λ, β) + inner(f, z) − normalization(a)`, with disjoint
   parameter sets.  Verified numerically: the exact-graph upgrade changed the
   combined value by exactly its outer delta, 122.2010 bits (§3.2).  There is
   therefore **no outer/inner saddle**; staged optimization is legitimate, and
   "joint outer+inner" search is unnecessary — optimize each side separately
   to its own optimum.  **PROVED** (algebra + numerical check).
5. The prompt's suspicion that the outer tilt was frozen is **architecturally
   correct but numerically immaterial here**: `outer_log2 =
   219479.65441553813` is bit-identical across all four attack stages — the
   tilt was never optimized at this profile — but my bounded convex
   optimization of the exact-graph total-spectrum outer recovered only
   **1.35 bits** (the inherited tilt from the nearby seed profile was already
   near-optimal; §3.3).  The linear-BL exact-graph family is ~182k bits worse
   at this profile and is not a candidate.
6. **Answer to the narrow question**: the most principled and economical
   elimination of the bottleneck is not a new inequality, a new outer family,
   or more Powell depth.  It is (a) converged/adaptive Collatz evaluation in
   every scoring path, (b) a converged re-hardening + re-ranking sweep of the
   atlas and census, and (c) freezing the Collatz vector in witness artifacts
   so certified margins stop depending on an iteration budget (§6).  Of the
   prompt's four hypotheses, the evidence assigns the current bottleneck to
   **(2) optimization/evaluation failure within an adequate witness family**
   — with (3) a secondary contributor and no observed instance of (1) or (4).

---

## 2. Soundness findings: BSP architecture vs implementation (Q1)

All findings below are from direct code reading of
`scripts/certify_packet_group_g4_cell_bsp.py` (verifier, 415 lines) and
`scripts/probe_packet_group_g4_cell_bsp.py` (producer, 484 lines).

### 2.1 Confirmed (PROVED, code-as-math)

- **Per-cell nonconformity is harmless.**  The verifier audits the *entire*
  source stratum mesh first (`audit_source_stratum`,
  certify_packet_group_g4_cell_bsp.py:157-205, calling
  `g4.verify_support_stratum_mesh` on all covered cells), then treats the BSP
  as interior decoration of one already-certified root simplex.  No BSP facet
  is ever matched against a neighboring cell, so hanging facets cannot arise
  by construction.
- **The split tree covers each parent exactly.**  Children are the parent
  intersected with the two *closed* half-spaces of the same hyperplane
  (`affine_constraint`, :126-132; `visit`, :227-258).  Closed∪closed ⊇ parent
  with overlap only on the cut plane — an overcount-safe boundary, consistent
  with the cell-local ledger.  The chart substitution `a4 = M − Σa_{0..3}`
  (:128-129) is algebraically correct, and child nonnegativity of `a4` is
  inherited from the root simplex constraints.
- **Leaf vertex enumeration is exhaustive over exact rationals.**
  `enumerate_vertices` (:90-104) solves every `C(n,4)` active constraint set
  with exact `Fraction` elimination and filters by all constraints; every
  vertex of a 4-polytope arises from some 4-subset of active facets, so no
  vertex can be missed.  Degeneracy only produces duplicates, removed by the
  set.  The producer's vertex digests (:244-245) are cross-checked but not
  load-bearing — the verifier re-enumerates.
- **One fixed witness per leaf, held over the whole leaf.**  Each leaf has a
  single `owner_witness` (:239-241) evaluated at every enumerated vertex
  (:334-361); there is no vertex-varying selection anywhere in the replay.
  Convexity of `F` (unchanged Lemma 1) then covers the leaf.
- **Support eligibility cannot leak.**  Every used owner is checked finite on
  every coordinate of the stratum support (`require_support_eligible`,
  :325-327) — stratum-level, stricter than leaf-level, safe.
- **Boundaries are only overcounted.**  Shared cut planes (closed children)
  and shared source-cell facets are counted at least once each; the union
  term is `max(leaf vertex values) + log2(root cell lattice count)` (:373-375),
  with the root count independently recomputed and matched against both the
  source ledger and the BSP artifact (:300-304).
- **No discovery sample can stand in for a certified region.**  The verifier
  ignores the producer's diagnostic leaf statuses and scores; it re-derives
  the partition from the root constraints and re-evaluates every leaf vertex
  outward.  A tree containing `pointwise_atlas_failure` leaves fails the
  outward gate loudly (:410-411).  The producer/verifier do share
  `solve_linear`/`enumerate_vertices` *implementations by duplication* (not
  import), which is acceptable; the exact-rational algebra was checked
  directly rather than inferred from agreement.
- The strict-cut requirement (:253-254) plus full-dimensionality check
  (:235-236) exclude degenerate or empty children.
- The previous report's underflow finding (F3) is fixed and verified:
  `packet_group_drive_stratified.py:104-108` raises on structurally positive
  products that underflow to zero (fails closed).

### 2.2 One scope note

`certify_packet_group_g4_cell_bsp.py` is honest about certifying **one cell
only** (:377-378).  Nothing yet prevents the *diagnostic* batch totals
(`+17712.39`, `+12078.98`) from being mistaken for certified quantities in
prose; they are binary64 replays of frozen constants.  Keep the current
labeling discipline.

### 2.3 Missing verifier invariants (named, per request)

1. **Global assembly verifier (blocking for the theorem).**  There is no
   artifact schema or checker that assembles the final union: every cell of
   all 31 strata appearing **exactly once**, tagged `mixture` or `bsp`, with
   the per-cell outward terms and one final outward log-sum-exp against
   `−40`.  Required checks: cell-id set equality with the audited mesh per
   stratum, no duplicates, digest binding of each per-cell artifact, and
   re-evaluation (not trust) of every term.
2. **Frozen Collatz vector in witness artifacts (blocking for thin margins).**
   The hardener re-*proposes* its Collatz vector with `--iterations`
   (default 40).  Any positive vector yields a valid bound, so this is sound
   — but the certified constant then depends on an unrecorded iteration
   budget, and §3 shows 40 iterations can sit 300+ bits above the converged
   bound.  Witness rows should store the frozen vector `v` (exact dyadics)
   and the verifier should check *that* vector's Collatz inequality, making
   certificates deterministic and iteration-independent.  Until then, at
   minimum record and enforce a per-witness minimum iteration count.
3. **Outward implementations of the exact-graph outers (blocking).**  The
   graph-upgraded witnesses use `exact_graph_total_spectrum` /
   `exact_graph_linear_bl`, which are binary64 diagnostics by their own
   docstring (probe_packet_group_exact_graph_linear_bl.py:12-14).  The
   outward hardener (`certify_packet_group_triangle_ledger._parameters`)
   supports only `linear_bl`, `total_weight`, and
   `exact_graph_puncture_total_weight`.  No graph-upgraded witness can
   currently enter a certificate.  The outward port needs the stated
   premises as a lemma: 128 holes in distinct band-zero groups (layout),
   conditional independence of hole-group factors given the graph word, the
   digest-bound graph24 spectrum, and the same α/β extraction argument as
   the existing total-weight lemma.  **PLAUSIBLE → PROVED** requires exactly
   this port plus the layout lemma reference.

---

## 3. Forensic decomposition of `s1fr032620` (Q2)

Profile `[444306, 21517, 37620, 17842, 3003]`; uniform per-profile target
`−111.41506501642425`; normalization `633526.1426` (common to all rows).
All numbers below are measured from the stored artifacts or reproduced
locally; my 40-iteration re-evaluation reproduces the stored
`inner_probability_log2 = −219168.6159…` exactly, confirming provenance.

### 3.1 The attack sequence, decomposed

| stage | outer branch | outer_log2 | inner_log2 | combined | miss |
|---|---|---:|---:|---:|---:|
| coordinate tuning | total_spectrum | +219479.6544 | −212384.0023 | +7095.65 | 7207.07 |
| cross-profile Powell | total_spectrum | +219479.6544 | −218844.9866 | +634.67 | 746.08 |
| two continuations | total_spectrum | +219479.6544 | −219168.6159 | +311.04 | 422.45 |
| exact-graph upgrade | exact_graph_total_spectrum | +219357.4534 | −219168.6159 | +188.84 | 300.25 |
| **this audit: converged inner** | exact_graph_total_spectrum | +219357.4534 | **−219491.3437** | **−133.89** | **closes, +22.48** |

Two structural facts fall out of the table:

- **The outer tilt never moved.**  `outer_log2` is bit-identical
  (`219479.65441553813`) across the first three stages; every one of the
  ~6,900 bits of progress came from the inner.  The joint Powell optimizes
  inner parameters only, with `outer_charge`/`outer_constant_log2` copied
  from the target row (probe_packet_group_joint_inner_opt.py:118-123), and
  the exact-graph upgrader re-evaluates the new outer at the *inherited*
  tilt (upgrade_packet_group_exact_graph_outer.py:35-56).
- **Exact separability.**  `Δcombined = Δouter = 122.2010` bits between rows
  3 and 4 with the inner untouched.  Outer and inner have disjoint parameter
  sets and interact only through the fixed profile.  **PROVED.**

### 3.2 Where the 300.25 bits actually lived

- **Outer relaxation slack: ≈ 1.35 bits.**  I optimized the exact-graph
  total-spectrum outer over its six parameters (five log-variables with
  anchor gauge, log β) at this profile: Powell from the frozen tilt converged
  at `+219356.1049` (recovered 1.348).  The function is convex in these
  parameters (log-sum-exp moments, max-over-`w` extraction, affine charge),
  so this is the family optimum; the large *symmetric* finite-difference
  slopes at the seed (−4k to −57k bits/unit) with near-zero one-sided
  improvement are the signature of a kinked (max-selection) valley floor,
  not a plateau artifact.  **DIAGNOSTIC** (binary64), promotable by an
  outward evaluation at the optimized tilt.
- **Alternative outer family: not competitive.**  `exact_graph_linear_bl` at
  the same tilt: best scanned value `+401075` (band1 = 0.6) — ~181,718 bits
  worse.  Consistent with the known total-spectrum dominance on zero-heavy
  profiles.
- **Inner evaluator non-convergence: ≈ 318–323 bits — the entire remaining
  miss.**  At the stored `(pole = 0.57338…, fugacities)`:

  ```text
  iterations   inner_log2        combined      margin      λ̂ (log2/block)
  40           −219168.6159      +188.84       −300.25     0.010284
  80           −219466.3460      −108.89       −2.52       0.001094
  160          −219491.3437      −133.89       +22.48      0.000135
  320/640      −219487.04        −129.58       +18.17      0.000124
  ```

  Every row is a *valid* Collatz upper bound (any positive test vector
  certifies); the 160-iteration vector happens to certify the most and may
  legitimately be frozen.  `λ̂ × 32768` blocks is the amplifier: 0.010284
  bits/block ⇒ 337 bits of avoidable slack at 40 iterations.
- **Graph/puncture correction: exactly 122.2010 bits** (separable, above).
- **Lattice-count charge:** the cell count enters after the pointwise bound;
  with the profile closing at +22 bits pointwise, the cell term needs the
  BSP leaf structure re-replayed (below), not more pointwise work.
- **Rounding gap (rational leaf vertex vs rounded integer profile): not yet
  measured.**  The census row's exact rational vertex is within ~0.5 of the
  rounded profile per coordinate; local affine slopes are tens of bits per
  unit, so the vertex value may differ by tens of bits — material against a
  +22-bit margin.  This is the one measurement I could not complete
  (requires re-evaluating the fractional-normalization path at the stored
  rational vertex) and it belongs in the primary experiment.

**Answer to Q2's key question:** "300 bits short" was **neither** the local
proof-content gap **nor** a meaningful rounded-profile diagnostic — it was
evaluator noise.  The actual local status: closes by ≈ +18–22 bits
(binary64) at the rounded profile with existing machinery; thin, so the
leaf-vertex and outward confirmations are mandatory before celebrating.

### 3.3 Optimization audit (Q3)

- **Separability makes the staged design correct in principle** — there is
  no joint saddle to miss (headline 4).  The suspected "saddle-point
  mismatch" from freezing the outer under the exact-graph upgrade is real in
  the code but cost only 1.35 bits here, because the tilt inherited from the
  nearby seed profile `[449873, 20797, 28562, 22088, 2968]` was almost
  optimal.  This will not always be true: **the outer tilt should be
  (convex-)optimized per target profile as a standard step** — it is a
  six-parameter convex problem at ~0.01 s/evaluation, essentially free.
- **The false plateau has a specific, fixable cause**: the Powell objective
  evaluates `inner_probability(…, iterations=16)`
  (probe_packet_group_joint_inner_opt.py:86-87) and final scoring uses 40.
  Measured objective noise from iteration starvation: 0–450 bits per witness
  (§3.4) — larger than every gap the optimizer was chasing in the endgame.
  Powell on an objective with state-dependent bias of that size produces
  exactly the observed pathologies: stalls, seed-dependent "regressions"
  (the 346.62 row is a different seed scored at a different effective
  convergence — not evidence about the landscape), and unrepeatable
  continuation gains.
- **Falsifiable test delivered and run** (per the request to prefer these
  over "more seeds"): re-evaluate the incumbent at converged iterations; if
  the miss survives, tune; if not, the plateau was noise.  It was noise.
- **Gauge and parameterization**: anchor-normalized fugacities are a genuine
  gauge fix (the combined bound is invariant under `f → c·f` since
  `16 × GROUPS = M`); pole bounds `(0.01, 0.95)` are inactive at the
  optimum; no other gauge freedom found.  The outer's max-over-`w` and the
  transport greedy introduce nonsmoothness — harmless for correctness, but
  prefer multistart Nelder-Mead/subgradient over Powell near kinks.
- **Recommended objective change**: adaptive escalation — double iterations
  until `λ̂_log2` changes by `< 1e-6`/block, then score.  At the measured
  0.5–1 s per converged evaluation this is affordable everywhere it matters.

### 3.4 Census-wide blast radius (measured)

`out/g4_cell_bsp_top1024_exact_graph_round1_failures.json`: 722 residual
profiles across 276 cells; gap distribution: max 7280.5, median 1033.9,
min 9.8; 56 rows ≤ 100 bits, 213 rows ≤ 500 bits.  Re-evaluating the stored
best witness of the 24 smallest-gap rows at 320 iterations (nothing else
changed): **19/24 close outright**; convergence gains ranged 0–450.7 bits
(a few witnesses were already converged, gain ≈ 0).  Extrapolation
(**DIAGNOSTIC**): a large fraction of the ≤ 500-bit tail closes for free, and
every larger gap shrinks by a witness-dependent amount before any retuning.
The frozen atlas constants themselves carry the same bias, so a converged
re-hardening also lowers every BSP leaf bound and both frontier ceilings.

---

## 4. Challenge to the outer families (Q4)

- **No missing outer inequality is implicated by the current frontier.**
  The evidence: the leading obstacle closed inside the existing family;
  19/24 sampled census failures close without touching the outer; the
  exact-graph total-spectrum family was at its convex optimum at the one
  profile examined in depth.
- Ranked answers to the specific candidates, with the lemma each requires:
  1. **Outward port of `exact_graph_total_spectrum`** — not a new bound but
     the certificate-path requirement (§2.3.3).  Do this first regardless.
  2. **Shaped/band-split spectrum with per-band poles** (the prompt's
     "graph conditioning + shaped spectrum") — legitimate upper bound by the
     same α/β extraction applied per band with the exact 42/86 band split
     spectra; moderate expected gain; build **only when** some profile
     reaches the redesign bar of §6.  **PLAUSIBLE.**
  3. **Combining two outer bounds before Cauchy extraction** (Hölder across
     group decompositions) — mathematically legitimate
     (`E[Π X] ≤ Π E[X^{1/θ}]^{θ}` with exponents summing appropriately), but
     it interpolates *between* families and cannot beat the better endpoint
     by much unless their tight regimes differ within one profile; no
     current need.  **PLAUSIBLE / low priority.**
  4. **Realizability exclusion** — the natural parity test fails: ordinary
     blocks force even weight, but the 128 punctures and graph replacement
     bits break the invariant, so total physical weight parity excludes
     nothing.  A support/multiplicity exclusion would need genuinely new
     structure; not worth it while obstacles keep dissolving.
     **SPECULATIVE.**
  5. **"Does the exact BCH spectrum already contain the needed
     information?"** — yes in the sense that (2) is the systematic way to
     expose it; the single-β total-spectrum charge is the only current
     compression, and it was not the binding constraint at `s1fr032620`.

---

## 5. Ranked proof refinements

1. Converged/adaptive Collatz evaluation in every scoring and tuning path
   (evaluator fix; hours of work; measured hundreds of bits everywhere).
2. Frozen-Collatz-vector witness schema + verifier consumption (§2.3.2);
   makes certified margins deterministic; required for thin-margin cells.
3. Outward port of the exact-graph outer branches with the layout/
   independence lemma (§2.3.3); unblocks every graph-upgraded witness.
4. Per-profile convex optimization of the outer tilt as a standard tuning
   step (free; removes the inherited-tilt risk that happened to be benign
   here).
5. Global assembly artifact + verifier (§2.3.1); mechanical but blocking.
6. Shaped band-split-spectrum outer — build only on the §6 redesign trigger.

---

## 6. The single best bounded experiment (Q5)

**Primary: the converged re-hardening and re-ranking sweep.**

- **Inputs**: `out/g4_cell_bsp_top1024_exact_graph_round1_failures.json`
  (722 rows), the frozen witness artifacts it references, the frozen atlas
  used by the top-1024 batch.
- **Code path**: add adaptive-iteration convergence (double until
  `Δλ̂_log2 < 1e-6`/block) to `inner_probability` call sites in
  `probe_packet_group_joint_inner_opt.py` and the census scorer; re-score
  every census row's referenced witness and the full-atlas envelope; then
  replay the top-1024 BSP against re-hardened constants
  (`probe_packet_group_g4_cell_bsp_batch.py`).
- **Cost**: measured 0.5–1.5 s per converged evaluation ⇒ the 722-row
  re-score is minutes on 16 workers; the atlas re-harden plus BSP replay is
  the existing batch runtime.  Bounded; no new geometry; single run.
- **Success**: census failure count drops by roughly half or more (sampled
  rate 19/24 on the small-gap tail); `s1fr032620` and most small-gap cells
  close; both frontier numbers (`+17712.39` processed, `+12078.98` ceiling)
  fall by thousands of bits.
- **Negative result teaches**: the survivors are the *real* obstacle list —
  profiles whose converged inner and convex-optimal outer still miss.  That
  list (expected small) is the redesign queue for §4.2, and reaching it
  would be the first genuine evidence of hypothesis (1) from the prompt.
- **Certificate feed**: mostly diagnostic, except item: rerun
  `certify_packet_group_g4_cell_bsp.py` on `s1fr032620` with
  `--iterations 320` (or the frozen-vector schema once implemented) — that
  output is outward and feeds the certificate directly.
- **Fallback A**: converged-objective Powell continuation on the surviving
  profiles only (bounded, per-profile).
- **Fallback B**: prototype the outward exact-graph total-spectrum port
  (needed for the certificate anyway, so no work is wasted).

Also fold in the one unfinished measurement from §3.2: evaluate the
re-hardened winning witness at the exact *rational* failing vertex of
`s1fr032620`'s leaf `n000001`, not only at the rounded integer profile, and
report the rounding gap explicitly.

---

## 7. Stop/go for `g=4` and the path to the outward certificate (Q6)

**Go.**  Positive evidence beyond "no counterexample found":

- Every closed obstacle in the program's history closed with large margins
  once evaluated correctly; the sequence of "hard" cases has been, in order:
  geometry slack (fixed by support stratification), atlas reuse (fixed by
  retuning), the adversarial hole tax (fixed by one true lemma), and now
  evaluator noise (fixed by convergence).  None was a wall; exactly one
  required new proof content.
- The frontier ceilings decline monotonically (`+67252 → +32688 → +22142 →
  +19154 → +12079`), and the measured convergence correction lowers both
  frontiers further before any new tuning.
- Direct measurement: the current worst cell's residual closes, and 19/24 of
  the smallest census failures close, with zero new proof content.

All of this is **DIAGNOSTIC** for the truth of the `g=4` statement — the
first moment is what it is; what has strengthened is the evidence that the
*bounding machinery* reaches it.  Distinct remaining obstacle classes after
`s1fr032620`: (i) whatever survives the converged sweep (expected: few,
possibly zero, genuinely hard profiles); (ii) the outward exact-graph port;
(iii) global assembly + the 30 low-dimensional strata (machinery exists);
(iv) thin-margin leaf closure where +20-bit pointwise margins meet rational
leaf vertices.  No evidence of a qualitative wall in any of them.

**Stop-tuning / redesign trigger** (the prompt's "at what observation"):
a profile that misses with (a) a converged, multistart-stable inner
(≥ 8 diverse seeds, adaptive iterations) **and** (b) a convex-optimized
outer within every implemented family.  That is the certified signature of
hypothesis (1) — a genuinely insufficient inequality — and the answer is
then the shaped-spectrum outer (§4.2) or a cap-histogram inner, not more
search.  **No profile has yet reached this bar** — `s1fr032620` failed (a).

**Ranked path to the outward end-to-end certificate, with gates:**

1. Converged sweep (§6).  Gate: census collapses; surviving list stable
   across two rounds.
2. Frozen-vector schema + outward exact-graph port (§2.3.2, §2.3.3).
   Gate: `s1fr032620` certified end-to-end by
   `certify_packet_group_g4_cell_bsp.py` with a graph outer.
3. Batch single-cell outward BSP certification for every cell whose mixture
   bound fails.  Gate: all top-1024 cells hold certified terms.
4. Global assembly artifact + verifier over all 31 strata (§2.3.1).
   Gate: every cell exactly once; outward LSE `≤ −40`.
5. Independent replay on a second machine; digests recorded in
   `PROOF_STATUS.md`.

---

## 8. Implications for `g=8`

- **Transfers directly**: separability of outer/inner (same algebra), the
  BSP-inside-frozen-mesh architecture, census-driven tuning, the assembly
  schema, and — most cheaply and most importantly — the evaluation
  discipline.  The λ̂ amplification is `×B = 32768` at `g=8` as well; adopting
  adaptive convergence and frozen vectors from day one avoids repeating this
  entire detour at nine profile classes and 510 supports.
- **Still expected to strain**: witness-count and chamber/BSP growth in
  dimension 8 (the `O(W⁴)` concern from the previous report stands), the
  per-profile tuning volume (722 census rows at `g=4` will multiply), and
  the middle-class point-cap obstruction documented for large `g`.  The
  converged sweep materially improves the outlook by showing much of the
  apparent tuning volume was noise — but dimension-8 geometry remains the
  honest risk, and the §7 redesign trigger should be carried over verbatim.

---

## 9. Missing artifacts or facts that limited this audit

1. The post-tuning top-1024 *replay* artifacts behind the quoted
   `+17712.39` / `+12078.98` figures were not present in `out/` (Peach-side);
   both numbers are quoted from `PROOF_STATUS.md` and could not be
   independently rechecked against JSON.
2. Witness rows do not record the iteration count that produced their
   `inner_constant_log2`; convergence status per stored witness is only
   observable by re-evaluation (as done here for 25 witnesses).
3. No outward implementation of `exact_graph_total_spectrum` /
   `exact_graph_linear_bl` exists, so no graph-upgraded value in this report
   can currently be promoted past DIAGNOSTIC.
4. The exact rational leaf-vertex value at `s1fr032620:n000001` under the
   re-hardened witness (the rounding-gap measurement) remains to be done.
5. No global assembly schema exists yet; the end-to-end union is currently a
   prose claim over per-cell artifacts.

---

## Appendix: labeled conclusions register

| # | Conclusion | Label | Promotion requirement |
|---|---|---|---|
| 1 | BSP verifier replay is sound (partition, vertices, owners, eligibility, overcount-only boundaries) | PROVED | — |
| 2 | Source-mesh audit is independent of BSP discovery | PROVED | — |
| 3 | Combined bound separates: outer(λ,β) + inner(f,z) − norm(a) | PROVED | — |
| 4 | Outer tilt never optimized at the stubborn profile (bit-identical across stages) | PROVED (artifact fact) | — |
| 5 | Exact-graph total-spectrum outer is family-optimal ±1.5 bits at this profile | DIAGNOSTIC | outward evaluation at optimized tilt |
| 6 | `s1fr032620` residual profile closes at stored witness, margin ≈ +22.5 bits (160 it.) / +18.2 (converged) | DIAGNOSTIC | single-cell outward rerun with converged/frozen vector |
| 7 | The 300.25-bit miss was Collatz-vector non-convergence (λ̂ 0.010284 → 0.000124 bits/block) | PROVED (measured mechanism) | — |
| 8 | 19/24 smallest census gaps close by convergence alone (gains 0–450.7 bits) | DIAGNOSTIC | full census re-score |
| 9 | Powell "regression" to 346.62 was incomparable evaluation, not landscape | PLAUSIBLE (strongly indicated) | rerun both rows at converged iterations |
| 10 | Linear-BL exact-graph family ~182k bits worse at this profile | DIAGNOSTIC | — (family comparison only) |
| 11 | No missing outer inequality implicated by current frontier | DIAGNOSTIC | converged sweep survivor list empty |
| 12 | Parity-based realizability exclusion fails (punctures/graph break even-weight) | PROVED | — |
| 13 | Frozen-vector witness schema removes iteration dependence from certificates | PROVED (design fact) | implement |
| 14 | `g=4` will close end-to-end | PLAUSIBLE | gates 1–4 of §7 |
