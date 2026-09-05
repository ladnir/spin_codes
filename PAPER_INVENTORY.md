# Paper Inventory

Generated: 2026-06-12

Purpose: classify the current manuscript before restructuring. This inventory
is a cleanup guide, not a proof artifact.

Implementation update: the split has been applied.  The current compiled
main spine is `intro.tex`, `prelim.tex`, `framework.tex`, `innerAcc.tex`,
`randomDenseConstruction.tex`, and `localCodeStructured.tex`.  The random
dense file is now a real compiled parent section; it inputs `outerDense.tex`,
`innerDenseScalar.tex`, and `integrationDense.tex` as subsections.  The
structured file is also a real compiled parent section; it inputs
`localCodeOuter.tex`, `localCodeInner.tex`, `localCodeCertificate.tex`, and
`localCodeProjections.tex` as subsections.
The original `innerDense.tex` and
`integration.tex` remain preserved sources.  The full-split certificate
development ledger now lives in
`explorations/fullsplit_certificate_development.tex`, with the compact main
certificate in `localCodeCertificate.tex`.
A follow-up cleanup pass removed the accumulator TODO subsection from the
compiled main path and preserved it in
`explorations/accumulator_warmup_notes.tex`; it also renamed the dense
integration and structured local-code sections to match the intended paper flow.
A second spine-tightening pass moved scalar dense-inner positional refinements
out of `innerDenseScalar.tex` and into
`explorations/scalar_dense_inner_refinements.tex`, leaving the compiled dense
inner section focused on the theorem-facing envelopes consumed by
`integrationDense.tex`.
The random sliding outer section has also been trimmed so fixed-tap tiny-weight
diagnostics live in `explorations/random_sliding_outer_low_weight_notes.tex`;
the compiled `outerDense.tex` now keeps the span generating law, geometric
spectrum envelope, and linear-weight exponent used by `integrationDense.tex`.
The structured local-code certificate section has been tightened so the main
paper states the RM/EBCH finite certificate and compact row-cover facts, while
the full post-prefix row table remains in
`explorations/fullsplit_certificate_development.tex` and the audit packet.
The row-level RM-prefix/post-prefix wrapper statements have now been collapsed
into a manifest-backed ledger summary in `localCodeCertificate.tex`; the former
formal wrappers remain preserved in the full-split exploration note.
The compiled main spine has also been language-audited so cleanup/provenance
phrases no longer appear as manuscript prose; remaining projection language is
intentional status marking for non-theorem BCH spectrum-model rows.
A main-spine audit removed compiled `\stan{...}` comments, stale
working-save-point wording, and the last layout warning in the framework
first-moment proof.  The accumulator section now includes the promised combined
warmup theorem: the random sliding dense outer plus accumulator inner has
linear distance for sufficiently small constants.  Remaining uses of
"diagnostic", "projection", and "heuristic" are intentional status boundaries
for non-theorem projection rows.
The scalar dense inner section has begun its proof-facing compression pass:
the run-tail exponent proof is now shorter in the compiled spine, the
monotonicity sign in that proof has been corrected, and the expanded Stirling
calculation is preserved in `explorations/scalar_dense_inner_refinements.tex`.
Unused fixed-tap finite-diagnostic bookkeeping lemmas have also been moved out
of `innerDenseScalar.tex`; exact final-survivor averaging, episode budgeting,
multi-candidate termination, and isolated-pair candidate counts now live in the
same scalar dense-inner refinement note.
First-start-aware no-OFF tails, zero-gap survival, and the single-run base
factor were also demoted to that note; they remain available for finite
diagnostics but no longer interrupt the compiled dense theorem spine.
The unused `OffAfter(T)` declaration has been removed from the compiled scalar
dense section, and the original `innerDense.tex` and `integration.tex` files
now have visible archival/provenance headers.
The dense integration section has also been compressed: baseline admissible
region and memory-scaling corollaries now live in
`explorations/scalar_dense_finite_diagnostics.tex`, leaving the compiled
integration spine focused on the qualitative theorem, linear-weight criterion,
and concrete `0.109` theorem.
The local direct first-moment restatement and integration-orientation
paragraphs have also been moved to
`explorations/scalar_dense_finite_diagnostics.tex`; the compiled integration
section now cites the framework theorem directly.
The structured local-code projection material has also been demoted: the
compiled paper now keeps only a compact BCH projection status note, while the
detailed BCH256/BCH512 comparison table and low-weight inflation stress tests
live in `explorations/outer_spectrum_comparison.tex`.
The random sliding outer section has been renamed to match the main-paper flow,
and duplicated interpretation/summary prose now lives in
`explorations/random_sliding_outer_low_weight_notes.tex`.
The scalar dense inner section was tightened again by moving the unused
polynomial-prefactor geometric envelope, the ON-to-OFF explanatory remark, and
the parameter-use summary into
`explorations/scalar_dense_inner_refinements.tex`; the compiled source keeps
the constant-factor geometric envelope consumed by the dense+dense theorem.
Another scalar dense-inner polish pass removed remaining commentary paragraphs
about non-main refinements, linear OFF-budget directions, and envelope
interpretation from the compiled source; the useful interpretation text now
lives in `explorations/scalar_dense_inner_refinements.tex`.
The framework section has now been shortened as well: optional IOWE and
piecewise/single-envelope bookkeeping variants were moved to
`explorations/framework_flexible_corollaries.tex`, leaving the compiled source
focused on the first-moment theorem and reusable envelope corollary.
The dense integration proof now handles the tiny-window late-placement summand
directly with the outer generating function \(W_{\mathrm{out}}(\lambda)\), while
fixed-tap band-cluster diagnostics remain outside the compiled theorem path.
The concrete `0.109` theorem now explicitly chooses a valid tiny-window
\(\xi_{\mathrm{tiny}}=1/4\), keeping it separate from the linear-window
\(\xi=8\) numerical gap check.
Its low-weight window also now uses the global outer envelope \(A_h\le3^h\),
avoiding the loose \(n^{O(1)}\) prefactor in the low-linear outer bound.
The random dense section has also been physically wrapped: compiled prose now
has one parent random sliding dense construction section with outer, inner, and
integration subsections.

## Status Labels

- `MAIN_THEOREM`: theorem-facing material that should remain in the main
  theory paper.
- `MAIN_CERTIFICATE`: finite certificate material needed for the current
  RM/full-split result.
- `SUPPORTING_LEMMA`: proof machinery that can stay near a main result or move
  to appendix.
- `DIAGNOSTIC`: useful finite evidence or sanity checks, not theorem text.
- `HEURISTIC_PROJECTION`: BCH/random-like/spectrum-model projections.
- `SUPERSEDED`: older paths replaced by later evidence.
- `OPEN_PROOF_DEBT`: unresolved theorem target or missing rigorous envelope.

## Top-Level Files

| file | current role | recommendation |
| --- | --- | --- |
| `intro.tex` | Short motivation and paper framing. | Keep main, later rewrite after theorem/certificate split stabilizes. |
| `prelim.tex` | Basic notation and combinatorial preliminaries. | Keep main. |
| `framework.tex` | Modular first-moment serial-concatenation framework. | `MAIN_THEOREM`; keep main with optional variants demoted. |
| `innerAcc.tex` | Accumulator exact enumerator, contraction, and combined warmup theorem with the log-memory dense outer. | `MAIN_THEOREM` warmup; keep main immediately after the framework. |
| `randomDenseConstruction.tex` | Parent section for the random sliding dense construction. | `MAIN_THEOREM` spine; compiled by `main_permConv.tex`. |
| `outerDense.tex` | Random sliding dense outer construction and spectrum bounds. | `MAIN_THEOREM` plus `SUPPORTING_LEMMA`; compiled as a subsection of `randomDenseConstruction.tex`. |
| `innerDenseScalar.tex` | Extracted scalar dense inner theorem material. | `MAIN_THEOREM`; compiled as the random dense inner subsection. |
| `integrationDense.tex` | Extracted analytic dense+dense integration material. | `MAIN_THEOREM`; compiled as the dense+dense asymptotics subsection. |
| `localCodeStructured.tex` | Parent section for the local-code independent structured construction. | `MAIN_CERTIFICATE` spine; compiled by `main_permConv.tex`. |
| `localCodeOuter.tex` | Direct-sum local-code outer interface and finite certificate template. | `MAIN_CERTIFICATE` interface; compiled as a subsection of `localCodeStructured.tex`. |
| `localCodeInner.tex` | Full-split local-code recursive inner interface and one-step split law. | `MAIN_CERTIFICATE` interface; compiled as a subsection of `localCodeStructured.tex`. |
| `localCodeCertificate.tex` | Manifest-backed RM/EBCH finite certificate. | `MAIN_CERTIFICATE`; compiled as the proved instantiation subsection. |
| `localCodeProjections.tex` | Short BCH/local-code projection status note. | Brief `HEURISTIC_PROJECTION`; compiled as the upgrade-path subsection. |
| `innerDense.tex` | Original monolithic dense-inner workspace. | Preserve as provenance; active content has been split into compiled spine and exploration notes. |
| `integration.tex` | Original monolithic dense integration workspace. | Preserve as provenance; active content has been split into compiled spine and exploration notes. |
| `outerExpandAcc.tex` | Alternate expander/accumulator outer construction. | Preserve as a separate exploration/construction note; do not keep in the main paper spine for this cleanup pass. |
| `innerSparse.tex` | Sparse recursive inner with refresh. | Preserve as a separate exploration/construction note; do not keep in the main paper spine for this cleanup pass. |
| `AUDIT.md` | Dense+Dense audit packet and proof-state ledger. | Preserve as audit packet; use as source, not main-paper prose. |

## `framework.tex`

Classification:

- `MAIN_THEOREM`: SC framework theorem, first-moment reduction, envelope
  corollaries.
- `SUPPORTING_LEMMA`: optional IOWE viewpoint and flexible corollaries.

Destination:

- Keep in main paper before construction-specific sections.
- Later cleanup should shorten explanatory prose, but no demotion is needed.

## `outerDense.tex`

Classification:

- `MAIN_THEOREM`: systematic dense/banded outer construction, fixed-word parity
  law, global geometric spectrum envelope, linear-weight outer exponent.
- `SUPPORTING_LEMMA`: low-linear envelope and high-weight exponent corollaries.
- `DIAGNOSTIC`: any discussion that only exists to motivate old finite scalar
  scans should move with scalar diagnostics if it becomes lengthy.

Destination:

- Keep theorem statements and proof machinery in main or theorem appendix.
- In final cleanup, expose only the interface consumed by `integration.tex` in
  the main flow.

## `innerAcc.tex`

Classification:

- `MAIN_THEOREM`: accumulator exact enumerator and low-output contraction as a
  clean first instantiation of the framework.
- `SUPPORTING_LEMMA`: optional sharpenings and TODOs.

Destination:

- Keep in the main paper as a warmup after the framework and outer interface.
- Present it as "log-memory outer plus accumulator inner gives linear distance,
  but with a weak constant." This demonstrates the framework before the harder
  dense recursive inner.
- Do not position it as the performance construction.

## `integration.tex`

Classification by section cluster:

- `Direct first-moment reduction`: already handled by `framework.tex`; the
  duplicate local restatement has been demoted.
- `A baseline linear-distance theorem`: `MAIN_THEOREM`; keep.
- `A concrete 10.9-percent dense+dense theorem`: `MAIN_THEOREM`; keep and
  cross-reference `RESULTS_LOCK.md` verification note.
- `A linear-weight exponent criterion`: `SUPPORTING_LEMMA`; keep near theorem
  or move to appendix depending on final length.
- `What must be sharpened next`, `How the enumerator enters`, and
  `Dense+dense estimates (working save point)`: mostly `DIAGNOSTIC` and
  `OPEN_PROOF_DEBT`; move to scalar dense finite diagnostics note.
- `Sigma-margin scan for the terminated-tail model`: `DIAGNOSTIC`; preserve as
  comparison evidence, not main theorem.
- `A block-outer finite-n upgrade` and old RM checkpoint in `integration.tex`:
  mixed `DIAGNOSTIC`, `SUPERSEDED`, and historical bridge. The current main RM
  certificate now lives more cleanly in `innerDense.tex`; demote this older
  version unless it is needed as background.
- BCH spectra and affine-coset material: `HEURISTIC_PROJECTION` and
  `OPEN_PROOF_DEBT`; move to `explorations/bch_outer_spectra.tex`.
- Exact isolated statistics and residual fit material: `DIAGNOSTIC`; move to
  scalar dense diagnostics or a proof-target note.

Destination:

- Main paper keeps the analytic dense+dense theorem line and a brief note that
  finite scalar diagnostics motivated later construction changes.
- Long finite scalar scans, BCH spectra exploration, and old RM-block history
  move to exploration notes.

## `innerDense.tex`

Classification by section cluster:

- Scalar dense recursive inner definition and ON/OFF lemmas:
  `MAIN_THEOREM` / `SUPPORTING_LEMMA`; keep for analytic dense+dense line.
- General three-term tail bound and dense inner interface theorem:
  `MAIN_THEOREM`; keep.
- Forced-termination contraction, isolated-slice, early-start, and pair
  statistics:
  `DIAGNOSTIC` / `OPEN_PROOF_DEBT`; preserve in scalar dense finite
  diagnostics, not main theorem text.
- `Prospective block-recursive BCH inner`:
  `DIAGNOSTIC`, `HEURISTIC_PROJECTION`, and `SUPERSEDED` relative to the
  full-split path; move to `explorations/block_recursive_bch_inner.tex`.
- Full-codeword random-split variant:
  `MAIN_CERTIFICATE`; keep as the current finite construction definition.
- Certified full-split termination atoms, first-active placement coverage,
  exact RM support split, small-prefix rows, post-prefix rows, fixed-pole
  wrappers, and current finite checkpoint:
  `MAIN_CERTIFICATE`; keep, but later compress into a clean certificate
  section plus appendix proof details.
- `Spectrum-model outer projections`:
  `HEURISTIC_PROJECTION`; move to an outer comparison or projection note, with
  only a short main-paper mention.
- Dense-inner TODOs:
  `OPEN_PROOF_DEBT`; move to inventory/proof-debt notes.

Destination:

- Split this file first during Phase 3. Keep scalar theorem and full-split
  certificate in compiled paper; move research-log diagnostics out.

## `outerExpandAcc.tex` And `innerSparse.tex`

Classification:

- `DIAGNOSTIC` / `OPEN_PROOF_DEBT`: useful alternate construction lines, not
  part of the cleaned main paper spine.
- Some local lemmas may later become `SUPPORTING_LEMMA` in their own notes.

Destination:

- Preserve intact during the first cleanup pass.
- Move out of the main compiled flow into separate construction notes or
  appendices after the dense+dense split is stable.
- Do not delete or rewrite them as part of the dense+dense consolidation.

## `AUDIT.md`

Classification:

- `MAIN_CERTIFICATE`: authoritative list of checkable full-split ledger pieces.
- `OPEN_PROOF_DEBT`: explicit proof debt and "do not trust yet" lists.
- `DIAGNOSTIC`: command histories and implementation audit details.

Destination:

- Preserve intact as the audit packet.
- Use it to maintain `RESULTS_LOCK.md` and drive future external audits.
- Do not copy long command transcripts into the main paper.

## Proposed Exploration Notes For Phase 3

Do not move prose until this inventory is reviewed. Proposed destinations:

- `explorations/scalar_dense_finite_diagnostics.tex`
  - old fixed-tap dense finite scans;
  - isolated-slice and early-pair diagnostics;
  - sigma-margin scans and scalar proof-target notes.
- `explorations/scalar_dense_inner_refinements.tex`
  - exact gap-composition laws, isolated-slice order statistics, paired-cover
    bounds, and linear OFF-budget refinements moved out of the compiled spine.
  - polynomial-prefactor geometric envelope and ON-to-OFF bookkeeping remark
    demoted from the compiled scalar dense-inner source.
- `explorations/random_sliding_outer_low_weight_notes.tex`
  - fixed-tap no-weight-one and band-cluster low-weight parity law for scalar
    dense diagnostics.
  - demoted dense-outer interpretation and theorem-interface summary notes.
- `explorations/accumulator_warmup_notes.tex`
  - optional accumulator sharpness improvements and non-main TODOs.
- `explorations/fullsplit_certificate_development.tex`
  - development history of the full-split certificate;
  - row-family derivation notes too detailed for main paper.
- `explorations/block_recursive_bch_inner.tex`
  - block-recursive systematic/permute BCH inner attempts;
  - column-profile and exact-turnoff diagnostics.
- `explorations/bch_outer_spectra.tex`
  - BCH spectrum search, MacWilliams/affine-coset program, public tables.
- `explorations/outer_spectrum_comparison.tex`
  - banded versus direct-sum block enumerator comparison;
  - random-block, RM, BCH-like, and sliding-structured outer sanity checks;
  - detailed BCH256/BCH512 projection table and low-weight inflation stress
    tests moved out of the compiled local-code section.
- `explorations/outer_expand_accumulator.tex`
  - preserved expander/accumulator outer material from `outerExpandAcc.tex`.
- `explorations/inner_sparse.tex`
  - preserved sparse recursive refresh material from `innerSparse.tex`.

## Main-Paper Target Shape After Review

Recommended compiled paper outline:

1. `intro.tex`, `prelim.tex`, `framework.tex`.
2. Accumulator warmup:
   - log-memory outer interface;
   - exact accumulator enumerator;
   - accumulator contraction;
   - simple linear-distance theorem with weak constant.
3. Random sliding dense construction:
   - concise outer interface;
   - dense inner interface;
   - dense+dense integration theorem;
   - concrete `0.109` asymptotic theorem.
4. Structured local-code construction:
   - local-code independent outer interface using `W_loc(z)^B`;
   - local-code independent full-split recursive inner interface;
   - finite certificate template;
   - proved RM outer plus EBCH inner instantiation at `N=2^21`,
     `delta=.09`;
   - compact result table.
5. BCH and better local-code instantiations:
   - low-weight spectrum bottleneck;
   - BCH spectrum requirement;
   - BCH256/BCH512 projections, explicitly heuristic;
   - sliding versus direct-sum block outer comparison.

## Phase 3 Guardrails

- Do not delete exploratory material on the first restructure pass.
- Keep labels stable where possible.
- If labels move, add a reference map.
- Do not present BCH/random-like rows as theorem text.
- Do not promote old scalar finite diagnostics to the current main finite
  theorem.
- Keep generated scratch files and `__pycache__` out of the paper inventory
  unless they are canonical artifacts listed in `RESULTS_LOCK.md` or
  `AUDIT.md`.
