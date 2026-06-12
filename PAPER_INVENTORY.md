# Paper Inventory

Generated: 2026-06-12

Purpose: classify the current manuscript before restructuring. This inventory
is a cleanup guide, not a proof artifact.

Implementation update: the first split has been applied.  The current compiled
main spine is `intro.tex`, `prelim.tex`, `framework.tex`, `outerDense.tex`,
`innerAcc.tex`, `innerDenseScalar.tex`, `integrationDense.tex`, and
`localCodeStructured.tex`.  The original `innerDense.tex` and `integration.tex`
remain preserved sources.  The full-split certificate development ledger now
lives in `explorations/fullsplit_certificate_development.tex`, with the compact
main certificate in `localCodeStructured.tex`.
A follow-up cleanup pass removed the accumulator TODO subsection from the
compiled main path and preserved it in
`explorations/accumulator_warmup_notes.tex`; it also renamed the dense
integration and structured local-code sections to match the intended paper flow.

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
| `framework.tex` | Modular first-moment serial-concatenation framework. | `MAIN_THEOREM`; keep main, possibly tighten prose. |
| `outerDense.tex` | Dense/banded outer construction and spectrum bounds. | `MAIN_THEOREM` plus `SUPPORTING_LEMMA`; keep main or later appendix split. |
| `innerDenseScalar.tex` | Extracted scalar dense inner theorem material. | `MAIN_THEOREM`; current compiled source for the random dense inner interface. |
| `integrationDense.tex` | Extracted analytic dense+dense integration material. | `MAIN_THEOREM`; current compiled source for random sliding dense construction and asymptotics. |
| `localCodeStructured.tex` | Compact local-code interface, RM/EBCH certificate, and BCH projections. | `MAIN_CERTIFICATE` plus `HEURISTIC_PROJECTION`; current compiled source for the structured construction line. |
| `innerDense.tex` | Original monolithic dense-inner workspace. | Preserve as provenance; active content has been split into compiled spine and exploration notes. |
| `integration.tex` | Original monolithic dense integration workspace. | Preserve as provenance; active content has been split into compiled spine and exploration notes. |
| `outerExpandAcc.tex` | Alternate expander/accumulator outer construction. | Preserve as a separate exploration/construction note; do not keep in the main paper spine for this cleanup pass. |
| `innerAcc.tex` | Accumulator inner instantiation with exact enumerator. | Promote as a warmup theorem/example after the framework. |
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

- `Direct first-moment reduction`: `MAIN_THEOREM`; keep.
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
  - random-block, RM, BCH-like, and sliding-structured outer sanity checks.
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
