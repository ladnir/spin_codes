# Paper integration review

Completed 2026-09-07 in the authoritative GitHub workspace.
The manuscript is `main.tex`; the built output is
`../output/pdf/spin_codes_draft.pdf` (48 pages).
PDF SHA-256:
`6b648fd31e9f68964f49f60f4b70da06246861125f8588220d33d96618f06167`.

Artifact follow-up: added a repository README, a reader-facing artifact guide,
paper-to-code map, standard-library reproduction commands, evidence inventory
and ZIP packager, and Windows/Linux quick-check CI. The paper now links to
the repository, guides, figure generator, and encoder instructions. The quick
checks pass in an isolated compact tree; eight artifact tests pass. Missing
compact asymptotic inputs were restored and all 31 original manifest entries
now authenticate, with an explicit relocated selected-map input and pinned
LF checkout rules. The 788-entry BCH set is packaged locally in a 63.5 MiB
verified ZIP; see `../artifact/VALIDATION.md` for its checksum and release gaps.
No numerical proof search, benchmark, commit, push, or public release occurred.

The artifact PDF pass has no TeX warnings or unresolved references. The new
linked passages on pages 6, 21, 27, and 45 were inspected in rendered form.
The final parameter/theorem layout on pages 21--24 was re-rendered and checked;
the first plot follows its introduction and the finite theorem is not split
across pages. Main-branch links still require integration before publication.

Follow-up PDF comment 1: renamed Section 1.1 to “Three SPIN constructions.”
Rebuilt the PDF and checked the rendered heading and contents entry.

Follow-up parameter-plot revision: Section 7.1 now explains three controlled
slices before the selected construction: K/B at fixed (t,s)=(64,20), s/t
at fixed K=2^20, and K/s at fixed t=64. Figures 1--3 use the retained
BCH-64/128 nested-map study; Figure 4 remains the separate BCH-256 certified
curve. Weak bounds and first-moment obstructions are visible without being
presented as actual code failures. No unknown or weak full point is bridged
by a curve. The figures are native vector PGFPlots inputs.

`build_parameter_figures.py --check` validates 130 unique geometries (77
useful full bounds, 7 weak full bounds, 46 first-moment obstructions), exact
geometry coverage, duplicate consistency, and generated-source freshness.
Its source is the tracked, rounded `CURRENT_RESULTS_TABLES.md`, not an
unavailable raw numerical export; full values have six decimal places.
Q1 overlays use full margin plus the retained loss only where available.
This pass did not regenerate or independently replay the small-BCH bounds.
The original producer and input fingerprints remain recorded in the source
table and in the landscape study documentation.

The latest build has no undefined references/citations, box warnings, or
LaTeX/package warnings. Rendered pages 20--25 were inspected at reading
resolution, including all four figures, their captions, the construction,
theorem, and section transitions. Legends were moved away from plotted
curves, and a float barrier keeps the three diagnostic plots ahead of the
selected construction. The earlier whole-document visual review is recorded
below; this follow-up rechecked the changed main-text pages.

## Plan requirements and evidence

| Requirement | Implemented evidence |
|---|---|
| Shared architecture and clean asymptotic/finite framing | Section 6.8 transitions from realized-spectrum control for the growing outer to deterministic spectrum inequalities for the fixed BCH constituent. |
| Brief randomized-outer context | One sentence acknowledges the separately certified expander--accumulate example; no new subsection, formal theorem, timing, or historical inner discussion. Its internal source is `workstreams/finite_asymptotic_theory/SPARSE_EA_FREEZE.md`. |
| Preserve the generic construction and asymptotic result | The source diff preserves the construction and Theorem 6.3. The routing explanation now states the actual support constraint; the zero-extension exponent has its stray comma removed. |
| Dedicated finite section after structured SPIN | Section 7 defines the outer, selected maps, setup distribution, recurrence conventions, and Theorem 7.1. Theorem scope is the five listed message lengths at relative distance greater than 0.10 and setup failure below 2^-40. |
| Actual BCH constituent, not arbitrary dimensions | Appendix B.1 specifies the implemented polynomial generators and proves common spectra for the dimension-128 intermediate codes. |
| Mathematical finite appendix | Appendix B gives spectrum constraints and rational dual checking, the exact inner map, transfer matrices, Q1 coefficient domination, sparse/composition bounds, the shuffled Bernoulli lemma, two-tilt reduction, rectangle coverage, and exact final union. |
| Parameter effects and certified engineering curve | Figures 1--3 show matched/nested small-BCH slices in K, B, s, and t. Table 2 and Figure 4 use the five exact-ledger BCH-256 margins. Text explains the local slope, state plateau, and higher-occupancy constraint. Heuristic spectrum uncertainty is separate; the older (64,20) curve is not plotted as a matched estimate. |
| Replace obsolete extrapolation guidance | Section 8.3 replaces the old regression and persistence staircase with finite selection using complete bounds and matching selected maps. |
| Measured implementation | Section 9 defines the batched transposed workload and reports the six selected/reference timing cells, setup, memory, hardware, compiler, variation, and serial methodology. It makes no measurement claim at exponents 22 or 24. |
| Front matter and organization | Abstract and introduction include the finite theorem and measured encoder. Obsolete finite-result placeholders are removed. |
| Build and visual verification | latexmk/BibTeX/TeX Live 2026 complete successfully. Final TeX log has no undefined references/citations, overfull or underfull boxes, or LaTeX/package warnings. All pages inspected in rendered overview; main finite pages and dense appendix pages also inspected at reading resolution. Final local edits re-rendered and inspected. |
| Deliverable and preservation | Fresh PDF and manuscript sources remain in this workspace. No merge, push, benchmark campaign, or frozen proof-source modification. |

## Verification performed

`python -B paper/check_finite_integration.py` passes. It checks the five
table margins, plot coordinates, all three component bounds at the upper
rungs, the Q1-dominance claim, six timing entries, setup and memory figures,
and the selected-map transcription. It also enumerates all 524,288 inner
states, checks rank and CA=0, and reconstructs the kernel spectrum exactly.

`python -B workstreams/bch_rm2sub_bridge/verify_paper_milestone.py
--require-local-evidence` passes: 788 authenticated pins, zero missing local
dependencies, and exact retained union checks. Of those dependencies, 515
are not in Git. This check is not a new interval-arithmetic replay. The paper
accurately describes the retained 512/768-bit replay evidence and the
separate external artifact packaging obligation.

`git diff --check` passes. Changed tracked sources are confined to the paper
and the historical-plan pointer. A PDF text check found no unresolved `??`
references or TODO markers. No claim of new encoder benchmarking or complete
proof-assistant verification is made.

Temporary page renders and layout-check intermediates remain under
`tmp/pdfs/` because the cleanup command was rejected by the execution policy.
They are not manuscript deliverables and should not be committed.

## Recommended next review

Review the parameter slices on pages 21--22 first for engineering clarity,
then the BCH-256 theorem and Appendix B for claim scope and proof
exposition. A submission pass can decide how much artifact detail to retain
in print and prepare the separately distributed replay bundle.
