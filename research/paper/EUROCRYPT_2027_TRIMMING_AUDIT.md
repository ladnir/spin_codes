# EUROCRYPT 2027: revised outline and trimming audit

For the current submission status, checked requirements and fresh-reader review,
see [the readiness audit](EUROCRYPT_2027_READINESS_AUDIT.md). Page counts below
are historical snapshots unless explicitly marked as the latest pass.

Date: 2026-09-17. Baseline: commit `4ec2bf06`, with the local
application-focused title and revised AI disclosure. This document records
the editorial plan and completed relocation passes; it is not a new
correctness audit of the mathematical arguments.

## Submission constraint and measured baseline

The [official instructions](https://eurocrypt.iacr.org/2027/papersubmission.php)
allow at most 27 main-paper pages excluding references, followed by clearly
marked supplementary material. Reviewers need not read the supplement;
the main paper must contain the core technical contribution. The proceedings
limit is 30 pages total, so a shorter body also leaves more room for references.

The anonymous audit baseline (`submission.tex`) has **70 pages**:

- Main text starts on page 1 and continues onto page 43.
- References start partway down page 43 and end on page 45.
- Supplementary material starts on page 46; the AI disclosure is on page 70.

Thus 16 numbered pages currently extend beyond the main-paper limit.
Removing the author block and draft contents pages has already happened in
this measurement. It is not an additional saving available to the cut.
The anonymous PDF is an audit baseline, **not a submission-ready version**.

The estimates below use section destinations in the compiled anonymous PDF,
including their vertical positions, and the reference heading as the endpoint.
They are rounded layout estimates, not independent page counts: floats and
page breaks will move when text is relocated.

| Current material | Starts on page | Approximate page equivalents |
|---|---:|---:|
| Title and abstract | 1 | 0.9 |
| 1. Introduction | 1 | 5.5 |
| 2. Preliminaries | 7 | 1.3 |
| 3. Common first-moment framework | 8 | 1.9 |
| 4. Accumulator SPIN | 10 | 3.9 |
| 5. Random SPIN | 14 | 2.4 |
| 6. Structured SPIN | 16 | 6.9 |
| 7. Finite-length Structured SPIN | 23 | 7.8 |
| 8. Finite lengths, asymptotics, and complexity | 31 | 2.7 |
| 9. Encoder implementation and performance | 34 | 4.3 |
| 10. Polynomial commitments | 38 | 3.0 |
| 11. Flock | 41 | 1.7 |

## Editorial decision

Center the main paper on the result promised by the new title:
fast codes for correlation generation and polynomial commitments, with an
explicit structured construction, a credible distance proof, finite
certificates, and measured application benefits.

The current text spends roughly six pages deriving the two supporting
construction variants before reaching the implemented structured family.
At the same time, the main structured asymptotic theorem points to Appendix A
for its substantive distance proof. Cutting only exposition would preserve
the wrong balance. **Move supporting derivations out and bring the central
proof mechanism in.**

Accumulator SPIN and Random SPIN serve primarily to build intuition, not as
coequal contributions. Accumulator SPIN is the minimal example and may also
be interesting independently for its simplicity. Random SPIN isolates the
effect of randomized recursive mixing. Retain this progression into structured
SPIN; move standalone theorem development and detailed proofs, not the
reasoning that motivates the main construction.

PCG, PCS, and Flock demonstrate the usefulness of SPIN through integration
into existing systems. Present them together as application evidence, not
as three new protocol contributions. Cite existing machinery and explain
only the required interface, the substitution, its parameter scope, and
the resulting performance and tradeoffs.

## Revised main-paper outline and 24.5-page budget

These allocations reflect the later decision to put the engineering study
and all but two tables in the supplement. The pass log records measured
page counts; these are approximate section budgets.

| Main-paper section | Budget |
|---|---:|
| 1. Introduction (including title, abstract, contributions, and related work) | 3.0 |
| 2. Preliminaries and the first-moment framework | 2.5 |
| 3. Building intuition: Accumulator SPIN and Random SPIN | 1.5 |
| 4. Structured SPIN: construction, routing, IMT, and transpose | 4.0 |
| 5. Asymptotic distance and linear-time encoding | 6.5 |
| 6. Selected finite instantiations and guarantees | 2.5 |
| 7. Encoder implementation and performance | 1.5 |
| 8. Applications: PCG, PCS, and Flock | 2.5 |
| **Target body** | **24.5** |

The remaining 2.5 pages are reserve, not space that must be filled. Use
them only where construction or proof clarity needs more room. Do not use
smaller fonts, narrower margins, or compressed spacing to meet the limit.
The larger central-proof allocation reflects pass 3. Savings now come from
consolidated applications and relocating the implementer-oriented study.

### Section-level contents

1. **Introduction.** State the code-design problem, main guarantees, and
   measured benefits. Present supporting variants as intuition and the
   applications as evaluations of the main construction.
2. **Preliminaries and framework.** Establish the encoder/transpose interface,
   setup randomness, and the first-moment calculation used throughout.
3. **Building intuition.** Start with the minimal accumulator construction,
   then explain the role of randomized recursive mixing. Keep defining
   recurrences and short calculations that teach the mechanism. End with the
   performance motivation for structure, not a second contribution summary.
4. **Structured SPIN.** Define the outer, randomized routing, IMT, and
   transposed encoder. Explain each design choice once.
5. **Asymptotic guarantee.** State the theorem and give the connected distance
   argument and linear-work justification in the main text.
6. **Finite instantiations.** Give the selected BCH guarantees, including
   quarter rate, and the finite proof overview. State the headline margins
   in prose. Put the complete margin tables in C and the whole parameter
   explanation, including all plots, in D; provide a main-text pointer.
7. **Encoder implementation and performance.** Keep the core transposed
   comparison as the principal table. Explain the implementation and timing
   policy concisely; summarize ordinary encoding and quarter rate in text.
   Full timing tables and methods go in E. Keep the PCG projection in Applications.
8. **Applications.** Use exactly three subsections, in this order:
   - **PCG:** required code interface, SPIN substitution, parameter scope, and
     projected correlation/OT throughput based on measured encoding costs.
   - **PCS:** required code interface, SPIN substitution, conditional security
     parameters, and measured commitment/opening results.
   - **Flock:** integration of the preceding PCS, measured end-to-end results,
     and proof-size/verification tradeoffs against the implemented baseline.
   Keep the standalone PCS comparison as the second main-body table.
   Report headline Flock results and tradeoffs in prose, with its complete
   table in F. Share methodology rather than repeating it; cite existing
   protocol descriptions instead of rederiving them.

## Appendix outline and content destinations

The letters below are planned destinations, not current appendix numbers.
Every appendix starts by identifying the main-text section it supports.
The main body retains the central proof argument; the supplement supplies
the complete supporting derivations, witnesses, and experimental detail.
No existing result is discarded merely because it is not central.

### A. Supporting constructions and their proofs

- **A.1. Uniform-interleaver analysis and random-block enumerators.** Keep
  the full uniform-interleaver specialization from `framework.tex`, the
  random-block lineage and enumerators from `accumulator_spin.tex`, and
  the injection law used by both supporting variants. This is the shared
  prerequisite for A.2 and A.3; retain all attribution.
- **A.2. Accumulator SPIN.** Full construction, exact accumulator enumerator,
  tail and sharp contraction bounds, asymptotic theorem, and proof from
  `accumulator_spin.tex`. Preserve its independent interest as a minimal
  construction, without giving it equal billing in the main body.
- **A.3. Random SPIN.** Full recurrence, exact transfer, tail bound, and
  asymptotic theorem from `random_spin.tex`, followed immediately by the
  exponent closure currently in `proof_appendix.tex`: limiting convolution
  growth, linear-weight closure, and uniform sparse-weight closure. Include
  the Silver/Expand-Convolute lineage and direct random-family cost analysis.
- **A.4. The BMS connection and termination.** Preserve the complete
  buffer/flush realization and qualified interpretation from
  `scaling_complexity.tex`. Keep the short significance statement in main
  Section 3; the termination qualification remains there too.

### B. Complete asymptotic Structured SPIN analysis

- **B.1. One sampled outer and spectrum control.** The Golay--BA-3 spectrum,
  good-outer selection event, and reuse/conditioning argument from
  `structured_appendix.tex`. The main proof retains the logic of selecting
  one outer; this subsection retains its complete estimates.
- **B.2. Structured routing.** Full route-domination argument and measure
  comparisons from `structured_appendix.tex`.
- **B.3. IMT maps and positive transfers.** Exact selected-map table,
  transfer definitions and domination proofs from `structured_imt_appendix.tex`.
  Preserve `tab:imt-maps` as the canonical definition used by finite proofs.
- **B.4. Occupancy regimes.** Complete positive-occupation, uniform sparse,
  and fixed-occupation arguments, including the continuum limit, one/two-row
  cases, weight-coupled bound for fixed occupancy at least three, constants,
  and certificate coverage. Keep the original dependency order unless a
  revised order is explicitly checked.
- **B.5. Completion, length schedule, and cost.** Preserve the full union
  argument and any length/complexity details moved from `structured_spin.tex`
  and `scaling_complexity.tex`. Main Section 5 still states the actual
  schedule, requested-length wrapper, and linear-work justification. Preserve
  the distinction between analytic arguments and checked arithmetic witnesses.

### C. Finite certificates at rates one half and one quarter

- **C.1. BCH-256 constituent and spectrum inequalities.** The fixed code,
  construction data, and deterministic inequalities from `finite_appendix.tex`.
  Do not replace them with an assumed or heuristic exact spectrum.
- **C.2. Selected IMT maps and transfers.** Finite transfer specialization,
  referencing the shared map definition in B.3 rather than duplicating it.
- **C.3. One-row and sparse bounds.** One active row, sparse occupancies,
  constant-row compositions, and all-one-word treatment.
- **C.4. Dense coverage and short-length refinements.** Change of input law,
  finite rectangles, composition factors, routing refinements, retained shell
  mass, fixed-tilt bounds, and label splits. Keep the complete refinements
  needed for the two shortest selected instances.
- **C.5. Quarter-rate certificate.** Exact BCH-128 subcode construction and
  spectrum attribution, separate feedback-map table, both distance cutoffs,
  sparse and dense type coverage, and the corresponding proof.
- **C.6. Exact union and completion.** Component-margin table, complete
  occupancy coverage, rational summation, injectivity, and proof completion.
  Keep local receipt references; put the consolidated replay index in G.

This appendix retains all mathematical content of `finite_appendix.tex`.
It also receives the complete finite-margin tables removed from
`finite_certificates.tex`. The selected theorem statements, numerical targets,
and headline margins remain in the main finite section.

### D. Extended parameter study and engineering curves

- **D.1. Length and outer-size slices.** Preserve the complete BCH-256,
  BCH-64, and BCH-128 fixed/adaptive-step figure and full-refresh references.
  The main finite section points here without reproducing a plot.
- **D.2. Short-length cancellation mechanism.** Full local derivation,
  adaptive-step selection rule, and the additional-mixing observation from
  `finite_certificates.tex`; retain the entire explanation here for implementers.
- **D.3. State and step-size slices.** Both existing state/step and
  state/message-length figures, map-chain definitions, and interpretation.
- **D.4. Certified engineering curve.** Complete full-versus-Q1 plot for the
  selected construction, with its matched-map and tested-length scope.
- **D.5. Parameter selection and costs.** Detailed selection criteria,
  admissibility constraints, and parameter tradeoffs from
  `scaling_complexity.tex`. Keep exploratory Q1 scores distinct from full
  certificates and measured performance throughout.

All four currently included engineering figures remain in the compiled
paper or supplement. A compact main-text extraction does not replace its
full source figure. Figures already inactive before this rewrite are not
silently restored or deleted.

### E. Encoder implementation and benchmark details

- **E.1. Kernels and layouts.** Full kernel choices, generated circuits,
  routing/index formats, tile choices, and memory footprints from
  `implementation.tex`; the core implementation explanation stays in Section 7.
- **E.2. Correctness and certificate binding.** Symbolic and dense-oracle
  checks, adjoint tests, exact map identity, and source/binary bindings.
- **E.3. Timing protocols.** Complete transposed and paired ordinary/transpose
  methodologies, setup exclusions, buffer policies, trial counts, and control
  measurements. Do not merge distinct measurement campaigns into one protocol.
- **E.4. External encoder baselines.** Dependency versions, native interfaces,
  actual lengths, permutation choices, distance provenance, original versus
  chosen-block BAA, the legacy EC profile, and RAA's separate reference length
  and setup-testing qualification. Keep qualifications required to read the
  main comparison beside that comparison too.

Only the cross-construction transposed comparison stays in the main body.
The selected-timing and ordinary/transposed tables are retained here with
full supporting details. Generators and receipt locations are indexed in G.

### F. Application details and additional comparisons

- **F.1. PCG projection.** Full accounting for the encoding/GGM estimate and
  contextual external timings from `implementation.tex`. Preserve which
  timings are measured and which throughput is projected; do not add an
  unperformed integrated protocol benchmark.
- **F.2. PCS interface and parameters.** Full row/column construction,
  binary packing, extension-field use, folded messages, parameter calculation,
  and conditional-security scope from `pcs.tex`. Keep the security scope in
  main Section 8.2 as well, not only here.
- **F.3. Standalone PCS methodology and breakdowns.** Input shapes, phase
  costs, allocation/serialization boundaries, Ligerito profiles, and trial
  protocols. The headline measured comparison stays in Section 8.2.
- **F.4. Flock adapter and evaluation details.** Weighted opening interface,
  two-opening parameters, checked relation and wrapper scope, shared prover
  optimizations, timing protocol, and phase breakdowns from `flock.tex`.
  Keep end-to-end results and proof-size/verification tradeoffs in Section 8.3.
- **F.5. Bolt comparisons.** Measured commitment comparison, opening-component
  measurements and calibrated projections, full opening table, and the
  optimistic Flock substitution. Preserve unmodeled costs and the distinction
  between measured results, projections, and runtime bounds.
- **F.6. Blaze compatibility.** Preserve the existing compatibility argument
  and its unimplemented/unparameterized scope; it is not a benchmark claim.

Keep the full Flock table here, including its measured/projected labels.
The main Flock subsection reports measured headline results and tradeoffs
in prose. No split or new generated table is needed.

### G. Artifact and reproduction guide

- **G.1. Claim-to-artifact index.** Consolidate the current introduction's
  artifact paragraph and source/receipt pointers scattered across the proofs,
  implementation, and applications. Map each theorem, certificate table,
  plot, and timing table to its producer, checker, and reproduction guide.
- **G.2. Proof replay.** Record pinned inputs, arithmetic semantics, coverage
  and exact-union checks, source authentication, and the distinction between
  higher-precision replay and an independent proof implementation.
- **G.3. Figure and timing regeneration.** Index the existing build scripts,
  campaigns, dependencies, and correctness checks. Retain the source-only
  checkout limitation: omitted generated inputs must be supplied or regenerated.
  Do not commit raw experiments merely to make this index complete.

This is a navigation/reproduction appendix, not a replacement for mathematical
definitions or proofs. B, C, E, and F retain local references to the specific
evidence they use. Author-identifying links follow the submission build's
anonymity policy; the nonanonymous artifact remains fully cross-referenced.

### H. Disclosure of generative AI use

Retain the concise `ai_disclosure.tex` unchanged by the trimming pass, along
with the prominent main-text reference to it.

Subsequent author-requested relocation: the disclosure now appears in the
main body immediately before the references. Its label is
`sec:ai-disclosure`, replacing the inventory's `app:ai-disclosure`; the
abstract reference follows the new location. The disclosure text is unchanged.

## Content-preservation ledger

The tables below specify the final destinations. The pass log records which
ones are already implemented. Paths are relative to `paper/`. A relocation
is complete only when its destination is included in the compiled manuscript.

| Current source/content | Main-body destination | Supplement destination |
|---|---|---|
| `abstract.tex`, `introduction.tex`: unique claims, lineage, related work | Abstract and Section 1; intuition in Section 3 | Full family comparison in A; detailed artifact navigation in G |
| `preliminaries.tex` | Section 2 | None: retain needed definitions in main |
| `framework.tex` | Section 2: common first moment and proof | A.1: expanded uniform-interleaver specialization |
| `accumulator_spin.tex` | Section 3: minimal example and intuition | A.1--A.2: all full enumerators, bounds, theorem, proof, attribution |
| `random_spin.tex` | Section 3: random mixing and motivation | A.3: full construction, transfer, theorem, proof, attribution |
| `proof_appendix.tex` (current Appendix C) | None separately | A.3: exponent closure; G.2: final artifact-interface paragraph |
| `structured_spin.tex` | Sections 4--5: construction, theorem, schedule | B.5: expanded length and complexity detail, if removed from main |
| `structured_appendix.tex` (current Appendix A) | Section 5: central outer and routing proof argument | B.1--B.2: complete derivations |
| `structured_imt_appendix.tex` (included within current A) | Section 5: transfer mechanism and occupancy union | B.3--B.5: maps, full regime bounds, constants, completion; G.2: replay index |
| `finite_certificates.tex`: selected construction, guarantees, tables | Main finite section: statements and headline margins | C: full margin tables; D: entire parameter narrative |
| `finite_appendix.tex` (current Appendix B) | Section 6: proof overview | C.1--C.6: complete existing proof; G.2: consolidated reproduction index |
| `scaling_complexity.tex`: native/requested lengths, asymptotics, cost | Main proof retains the native theorem, wrapper, and work argument | D: complete consolidated discussion |
| `scaling_complexity.tex`: BMS and random-family costs | Section 3: short qualified interpretation | A.3--A.4: costs and complete termination argument |
| `scaling_complexity.tex`: finite parameter selection | Section 6: explanation | D.5: full constraints and tradeoffs |
| `implementation.tex`: kernels, correctness, benchmark methodology | Section 7: core account and timings | E.1--E.4: full supporting details; G.3: reproduction index |
| `implementation.tex`: PCG implications | Section 8.1 | F.1: extended projection/context |
| `pcs.tex` | Section 8.2: substitution, scope, measured comparison | F.2--F.3: construction/parameters/methods; F.5--F.6: Bolt/Blaze |
| `flock.tex` | Section 8.3: integration, measured results, tradeoffs | F.4--F.5: adapter/methods/projection |
| `ai_disclosure.tex` (current Appendix D) | Existing prominent reference | H: full concise disclosure |

### Figures and tables: explicit destinations

| Existing figure/table label or source | Planned location |
|---|---|
| `fig:finite-k-b` / `figures/imt_parameter_k_b.tex` | Full three-outer plot in D; no main-text extraction |
| `fig:finite-s-t` / `figures/imt_parameter_s_t.tex` | D.3, unchanged data and scope |
| `fig:finite-k-s` / `figures/imt_parameter_k_s.tex` | D.3, unchanged data and scope |
| `fig:finite-bch-curve` / `figures/imt_certified_curve.tex` | D; main text states headline certified margins |
| `tab:spin-families` in `introduction.tex` | A, as a comparison including the structured endpoint; Section 1 keeps only the roadmap |
| `tab:imt-maps` | B.3, canonical shared definition |
| `tab:finite-bch-margins`, `tab:finite-quarter-margins` | C; theorem statements and thresholds stay in main |
| Unlabeled BCH generator-polynomial table near the start of `finite_appendix.tex` | C.1, preserving every construction constant |
| `tab:finite-quarter-feedback` | C.5 |
| Unlabeled component-margin table in `finite_appendix.tex` | C.6; add a stable label during relocation |
| `tab:finite-bch-performance` | E |
| `tab:ordinary-encoding` / `tables/ordinary_encoding.tex` | E; headline direction comparison in main prose |
| `tab:transposed-comparison` / `figures/transposed_comparison.tex` | Section 7, including distance/margin qualifications |
| `tab:spin-pcs-standalone` / `tables/pcs_standalone.tex` | Section 8.2 |
| `tab:pcs-opening` / `tables/pcs_opening.tex` | F.5, retaining SPIN measured rows as well as Bolt projections |
| `tab:spin-flock` / `tables/flock_comparison.tex` | Full table in F; measured results summarized in main prose |

Generated figure/table fragments must be changed through their producers,
not by hand. Preserve complete existing plots and data inputs when generating
the smaller main-text views.

### No-loss acceptance checks

1. Before moving text, save an identity-level inventory of all labels, theorem
   statements, figures, tables (including unlabeled ones), citations, and
   included source files. The current recursive `main.tex` include tree has
   **28 TeX files and 140 label occurrences**. These counts are a baseline
   sanity check, not a substitute for checking individual content.
2. Give every original substantive item a destination in the compiled main
   paper or supplement. A deleted duplicate must name the equivalent retained
   passage; never classify a unique assumption, proof step, result, limitation,
   or attribution as redundant merely to meet the page limit.
3. Preserve label identities where practical. Record explicit replacements
   for split figures/tables and moved section labels; check all references,
   duplicates, and theorem dependencies. Do not leave a proof depending on a
   definition that was removed during condensation.
4. Check row/point preservation for every split table or extracted plot.
   Keep all reported parameters, numerical values, evidence types, and
   certificate/timing bindings unchanged.
5. Verify that every appendix and moved float is reachable from the build
   and cited from the main text or its parent appendix. An unreferenced source
   file, old PDF, Git history, or repository-only proof is not preservation
   of the paper's content.
6. Compile both versions; inspect their contents, references, and moved
   floats. Check main-body length separately from supplement length. The
   relocation is complete only when this ledger is reconciled, not merely
   when the main body reaches the page target.

## Specific moves and protected content

### Introduction and framework

- In `introduction.tex`, combine the opening result summary and the five
  contributions: the same distances, lengths, and performance claims recur.
- Reduce "Three SPIN constructions" to a roadmap for the conceptual
  development below; move the full family table if it no longer earns its
  main-text space. Do not duplicate that development in the introduction.
- Merge "One distance calculation" into the framework rather than stating
  the same sum in two places. Shorten the organization/artifact paragraphs.
- Combine `preliminaries.tex` and the necessary parts of `framework.tex`.
  Keep the setup probability space, encoder/transpose convention, and the
  class-indexed first moment. Keep enough of the uniform-interleaver
  specialization to make the supporting examples concrete; move its extended
  derivation with their proofs.

### Intuition-building variants

- Retain approximately 1.5 pages of conceptual development after the common
  framework and before the structured construction. These are worked design
  steps, not three equally weighted standalone constructions. Their purpose
  here is understanding, not expanding the paper's list of contributions.
- For Accumulator SPIN, retain the outer--interleaver--accumulator interface,
  the accumulator recurrence, and the mechanism of its weight counting.
  Explain how this makes the common first-moment argument concrete. Its
  minimal structure may merit a short observation of independent interest,
  without restoring a full main-text treatment. Preserve the attribution to
  Block-Accumulate codes.
- For Random SPIN, retain the random-recursion definition and the key fact:
  fresh uniform feedback makes the next output bit uniform when the state is
  nonzero, whereas zero state behaves differently. Explain why tracking that
  state distinction supports a finite-state bound, and why coordinate-wise
  mixing under the analyzed schedule motivates a faster structured design.
- End with the transition to Structured SPIN: retain useful random mixing
  while exploiting constituent spectra, controlled routing, and batched inner
  operations. Explain what changes in the proof; do not imply the random-inner
  theorem automatically proves the IMT construction.
- Mention a supporting guarantee only where it advances the intuition, with
  its scope and a reference to the full theorem. Relocate standalone theorem blocks,
  enumerator derivations, full transfer matrices, and detailed proofs from
  `accumulator_spin.tex` and `random_spin.tex` to the supplement. A short
  calculation may remain when it teaches the mechanism used later.
- Preserve their theorem labels, attribution, and proof dependencies.
  In particular, the random proof uses the injection enumerator introduced
  with the accumulator construction; relocate them in dependency order.
- Retain the BMS connection, but keep its terminated-encoder qualification.
  Move the buffer-by-buffer streaming realization from `scaling_complexity.tex`
  beside the supporting proof rather than leaving a claim without its scope.

### Central construction and proof

- In `structured_spin.tex`, continue the conceptual development rather than
  restarting it. Compress repeated motivation and summaries of the dataflow,
  but preserve the reasons for the routing and inner choices. Define the route
  and recurrence once.
- Keep the exact setup law, reused-outer distinction, output-before-update
  convention, persistent state, and no-flush endpoint. These affect the proof.
- Preserve the exact nonzero-state transvection mixture and the transpose
  recurrence. Keep the finite/asymptotic outer distinction explicit.
- Promote a connected proof account from `structured_appendix.tex` and
  `structured_imt_appendix.tex`: select a good shared outer once; condition
  on its realized spectrum; derive the route/inner envelope; cover fixed,
  growing sparse, and positive occupations; conclude a summable union.
- Include representative transfer mathematics, not just a list of verifier
  files. In particular, explain why the low-dimensional state bounds dominate
  the true process and how the regimes fit together.
- Retain the theorem's actual length schedule and linear-work argument.
  Keep long matrix listings, exact witness tables, and detailed interval
  coverage in the supplement. A successful replay does not replace the
  analytic reductions.
- Absorb the nonduplicative parts of `scaling_complexity.tex` into this
  construction/proof account and the finite discussion. Do not retain it as
  another long section repeating all three families.

### Finite guarantees and engineering explanation

- Start with the selected BCH construction and certified guarantee, then
  explain the parameter choices. Currently the extended Q1 exploration
  precedes the selected finite theorem.
- Keep the five-length guarantee, unknown-spectrum method, all-one-word
  treatment, and complete occupancy coverage story. Move the detailed
  certificate table to C; state headline margins in main prose.
- Keep a compact quarter-rate statement and its two certified thresholds.
- Preserve the entire engineering explanation in D: short regions permit
  retained-state cancellations; shrinking `t` helps; increasing state alone
  does not remove that mechanism. Keep all four plots, the additional-mixing
  observation, and reproducibility sources. The main body needs a pointer,
  not a shortened engineering section.
- Keep Q1 distinct from full margins. The alternative small-step points
  must not inherit the selected `(128,19)` certificates or timings.

### Encoder evaluation and the compact Applications section

- Keep encoder design and standalone timings in Section 7, separate from
  the application evidence. Keep the same-host comparison and evidence boundaries: original
  versus chosen-block BAA, differing distance targets, and the RAA size caveat.
- Consolidate hardware/timing descriptions into a concise shared account;
  keep experiment-specific differences where they affect interpretation.
  Move exhaustive flags, tile choices, receipt names, and campaign details.
- Consolidate all application material into Section 8 with the three
  subsections PCG, PCS, and Flock. Move PCG implications out of
  `implementation.tex`; condense `pcs.tex` and `flock.tex` into the other two
  subsections. Do not repeat standalone encoder tables in Applications.
- Retain measured ordinary/transposed costs in the encoder evaluation.
  In the PCG subsection, distinguish projected OT throughput from an
  integrated PCG benchmark and state the additional costs in the projection.
- Keep the required PCS interface and conditional parameter scope. Cite the
  existing PCS construction rather than reintroducing its full machinery. Do not
  turn a fixed-matrix calculation into a composed security claim while cutting.
- Keep standalone PCS and Flock measured comparisons, including proof size,
  verification costs, and the shared baseline optimizations.
- Treat Flock as integration evidence for the preceding PCS, not another
  construction exposition. Move detailed application parameter calculations
  and integration choices to the supplement while preserving the assumptions
  needed to interpret the main-text results.
- Move the calibrated Bolt opening/Flock projections and the unimplemented
  Blaze discussion to the supplement. Keep necessary related-work credit.
  This focuses the main comparison on implemented, measured configurations.

## Execution order and acceptance checks

1. Finalize the structural plan and save the content inventory before moving
   scientific material. Use the appendix destinations and preservation ledger
   above to track every relocation.
2. Write the compact conceptual development, then relocate the supporting
   variants' full treatments and update cross-references. Compress duplicate
   motivation in the introduction/framework and structured opening. Build and
   record the actual page saving.
3. Rebalance the central proof, drawing from existing arguments rather than
   adding new claims. Integrate complexity and the finite guarantee.
4. Consolidate PCG, PCS, and Flock into the compact Applications section.
   Keep encoder design and standalone benchmarks separate. Condense the
   engineering plots and experimental narrative, preserving all details in
   clearly labeled supplementary sections.
5. Build the anonymous version and check: at most 27 pages of main text,
   no unresolved references, coherent first-use order, unchanged certificate
   and timing bindings, and no identifying metadata or artifact links.
   Reconcile every preservation-ledger entry against the compiled main paper
   and supplement; page count alone is not completion.

The initial anonymous build has empty author metadata and no GitHub URI
annotations. That is a limited check, not a complete anonymity audit.

## Pass log

### Pass 1: intuition and Appendix A (2026-09-17)

Implemented the first structural move, without changing any theorem,
certificate, numerical result, or timing.

| Ledger item | Completed destination / remaining work |
|---|---|
| Main intuition | New `spin_intuition.tex`: accumulator counting, state-dependent random mixing, minimal-construction significance, transition to structured encoding |
| Uniform-interleaver specialization | A.1, `supporting_framework.tex`; main framework keeps the common theorem and a short pointer |
| Shared outer lineage and enumerators | A.1, `supporting_outer.tex`, including the injection law before either supporting theorem |
| Accumulator construction and proof | A.2, `accumulator_spin.tex`; complete formal content retained |
| Random construction and proof | A.3, `random_spin.tex`, including direct encoding costs moved from `scaling_complexity.tex` |
| Random exponent closure | Included immediately within A.3 from `proof_appendix.tex`; its final artifact-interface paragraph remains there until G is assembled |
| BMS termination and interpretation | A.4, `supporting_bms.tex`; structured inner's linear-work explanation remains in the main complexity section |
| Family comparison table | Opening of A, same three rows and original label; introductory roadmap now points there |
| Existing structured/finite appendices | Still fully included, now B and C; their internal reorganization is pending |
| AI disclosure | Still fully included as current D; becomes H only when planned D--G exist |

The original label/source inventory is retained in
`TRIMMING_CONTENT_INVENTORY.md`. A recursive include-tree comparison found:

- all 140 original label identities present, with no duplicates (143 total
  after adding three navigation labels);
- all 101 original theorem, lemma, corollary, definition, proof, numbered
  equation, and align blocks preserved verbatim up to whitespace;
- all 120 original unnumbered display blocks and all 13 tabular blocks preserved;
- no missing citation keys; all four existing engineering figure sources
  remain included and unchanged.

The comparison is a content-preservation check, not a fresh proof review.
Only redundant navigation prose was replaced; no unique scientific claim
was discarded. The source inventory covers 28 original TeX inputs; the new
include tree contains 33.

Both PDF builds succeeded with no undefined references, duplicate labels,
or overfull boxes. The new intuition, appendix opening, moved proof pages,
and section transitions were visually checked. All 34 `test_imt*py` tests
passed, and `git diff --check` passed.

The anonymous build now has 71 total pages, with main text ending and
references beginning on page 36, versus page 43 in the baseline. Appendix A
starts on page 39. The new intuition occupies about 1.5 pages. Total length
has not been optimized: the full proofs are retained and the explanatory
section adds text. The main body is still above the submission limit.

Pending: condense duplicate introductory/construction motivation, promote
the central structured proof, consolidate applications, and implement
appendices D--G. Do not mark those ledger entries complete merely because
Appendix A is done.

### Pass 2: repeated motivation and introduction (2026-09-17)

Condensed `introduction.tex` and the opening of `structured_spin.tex`.
The introduction now groups the main results once, treats the two supporting
variants as intuition, and presents PCS/Flock as application evidence.
The structured opening explains the three design choices without repeating
the full random-construction story. No construction definition or formal
claim was changed.

| Condensed or relocated material | Retained content and destination |
|---|---|
| Introductory architecture/motivation repeated in the old roadmap | Introductory composition and interface; Section 4 intuition; structured design rationale |
| Accumulator and Random SPIN quantitative summaries | Intuition section and Appendix A.2--A.3 retain the rates, distances, schedules, costs, and proofs |
| Silver/Expand-Convolute history from the old roadmap | Introductory related work retains the refuted conjecture, coupling/concentration proof lineage, and proved-versus-empirical parameter distinction; A.3 retains the exact inner attribution |
| Repeated asymptotic and finite numerical claims | One introductory result paragraph; full statements remain in structured and finite sections |
| Quarter-rate distances/margins from the contribution list | Full selected thresholds remain in `finite_certificates.tex`; introduction points to that subsection |
| Introductory Bolt commitment speedup | The complete measured comparison and its same-volume/separate-campaign scope remain in `pcs.tex`, destined for F.5 |
| Original introductory first-moment formula and explanation | Moved intact to the start of A.1 in `supporting_framework.tex`, retaining `eq:intro-weight-first-moment` |
| BMS claim repeated in the contribution list | Short qualified statement in the introduction, intuition in Section 4, and full terminated-encoder argument in A.4 |
| Long organization/artifact account | One paragraph retains section navigation, source/guide/map links, and the distinction between value reproduction, proof replay, and new timing runs |
| Structured opening's repeated random-mixing motivation | Intuition section; short bridge at the start of the structured section |
| Structured design's nonduplicative rationale | Retained: low-weight multiplicities rather than minimum distance alone; fixed-message within-block randomization; late activation; distinct-region routing; batched IMT and its actual state law |
| Dataflow explanation repeated after the diagram | Short arrow-by-arrow interpretation remains; detailed setup, independence, maps, and constituent-sharing distinctions are unchanged |

The appendix destinations and original inventory remain in force. The
generic framework still permits other declared outer/interleaver/inner
combinations; condensation does not extend the guarantees to unproved ones.
No unique experimental result was dropped from the compiled manuscript.

Checks against the complete pass-1 include tree retain all 143 labels,
all 101 formal blocks, all 123 unnumbered displays, all 13 tabular blocks,
and all 26 citation keys. Formal blocks and displays were compared up to
whitespace. There are no duplicate labels. Numerical figure/table sources,
proof bodies, the abstract, and the AI disclosure were not edited.

Both PDF builds succeed without undefined references, duplicate labels,
or overfull boxes. The revised introduction, structured opening, and A.1
were visually checked. All 34 `test_imt*py` tests and `git diff --check` pass.
The anonymous main text ends on page 32 (previously 36; original baseline
43); references begin on page 32, and Appendix A starts on page 35.
The complete anonymous PDF has 67 pages; the author draft has 69.

Next: rebalance the central structured proof between the main body and
Appendix B. Compact applications and the extended engineering appendices
remain pending. The manuscript is still above the 27-page main-body limit.

### Pass 3: central Structured SPIN argument (2026-09-17)

Added structured_proof.tex as the main section "Asymptotic distance and
linear work," immediately after the construction. CWC guided the order:
select the shared outer, fix its realization, compare routing laws, bound
each occupation regime, then sum before applying Markov's inequality.
This pass exposes the proof dependencies; it does not change the theorem
or strengthen any certificate.

| Material | Main body / supplementary destination |
|---|---|
| Scalable BA-3 definition, native schedule, theorem, and requested-length wrapper | Moved intact from structured_spin.tex into Section 6 |
| Shared-outer tail event, selection lemma, and proof | Moved from Appendix B into Section 6.2; the realized spectrum is controlled before taking powers |
| Route-domination lemma, full proof, and concavity inequality | Moved from Appendix B into Section 6.3 |
| Positive-occupation argument | Section 6.4 explains the weighted-state transfer, positive-vector bound, reference-density change, exponent, and uniform negative bound |
| Sparse occupation | Section 6.5 explains why the dense conditioning loss fails, the fair-row comparison, tilted marks, polynomial witness, and summable bound above the explicit cutoff |
| Fixed occupation | Section 6.6 explains isolated impulses, empty-interval mixing, the zero/live continuum kernel, coefficient extraction, and the finite set of remaining bounds |
| Final union and linear work | Original theorem proof moved intact into Section 6.7 |
| Exact outer spectrum and majorant calculation | Retained in Appendix B |
| Map table, transfer entries, syndrome-fiber bounds, and numerical witnesses | Retained in Appendix B |
| Full continuum-limit comparison and product-norm derivation | Retained in Appendix B; the main section states their roles and resulting inequalities |
| Transition from growing to fixed outers | Moved to the opening of the finite-length section |

All 143 previous label identities remain, with no duplicates (146 total
after adding navigation labels). All 101 previous formal blocks, 123
unnumbered displays, 13 tabular blocks, and 26 citation keys are retained.
The comparison is up to whitespace for mathematical blocks, not a new
independent verification of the underlying certificates. The include tree
now contains 34 TeX sources. The original inventory remains applicable.

Both PDF builds succeed. The central proof, finite-section transition,
and revised Appendix B opening were visually checked in the anonymous and
author versions. All 34 IMT tests and the Git whitespace check pass.
No certificate, measured timing, figure data, or AI disclosure was changed.

This is a structural pass, not a net page reduction. The anonymous main
text now ends on page 37, versus 32 after pass 2; the complete anonymous
PDF has 70 pages and the author draft has 72. The new central section
occupies roughly 6.5 pages, above its tentative five-page allocation.
Detailed supporting arguments remain compiled in the supplement.

Next: consolidate PCG, PCS, and Flock into the planned compact Applications
section, preserving protocol details and full experiments in Appendix F.
Then rebalance the finite engineering discussion and remove the separate
complexity section's duplication. A final compression pass on the central
proof may still be needed. Appendices D--G and the 27-page main-body target
remain unfinished; neither proof promotion nor successful compilation
establishes submission readiness.

### Pass 4: two main tables and an implementer-oriented supplement (2026-09-17)

Following the revised editorial direction, the main body retains exactly
two tables: the cross-construction transposed-encoder comparison and the
standalone PCS comparison. The latter shows both the measured prover
benefit and the proof-size/verification tradeoff. Flock's measured results
are stated in prose. No engineering plot remains in the main body.

| Material | Completed destination |
|---|---|
| Half- and quarter-rate theorems, selected constructions, proof overview, headline margins | Main Section 7 |
| Complete half- and quarter-rate margin tables | Appendix C.7; original tables, labels, and values retained |
| Entire parameter explainer and all four engineering plots | Appendix D, engineering_appendix.tex |
| Remaining length, state, admissibility, and cost discussion | Appendix D.3 via scaling_complexity.tex |
| Core transposed comparison and concise implementation/method account | Main Section 8 |
| Selected timing ladder, ordinary/transposed table, kernels, memory, controls, full methods, baseline details | Appendix E, implementation_appendix.tex |
| PCG projection, PCS interface and conditional parameter scope, Flock headline results | Main Section 9, applications.tex; exactly three application subsections |
| Standalone PCS comparison | Main Section 9.2; unchanged generated table |
| Full PCG accounting, PCS details, opening table, complete Flock table, Bolt projections, Blaze discussion | Appendix F; applications_appendix.tex includes pcs.tex and flock.tex |
| AI disclosure | Unchanged, now Appendix G until a separate reproduction appendix is added |

The quarter-rate theorem now names the two parameter pairs introduced
immediately before it, instead of requiring the reader to look up a table
row. Its inequality, quantifiers, construction, and both parameter pairs
are unchanged. This is the only previous formal block with a wording
change in this pass; the other 100 formal blocks remain verbatim up to
whitespace. All 135 previous unnumbered displays, 13 tabular blocks, and
26 citation keys remain. All 146 previous label identities are retained,
with no duplicates; there are 156 labels and 38 included TeX sources.
All original figure/table fragments and data inputs are unchanged.

The main narrative follows CWC's reader-state and scope rules: numerical
targets remain beside their theorems, the PCS conditional-security scope
remains beside its result, and measured versus projected costs stay explicit.
Nothing relies on Git history as its only preservation destination.
The artifact's PAPER_MAP.md now explains the source layout.

Updated the two plot-location regression tests and the finite manuscript
checker to follow the relocated sources. The checks still authenticate
the same evidence; they also require the appendix inputs to be included
after the main paper. All 34 IMT tests pass. The finite integration check
passes for seven certificates, four timing cells, 57 map words, 746
authenticated files, 130 parameter geometries, and 110 adaptive-length
cells. The asymptotic IMT integration check and application-table/exact
conditional-parameter check also pass. These are authentication and
transcription checks, not a new interval replay.

Both builds succeed without overfull boxes, undefined references, or
duplicate labels. Changed main sections, moved appendix tables, engineering
plots, and appendix transitions were visually reviewed. The anonymous main
text ends on page 24, down from 37; references begin on page 24.
The complete anonymous PDF has 74 pages and the author draft has 76.
The larger supplement preserves the full experimental account while the
main paper is now within the 27-page limit, without typography changes.

Next: read the main body end to end for flow, notation, and claim scope,
then perform the anonymity/submission audit. A consolidated reproduction
appendix is still pending, but the local reproduction pointers remain in
the compiled supplement and the artifact guide. Meeting the page limit
does not by itself establish submission readiness.

### Pass 5: display-math economy before further content cuts (2026-09-17)

Reviewed short displays in the main body, following the preference to
recover space from presentation before cutting scientific content.
Converted 17 unnumbered display blocks to inline math:

- Preliminaries: binary entropy, the encoder's type, and its image code.
- Common framework: the composition defining the realized encoder.
- Structured construction: fixed map signatures, constituent-map signature,
  three routing index rules, inner input grouping, repeated inner-map
  signatures, and the one-step recovery identity in the bijectivity proof.
- Asymptotic section: base-encoder and Golay signatures, shared-outer
  assignment, the relative-weight window, and the row-profile definition.

The formulas, domains, quantifiers, and construction choices are unchanged.
The numerical bounds, central recurrences, transvection law, adjoint,
theorem statements, and substantial transfer calculations remain displayed.
Every numbered equation and all 156 labels are retained. The only proof
block with a formatting change is the bijectivity proof's inline recovery
identity. No appendix material, plots, tables, or scientific content was cut.
No font sizes, margins, display-spacing lengths, or document-class settings
were changed.

Both builds succeed without overfull boxes, undefined references, or
duplicate labels. The changed definitions, construction, proof transitions,
and final main pages were visually checked. All 34 IMT tests and both
asymptotic and finite manuscript integration checks pass.

The anonymous main body now ends on page 23; references start on page 24.
The complete anonymous PDF has 73 pages and the author draft has 75.
The preceding table/engineering relocation is preserved. Further trimming
is not needed merely to meet the 27-page limit: prioritize an end-to-end
readability and claim-scope review, then the submission/anonymity audit.

### Pass 6: restore Bolt as an explicitly estimated PCS comparator (2026-09-17)

The headline PCS table now includes Bolt-max in the amortized limit:
2,993 ms projected commitment plus opening, 11.1 MiB modeled communication,
and no verification-time claim. Asterisks identify both estimated cells.
The note states the GF(2^32)-versus-Boolean input distinction and the
negligible per-opening matrix-proof assumption. The comparison matches input
bytes, not polynomial dimension or a common composed-security theorem.

Appendix F.2 and `artifact/BOLT_PCS_ESTIMATE.md` give the components and
limitations. The estimate includes the calculated outer authentication,
sizes from the same verified inner-proof proxies used for timing, and explicit
messages. It is not a serialized complete Bolt proof or a communication bound.
The earlier paper-anchored size scenarios serve only as a secondary check.
No new benchmark or changes to the pinned measurement data were needed.

Six new arithmetic/table tests pass, including exhaustive small-tree checks.
The 34 existing IMT tests and both manuscript integration checks still pass.
Both PDFs build cleanly, and the updated table and appendix were rendered and
visually checked. The anonymous main body ends on page 24 (references start
on the same page); total lengths are 75 anonymous and 77 author-draft pages.
The main body remains below the 27-page limit. The older Brakedown-encoder
comparison is deferred for discussion, as requested.
