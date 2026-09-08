# Structured SPIN: paper revision plan

Updated: 2026-09-07. Status: manuscript revision implemented and PDF verified.
See `REVISION_REVIEW.md` for the requirement-by-requirement completion record.

This plan applies to `paper/main.tex` in the authoritative GitHub workspace.
The root-level `PAPER_RESTRUCTURE_PLAN.md` describes an older manuscript and
does not govern this revision. Numerical inputs and evidence are collected in
`../workstreams/bch_rm2sub_bridge/PAPER_HANDOFF.md`.

## Narrative

Present one structured SPIN architecture, then explain its asymptotic and
finite-length instantiations. Both use structured recursive inners from the
RM2Sub family. Introduce their common interface once; specify the selected
maps where each theorem needs them, without implying that identical dimensions
mean identical maps.

The asymptotic construction uses a growing Golay–BA-3 outer. Its role is to
establish the asymptotic distance guarantee through suitable outer spectral
control and the routing/inner analysis. Do not imply that outer minimum distance
or growth alone suffices. Nor should the finite BCH theorem be presented as an
asymptotic theorem for a fixed BCH constituent.

The finite construction uses a fixed algebraic BCH constituent and rigorous
spectrum inequalities. Exact spectra are particularly convenient, but are not
required: the selected BCH [256,128,d >= 38] constituent has an unknown exact
spectrum, and its certificate uses proved bounds instead.

Explain the choice of a fixed outer in the transition between these results:
reusing a sampled outer across rows requires control of its realized spectrum.
An ensemble mean alone need not control the higher powers that enter the
count. Concentration or selection can therefore consume finite failure budget.
This is a difficulty in obtaining sharp finite proofs, not evidence that
randomized outers are unsuitable.

Mention the existing finite randomized-outer construction briefly as evidence
that this route is viable. Keep this acknowledgment in the main narrative:
no separate subsection, formal claim, implementation comparison, or description
of its historical inner. Do not introduce a random-inner detour or discuss
the discarded Toeplitz design.

## Revision order and file map

1. **Preserve and clarify asymptotic structured SPIN.**
   In `structured_spin.tex`, retain the generic construction and existing
   asymptotic theorem. Correct any claim that transposition alone rules out
   late activation: state the support constraints actually established.
   Replace the obsolete finite-BA placeholder with the transition above.
   Retain the asymptotic proof in `structured_appendix.tex`.

2. **Give finite-length structured SPIN its own section.**
   Expand `finite_certificates.tex`, immediately after the asymptotic
   structured section, under the title “Finite-length Structured SPIN.”
   Motivate the fixed outer, define the selected BCH constituent and
   (t,s) = (128,19) inner, then state the finite theorem.
   Define the setup probability space and recurrence conventions once:
   independent row and region permutations, fresh nonzero field multipliers,
   zero initial state, output before update, state carried across regions,
   no flush, and no parity fanout. A single setup is shared by all messages.

   For each listed K, set N = 2K and H = floor(N/10). State that the probability
   of any nonzero message producing weight at most H is below 2^-40.
   These are five certified instances, not coverage of every intermediate K.

   | log2 K | Certified margin, bits | Measured online time, ms |
   |---:|---:|---:|
   | 16 | 53.9443672720 | 0.560 |
   | 18 | 52.3463883689 | 2.322 |
   | 20 | 50.4482033129 | 11.259 |
   | 22 | 48.4706838718 | Not measured |
   | 24 | 46.4762221823 | Not measured |

   Define margin as minus the base-two logarithm of the certified failure
   upper bound, not the true failure probability or total application security.
   Give a short proof overview: single active row, sparse occupancies, dense
   occupancies, explicit treatment of the all-one word, and the final union.

3. **Add a mathematical finite-proof appendix.**
   Specify the constituent and selected maps precisely, and connect the
   implemented constituent to the spectrum bounds used by the certificates.
   Explain the spectrum inequalities, transfer bounds, occupancy coverage,
   and certified summation. Cite reproducibility artifacts as evidence, not
   substitutes for the mathematical argument. Respect the distinct coverage
   methods used at smaller and larger rungs. Preserve frozen proof sources,
   ledgers, receipts, and benchmark kernels unchanged.

4. **Present the engineering curve with the finite result.**
   Plot the five certified margins against log2 K; connecting lines are visual
   guides. Explain the roughly one-bit loss per doubling over the upper rungs,
   the role of state size, and the remaining higher-occupancy constraints.
   This is a finite engineering observation, not an asymptotic extrapolation.
   Replace outdated principal guidance in `scaling_complexity.tex` rather
   than retaining conflicting parameter rules.

   Keep spectrum-based heuristic extrapolation brief and separate from the
   certified curve. Do not overlay the historical (64,20) heuristic as though
   it were a matched estimate for the selected (128,19) construction.

5. **Replace the implementation placeholder.**
   In `implementation.tex`, describe the measured no-fanout encoder and put
   workload, methodology, setup cost, and memory use together. The timings
   concern the transposed encoder: 2K input blocks to K output blocks, each
   128-bit block carrying 128 parallel binary instances. They are online
   encoder measurements, not full-application timings. The implementation
   currently supports exponents 16 through 20; do not imply measured or
   implemented support at 22 and 24. Use `PERFORMANCE.md` for provenance.

6. **Update the front matter last.**
   Revise `abstract.tex` and `introduction.tex` to describe the asymptotic
   theorem, finite BCH certificates, and measured implementation as connected
   contributions. Update organization and remove stale finite-result TODOs.
   Avoid repeatedly restating differences already defined in the construction.

## Acceptance checks

- Follow Controlled Writing for Cryptography: motivation, construction,
  property, and interpretation in that order; stable notation and probability
  scope; no unnecessary historical digressions or repeated footnotes.
- Check every reported margin, parameter, and timing against its authoritative
  artifact. Distinguish provenance checks from a fresh numerical proof replay.
- Check that no theorem assumes an exact BCH-256 spectrum, transfers a bound
  between different selected maps, or silently extrapolates beyond the five K.
- Preserve the asymptotic theorem's actual hypotheses and scope.
- Build a fresh `output/pdf/spin_codes_draft.pdf`; resolve references and
  citations, then visually inspect the revised pages and any layout warnings.
- Deliver revised sources, the fresh PDF, and a concise summary of changes.
  No new parameter search, benchmark campaign, merge, or push is needed.

## Next action

The approved parameter-plot follow-up is implemented: Section 7.1 now has
K/B, s/t, and K/s slices from the retained matched/nested BCH-64/128 study.
The BCH-256 certificate remains a separate fourth figure. The reproducible
generator, source precision, and check commands are described in README.md.

Author review of Section 7 and Appendix B: assess the main-text narrative and
the balance between printed derivations and artifact detail before submission.
External artifact publication remains a separate packaging step, not an
unreported prerequisite of this completed manuscript revision.
