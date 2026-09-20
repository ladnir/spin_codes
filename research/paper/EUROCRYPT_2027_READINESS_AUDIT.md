# Eurocrypt 2027 submission-readiness audit

Audit date: 2026-09-17. Base revision: `4ec2bf06`, plus the current uncommitted
paper work. This is the current decision record; the trimming audit preserves
historical page counts and earlier plans.

Follow-up, 2026-09-17: the author pointed out that the existing libOTe benchmark
already supplied the OT measurements. The paper now imports those results,
not the old projection. The initial audit identified an omission from the
paper, not missing performance evidence. Source and log pins are recorded in
`artifact/OT_PERFORMANCE.md`; the four medians were independently recomputed
from the retained logs. No benchmark was rerun. Build counts below describe
the original audit snapshot until the follow-up build record is appended.

## Outcome and scope

Artifact scope decision, 2026-09-17: the author requests only the core SPIN
implementation, build instructions, and correctness tests. Numerical proof
archives, parameter sweeps, and application/comparison experiments are outside
the artifact. Earlier recommendations to package all evidence are superseded.
The paper's mathematical claims and their supporting arguments are unchanged.

The central code contribution is identifiable and its asymptotic/finite
separation is clear. The highest-impact issues found in this pass concerned
application claim scope and presentation, not a demonstrated failure of the
distance theorems. The corresponding prose corrections are implemented.
This is a bounded review, not an independent reconstruction of every proof or
a declaration that all submission obligations have been met.

A fresh-context reviewer read the main body before the supporting evidence,
using Controlled Writing for Cryptography. It received no conversation history,
made no edits, and ran no benchmarks. The main agent independently checked
official submission rules, local evidence bindings, typography and anonymity.
The reviewer then checked the revised passages for overcorrection.

## Eurocrypt requirements

Sources checked live:
[submission instructions](https://eurocrypt.iacr.org/2027/papersubmission.php),
[call for papers](https://eurocrypt.iacr.org/2027/callforpapers.php), and
[AI policy](https://eurocrypt.iacr.org/2027/aipolicy.php).

| Requirement | Result / action |
| --- | --- |
| Deadline | September 17, 2026, 23:59 AoE; September 18, 04:59 PDT / 11:59 UTC. No submission was made by this audit. |
| At most 27 pages excluding references | Main text occupies pages 1--24. Pass. |
| Clearly separated supplementary material; core results intelligible without it | Supplement follows references with an explicit heading. The main text includes the construction, asymptotic proof and finite certificate framework; full transfers and witnesses remain supplementary. The cold review found the selected maps insufficiently introduced; corrected. |
| Current LNCS class, ordinary font/margins/spacing, visible page numbers | Updated the pinned class to the current publisher download; see `LNCS_TEMPLATE.md`. No margin or font compression. Both builds display page numbers. |
| Anonymous PDF | Submission has anonymous title block and empty author metadata. Full-PDF text/URI screening found no author artifact URL, local user path, or affiliation disclosure. References name prior authors normally. This is not an anonymity audit of a future source archive. |
| DOI links where possible; LNCS bibliography style | Uses `splncs04`. Added eleven verified DOI identifiers. The Allerton 1998 item and preprints retain their existing bibliographic form; no DOI was invented. |
| Prominent and proportionate AI disclosure | First-page footnote plus Appendix G cover all writing, proof development, implementations, iterative supervision, GPT-5.6/GPT-6.0, and author responsibility. Prompt-text expectation remains unresolved; see below. |
| Author list, conflicts, contribution and exclusivity declarations | Author-side HotCRP checks remain. Confirm automatic and additional conflicts, no prohibited parallel/duplicate submission, and the six-submission cap per author. |
| Accepted version and presentation | Proceedings limit is 30 pages total; prepare a separate short proceedings version if accepted. At least one author must present in person and consent to recording. |
| Artifacts | The CFP announces a later artifact-review call. The author selects a core-code artifact only. Any review-time code archive must be anonymous; full experiment/evidence packaging is outside the chosen scope. |

The AI policy asks for input-prompt text when substantial sections are generated.
The current appendix explicitly says the complete transcript is omitted because
of its volume. No exemption for volume was found on the official page. The
disclosure is prominent and candid, but this audit cannot certify full policy
compliance or predict enforcement. The policy also asks that the work accurately
represent the authors' underlying contributions rather than primarily the tool's
generative capabilities. Author-originated ideas and iterative supervision are
relevant facts, not a guarantee of how chairs will interpret that condition.
Do not fabricate prompts or silently weaken the disclosure. Any clarification
with chairs, or decision about available prompt material, belongs to the authors.

## Fresh-reader findings and disposition

| Priority | Finding | Implemented response / residual |
| --- | --- | --- |
| P1 | Headline 50M OT/s lacked reproducible accounting in the paper. | Resolved by importing the existing libOTe measurements: 19.531/20.625 ms without hashing, 22.099/21.841 ms with hashing. The main text and appendix state per-party timing boundaries and noise configuration. This does not claim a new structured-LPN security theorem. |
| P1 | "Conditional 100-bit PCS" could be read as a full security result conditioned only on a good code. | Abstract and application text now say interactive fixed-matrix testing. Setup failure, extraction and Fiat--Shamir reductions are excluded. The total composed bound remains unclaimed. Flock's adaptive-composition qualification remains. |
| P2 | The five image-weight shells appeared only when the proof needed them. | Introduced weights 48,56,64,72,80 and C's distinct nonzero weight-five spanning columns before the scalable theorem; explained their proof/computation roles. |
| P2 | New ingredients were not sharply separated from serial-concatenation background. | Added a concise attribution and the direct-work contrast with the internal Random SPIN example. A quantitative theorem-level comparison to the closest external construction remains a useful optional improvement. |
| P2 | Standalone PCS allocation boundaries differ. | Put the difference immediately after Table 2: reusable SPIN workspace versus timed Ligerito internal allocations. The comparison is of implementations, not an isolated encoder replacement. No allocation-cost ablation is claimed. |

The reviewer accepted the revised scope and found no new inconsistency. Its two
follow-up suggestions were applied: remove revision-history discussion of the
unsupported OT estimate, and say "direct encoding work" rather than imply an
implementation-independent lower bound for Random SPIN.

## Technical and measurement checks

No theorem thresholds, selected maps, numerical certificate data, timing inputs,
or generated tables were changed. The new map paragraph restates exact properties
already used by the appendix and checkers. A setup-failure upper bound is not an
attack probability estimate or an upper bound on actual security.

Commands run from the repository root:

```text
python -B -m unittest discover -s paper -p "test_*.py"
python -B paper/check_finite_integration.py
python -B paper/check_imt_integration.py
python -B paper/build_application_tables.py --check
python -B paper/build_imt_comparison.py --check
```

Results: all 52 tests passed; seven finite targets, four timing cells, 57 map
words and 746 evidence files authenticated; 130 Q1 geometries, 110 adaptive-length
cells and five matched full anchors checked; asymptotic 11% bindings and 38 map
words checked; application tables and exact conditional parameter calculations
passed. These are binding, transcription, regression and exact small-instance
checks, not a fresh full interval replay or a new full-grid certificate.

The earlier sequential proof audits added checks for the IMT adjoint, BA outer
tails and finite probability reductions. This pass ran those regressions too;
it did not reinterpret their small-instance tests as a full proof verifier.

The performance audit checked timing boundaries and evidence labels. No new
benchmark was run. Original BAA, chosen-block BAA, EC and RAA remain distinct;
the RAA reference margin is not transferred to unproved smaller sizes. Bolt
opening costs remain projections, not measurements. PCS byte-volume and
field/dimension differences remain explicit.

## Build and visual verification

Both `latexmk` builds completed with the pinned publisher class. The anonymous
PDF has 77 pages: main text 1--24, references 25--26, supplementary material
27--77. The author draft has 79 pages, including its two-page contents list.
There are no overfull-box, undefined-reference or multiply-defined-label
warnings. The existing amsmath `vec` warning and host Perl locale notice remain;
neither stopped a build. `git diff --check` passed on the tracked files changed
by this pass.

Rendered and visually inspected anonymous pages 1, 2, 12, 23--27, 70, 71 and 77,
and draft pages 1, 4, 5, 14, 15, 25, 27, 28 and 72. The changed prose, map
introduction, application table, bibliography, supplement boundary and AI
disclosure have no observed clipping or overlap. This is targeted visual QA,
not a claim that every unchanged appendix page was visually re-reviewed.
The full-PDF text and link screening found no TODO/FIXME tokens or identifying
artifact links. Ordinary bibliography citations to the authors' prior work remain.

## Remaining work, in priority order

1. Resolve author-side submission checks, including AI prompt-material policy.
2. If preparing the core-code artifact, check its build and usage instructions
   from a clean checkout and prepare anonymous access as needed. Do not expand
   this into reproducing all proof searches or application experiments.
3. Decide whether the application section should remain performance evidence
   or support a stronger security claim. The latter needs commitment extraction,
   Fiat--Shamir and adaptive Flock composition, not just more distance margin.
4. Optional: quantify standalone allocation overhead and strengthen the short
   comparison against prior theorem-level guarantees. Neither warrants delaying
   a correctly scoped codes paper for an open-ended implementation campaign.

The next recommended paper work is the final application-claim and submission
check. A clean core-code build is the bounded artifact task. No further OT
benchmark campaign, full evidence packaging, or inner optimization is needed
for this integration.

## OT follow-up verification

The measured OT results are now in the abstract, introduction, Section 9.1
and Appendix F.1 (Table 10). The source explicitly selects semi-honest mode;
the paper and artifact record that mode. All four medians were recomputed
from nine post-warmup trials per mode. All 52 paper tests and the application
table checks pass. Both PDFs rebuild without overfull boxes or unresolved
references. The anonymous PDF has 78 pages, still with 24 main-text pages;
the author draft has 81 pages. The extra OT methodology is supplementary.
Changed abstract, introduction, application and appendix pages were rendered
and visually checked in both builds, including the final mode-label correction.
No raw benchmark data or binaries were added to the paper tree.
