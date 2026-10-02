# ePrint readiness review — 2026-09-22

## Final author-version pass — 2026-09-23

Follow-up at Peter's request: manuscript filenames and document-level repository
links have been replaced by the repository citation. The only GitHub URL now
embedded in the author PDF is `https://github.com/ladnir/spin_codes`.
This resolves the broken-link issue below without publishing private measurement
archives. The author PDF was rebuilt (84 pages, 882,296 bytes); affected pages
were visually checked, and the IMT, finite-integration, and application-table
checks pass. Its current SHA-256 supersedes the earlier checksum below:

```
085d8a974f8b7d4155a92d7b35baca50b3105c769f856e660532c1887117eb8c
```

No raw manuscript filenames, unresolved citations, or overfull boxes remain.
Numerical results and benchmark timings are unchanged. No commit or upload has
been made; the next step is to freeze the reviewed source snapshot and upload
this author PDF.

Target confirmed by Peter: the author/ePrint PDF, not the anonymous conference
version. The final build is `output/pdf/spin_codes_draft.pdf`, 84 pages and
888,113 bytes, with SHA-256:

```
0dce7c90aae27f2acf3ceb392e68e2132fc97e9600b70a4688a6d378aea7d8dc
```

**One release issue remains:** the embedded precomputed-performance link in
Appendix E (page 73) returns HTTP 404:
`https://github.com/ladnir/spin_codes/blob/main/research/paper/PRECOMPUTED_PERFORMANCE.md`.
The target exists locally but is untracked. Publish that compact ledger with
the reviewed manuscript source before uploading, or remove the unavailable
hyperlink. No numerical archives or raw measurements need to be published.
The other eleven distinct repository links embedded in the PDF returned 200.
The reviewed manuscript changes are still uncommitted on top of `f2010e03`;
that commit alone does not identify this paper version.

With Peter's approval, Section 10 now says "substantial author-directed
revision" instead of claiming direct editing by the authors. This was the only
manuscript edit in the final pass. PCS and Flock results were left unchanged.

Checks performed:

- All five current timing/table/integration checks passed:
  `precomputed_results.py`, `build_imt_comparison.py --check`,
  `build_application_tables.py --check`, `check_imt_integration.py`, and
  `check_finite_integration.py`. The finite check covers eight selected
  certificates, seven timing cells, 57 map words, and the 746-file original
  evidence set, plus the additional small-length evidence.
- The abstract, application text, and tables agree on the 9.289 ms
  precomputed encoder and 89.0 million regular-noise hashed OTs/s. The
  encoding-only 185 million blocks/s and heuristic stationary result are
  separately scoped. The additional small-length configuration remains
  separate from the original finite theorem.
- All 27 bibliography entries were identified: 20 DOI records resolved
  through Crossref/publisher metadata; four ePrint-only entries matched
  official ePrint records; the remaining entries were verified against
  arXiv, the JST catalogue, and the public SPIN repository. This checks
  bibliographic identity, not every claim supported by each citation.
- All 84 pages were rendered and inspected in overview. Front matter,
  main performance tables, exact-map and margin tables, the main engineering
  plot, and the current OT table were additionally inspected at readable
  resolution. After rebuilding, the four changed pages (28–31) were
  reinspected. No clipping or overlapping content was observed.
- The build succeeds with no undefined references/citations or overfull
  boxes. The existing amsmath accent warning and underfull boxes remain;
  this is not a warning-free build. All PDF fonts are embedded. Metadata
  lists the four authors and the correct title. No review annotations,
  local filesystem paths, TODO markers, or unresolved `??` were found.
- Before the disclosure edit, all 41 local source inputs matched the
  latexmk build hashes after Windows text-mode newline normalization.
  The final successful latexmk build incorporates the approved edit.

Scope: final presentation, claim-consistency, bibliography, and evidence-binding
review, not an independent reproof or fresh interval-arithmetic replay. No
benchmarks, commits, pushes, or uploads were performed. The anonymous PDF was
not rebuilt. Next: publish the scoped paper/ledger snapshot, check the repaired
public link, and upload the author PDF identified above.

## Earlier review and follow-ups

Performance follow-up: the optimized precomputed campaign supersedes the
standalone encoder timings in the original review. At K=2^20, ordinary and
transposed encoding now take 9.031 ms and 9.289 ms, and the original-BAA
comparison is approximately 3.4x. OT, PCS, and Flock measurements are unchanged.
The finite integration checker has been updated to authenticate the new timing
campaign and the separate K=2^16 two-round certificate. All five current
performance/table/integration checks pass; no benchmark was run.

Follow-up: finding 1 is resolved. All ten manuscript deep links now use the
canonical `spin_codes` repository and `research/` paths; every distinct
destination returned HTTP 200. The introduction cites the repository through
the new `spinCodeRepository` bibliography entry (reference 25 in the rebuilt
author PDF). The 81-page PDF was rebuilt, its embedded GitHub links checked,
and the introduction and bibliography pages visually inspected. No undefined
references or citations remain. Finding 2 is now resolved by updating the build
and reproduction guides for the relocated tree. Finding 3 is withdrawn: Peter
confirmed that original BAA was measured on the server, and the earlier concern
relied on out-of-date notes. Finding 4 is also resolved: the undefined historical
control comparison was removed, retaining the ordinary and transposed timings
of 10.797 ms and 10.674 ms. The author PDF was rebuilt; pages 72–73 were visually
checked and the application-table check passed. All four findings are now
resolved or withdrawn. The original review and its verification scope are
retained below except for the corrected BAA finding.

Reviewed the working tree at `C:/Users/peter/repo/permute_conv-github-bch`.
The active manuscript is `research/paper/main.tex`. The existing author PDF
is `output/pdf/spin_codes_draft.pdf`, dated September 17, with 81 pages.
This review does not establish which PDF was actually submitted to Eurocrypt.

The manuscript is substantially prepared for an author release, but the public
navigation and build instructions need a migration pass. The technical checks
below pass. This was a readiness and claim-scope review, not an independent
verification of every proof or a fresh interval-arithmetic replay.

## Findings

### 1. Repair ten manuscript links after the repository move

Locations: `introduction.tex:152`; `engineering_appendix.tex:85,111,125`;
`finite_appendix.tex:475,620`; `implementation_appendix.tex:87,214`;
`structured_imt_appendix.tex:472,475`.

All ten deep links still use `ladnir/permute_conv/blob/main/` followed by
`artifact/` or `workstreams/`. The repository root redirects to
`ladnir/spin_codes`, but the linked paths now require `research/`.
The old artifact README and finite-results links return HTTP 404. Every
corresponding local target exists under `research/`.

Correction: use `https://github.com/ladnir/spin_codes/blob/main/research/...`
for the moved files, and update the repository entry link to the canonical
name. For a reproducible release, pin the final reviewed revision where
appropriate. Verify the public targets and the links embedded in the rebuilt
author PDF. The corrected public artifact guide is:
<https://github.com/ladnir/spin_codes/blob/main/research/artifact/README.md>.

This is a release defect, not a mathematical issue. Anonymous builds suppress
these hyperlinks, so it particularly affects the author/ePrint build.

### 2. Update the documented working directories and output paths

Locations: `research/paper/README.md:36,66,72,92,112,119`;
`research/artifact/README.md`, “Build and test the half-rate encoder”;
`research/artifact/REPRODUCING.md`, opening working-directory convention.

The paper guide still uses another author's absolute checkout path. Its
commands advertised as running from the repository root instead assume
`research/`. Its `../output/pdf` destination, when used from `research/paper`,
also differs from the current root README's `../../output/pdf` destination.
The core-artifact CMake command points to a nonexistent root `workstreams/`.

Correction: specify working directories consistently. From the current root,
the manuscript checks are:

```powershell
python -B research/paper/check_imt_integration.py
python -B research/paper/check_finite_integration.py
python -B research/paper/build_application_tables.py --check
```

The root README already gives the appropriate paper build layout:

```powershell
Set-Location research/paper
latexmk -pdf -interaction=nonstopmode -halt-on-error '-outdir=../../output/pdf' '-jobname=spin_codes_draft' main.tex
```

For historical research commands, an explicit `Set-Location research` or
`cd research` is preferable to blindly rewriting paths inside frozen producers.
Distinguish the published `spin/` library from the implementation used for
the submitted timing campaign.

### 3. Original BAA measurement — concern withdrawn

Peter confirmed on 2026-09-22 that BAA was run on the server. Earlier notes
suggesting that the comparison needs a new measurement or remains unresolved
are out of date. The abstract's approximately 3× comparison is retained.
`research/artifact/APPLICATION_RESULTS.md` now records this clarification.

### 4. Remove or identify the undefined “preceding inner” comparison

Location: `implementation_appendix.tex:133` and the following paragraph.

The ordinary-encoding discussion compares against “the preceding inner and
implementation” and reports improvements of about 5% and 8%. The manuscript
does not identify that historical implementation sufficiently for a reader to
know what changed. The paper guide explicitly says the current narrative
presents IMT without a preceding-inner comparison.

Correction: remove the historical control comparison and retain the directly
relevant ordinary/transposed costs, or identify the baseline configuration and
source explicitly. CWC-02, CWC-03, and CWC-05 apply: development history is not
reader context, and the comparison needs a stable referent.

## Evidence and scope worth preserving

- `main.tex` defaults to the author version. The inspected PDF has all four
  authors, affiliations, matching PDF author metadata, visible page numbers,
  and the contents list. `submission.tex` remains the anonymous entry point.
- The finite theorem is explicitly restricted to the selected five half-rate
  lengths. The quarter-rate theorem separates its two distance thresholds.
- The PCS text distinguishes conditional fixed-matrix testing error from setup
  failure, extraction, and Fiat–Shamir security. Preserve that qualification in
  both the abstract and the application tables.
- Bolt projections are visibly marked as estimates, with the comparison's
  input-field and allocation differences disclosed.
- The core artifact deliberately excludes numerical receipts and raw
  measurements. The manuscript discloses that boundary. Publishing a compact,
  separately versioned proof-evidence package would improve independent
  verification, but is an additional release choice, not an existing artifact
  promise or a demonstrated error in the proofs.
- The current standalone library includes post-submission optimizations.
  Its current timings must not silently replace the retained paper campaign.

## Verification performed

From `research/`, the following checks passed:

1. `python -B paper/check_imt_integration.py`: 38 map words and the retained
   asymptotic 11% evidence bindings.
2. `python -B paper/check_finite_integration.py`: 746 authenticated files,
   seven selected certificate targets, four matching timing cells, 57 map
   words, and exact unions. Its checks also cover the 130-geometry Q1 study,
   110-cell length study, and five matching full-certificate anchors.
3. `python -B paper/build_application_tables.py --check`: application tables
   and exact conditional parameter bounds.

The existing PDF text has no `??`, TODO, FIXME, or TBD markers. Rendered pages
1, 26, and 27 were inspected: the front matter and main encoder/PCS tables are
legible, with no observed clipping or overlap. The existing TeX log includes
an amsmath accent warning and underfull boxes, so the old blanket “no warnings”
review records should not be treated as validation of this PDF.

No manuscript source was changed, no PDF was rebuilt, no benchmark was run,
and no archive was uploaded. Existing working-tree edits were preserved.

## Recommended next turn

Perform the final whole-document ePrint check and record the final source
revision and PDF checksum before upload. A separate deeper review can focus on the analytic
steps that the transcription checks do not establish: transfer domination,
the fixed-occupation continuum comparison, and uniform asymptotic remainders.
