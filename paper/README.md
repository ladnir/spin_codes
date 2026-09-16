# SPIN manuscript (LLNCS)

The PCS application revision adds ordinary-encoder measurements, a standalone
SPIN–Brakedown section, and the optimized Flock comparison. See
[the application evidence map](../artifact/APPLICATION_RESULTS.md) for retained
measurements, timing definitions, and the conditional security scope.
Run `python -B paper/build_application_tables.py --check` from the repository
root to check its tables without benchmarking.

The asymptotic construction now uses [Independent-Map Transvection (IMT)](../workstreams/inner_design/IMT.md)
at 11% relative distance, with the same shared Golay--BA-3 outer and 39/4
growth constant. Its reviewed argument is in `structured_appendix.tex` and
`structured_imt_appendix.tex`; the latter includes the exact two maps.
The selected finite results now use verified IMT certificates and matching
timings: five half-rate lengths, plus two quarter-rate distance thresholds
at K=2^20. See the [integration ledger](../workstreams/inner_design/finite_migration/PAPER_RESULTS.md).
The three broader parameter-slice plots now use the authenticated 130-cell
IMT Q1 grid. The selected BCH-256 curve compares Q1 and full certificates at
five exactly matching configurations. The diagnostic maps differ from those
selected maps; no full-grid certificate is claimed.
The manuscript presents IMT without a preceding-inner comparison. Historical
source and receipt names remain unchanged in the artifact; Reed--Muller
expansion mathematics and the external RM-based BAA baseline remain in scope.
No library default changes with this manuscript update.

Run `python -B paper/check_imt_integration.py` from the repository root to
check all 38 shared map words and the retained 11% evidence bindings.
The finite checker authenticates seven selected certificate targets, four
timing cells, and the shared and quarter-rate map tables. These checks need
the local generated evidence; a clean source checkout does not include it.
The asymptotic check is additional to the finite checker and is not
included in the historical `artifact/reproduce.py quick` command.
Use `python -B paper/build_imt_comparison.py --check` for the current
comparison table. The preceding comparison generator remains historical.
The [review record](../workstreams/inner_design/imt_asymptotic/d11/PAPER_REVIEW.md)
describes the analytic review and the preserved evidence boundary.

The paper uses Springer's unmodified LLNCS 2.26 class and `splncs04.bst`,
vendored here from the [official CTAN package](https://ctan.org/pkg/llncs)
(copyright Springer, CC BY 4.0). The
[Eurocrypt 2027 submission instructions](https://eurocrypt.iacr.org/2027/papersubmission.php)
require default LLNCS fonts, margins and spacing, visible page numbers, and
anonymous submissions. No keywords are included. The text and appendices have
not been shortened to meet the length limits. The non-submission draft includes
a table of contents after the abstract; submission mode omits it. The draft
date was removed to use the conference front matter.

## Compile from PowerShell on this machine

LaTeX is installed inside WSL Ubuntu. These commands can be run from any
PowerShell directory. Bash preserves the `latexmk` arguments across WSL.

Author version:

```powershell
wsl -d Ubuntu --cd /mnt/c/Users/stani/OneDrive/Documents/SPIN_CODES/permute_conv/paper -- bash -lc 'latexmk -pdf -interaction=nonstopmode -halt-on-error -outdir=../output/pdf -jobname=spin_codes_draft main.tex'
```

Anonymous submission version (same full text):

```powershell
wsl -d Ubuntu --cd /mnt/c/Users/stani/OneDrive/Documents/SPIN_CODES/permute_conv/paper -- bash -lc 'latexmk -pdf -interaction=nonstopmode -halt-on-error -outdir=../output/pdf -jobname=spin_codes_submission submission.tex'
```

The outputs are `../output/pdf/spin_codes_draft.pdf` and
`../output/pdf/spin_codes_submission.pdf`. Run the same command after edits;
`latexmk` automatically runs LaTeX and BibTeX as often as needed.

The `\ifsubmission` macro is declared in `main.tex`. Its default is
`\submissionfalse` (authors and institutions shown); change this to
`\submissiontrue` for anonymous review. The `submission.tex` entry point forces
submission mode without changing the default. Submission mode also clears the
PDF author metadata and renders artifact links as plain descriptions so the
identifying GitHub repository is not exposed. Prior papers remain cited in the
third person. Both Chosen-Block BAA (ePrint 2026/1903) and the original
Block-Accumulate paper (CRYPTO 2026; ePrint 2025/1828) are cited in the
introduction's construction overview and related-work discussion.

With a native Windows TeX installation on `PATH`, the equivalent commands are:

```powershell
Set-Location 'C:\Users\stani\OneDrive\Documents\SPIN_CODES\permute_conv\paper'
latexmk -pdf -interaction=nonstopmode -halt-on-error '-outdir=../output/pdf' '-jobname=spin_codes_draft' main.tex
latexmk -pdf -interaction=nonstopmode -halt-on-error '-outdir=../output/pdf' '-jobname=spin_codes_submission' submission.tex
```

## Reproduction checks

Start at the repository's [artifact guide](../artifact/README.md) for the
one-command checks and [paper-to-code map](../artifact/PAPER_MAP.md).
The root command `python -B artifact/reproduce.py paper` checks the compact
inputs and builds the draft using the instructions below.

The current revision plan is [REVISION_PLAN.md](REVISION_PLAN.md). It covers
the asymptotic/finite structured-SPIN narrative, BCH-256 certificates, and
measured implementation. The root-level restructure plan is historical.

The manuscript entry point is `main.tex`. Build from this directory with
TeX Live (including PGFPlots), BibTeX, and latexmk:

```text
latexmk -pdf -interaction=nonstopmode -halt-on-error "-outdir=../output/pdf" "-jobname=spin_codes_draft" main.tex
```

The pre-restart root-level TeX draft is preserved in
`../old/paper_draft_pre_spin_2026-08-31/`. Its `ARCHIVE_MANIFEST.md` records
the moved files, byte lengths, and SHA-256 hashes.

Check finite integration from the repository root:

```text
python -B paper/check_finite_integration.py
```

Reproduce or check the IMT Q1 slices and matched certificate curve from the retained
pinned Q1 grid and full-certificate receipts (no grid search or numerical proof replay):

```text
python -B paper/build_imt_parameter_figures.py
python -B paper/build_imt_parameter_figures.py --check
```

The generator authenticates the no-constant IMT Q1 grid and its replay receipt,
checks all 130 geometries and nested map identities, then writes native PGFPlots
inputs under `paper/figures/`. Each input records the Q1 producer's SHA-256.
The selected BCH-256 plot compares Q1 with the full certificate for the exact
same maps at five lengths. It does not certify the different smaller-outer
diagnostic maps. See the [reproduction guide](../artifact/REPRODUCING.md) for
the evidence-package boundary and separate replay instructions.

The LLNCS revision was also compiled with TeX Live 2023 in WSL Ubuntu and writes its PDF to
`../output/pdf/spin_codes_draft.pdf`. The September 7 revision adds the finite
BCH-256 theorem, engineering curve, implementation results, and finite proof
appendix. The archived sources above remain unchanged.
