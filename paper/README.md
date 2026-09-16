# SPIN manuscript (LLNCS)

The PCS application revision adds ordinary-encoder measurements, a standalone
SPIN–Brakedown section, and the optimized Flock comparison. See
[the application evidence map](../artifact/APPLICATION_RESULTS.md) for retained
measurements, timing definitions, and the conditional security scope.
Run `python -B paper/build_application_tables.py --check` from the repository
root to check its tables without benchmarking.

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
python -B workstreams/bch_rm2sub_bridge/verify_paper_milestone.py --require-local-evidence
```

Reproduce or check the three parameter-slice figures from the retained
rounded tables (no grid search or numerical proof replay):

```text
python -B paper/build_parameter_figures.py
python -B paper/build_parameter_figures.py --check
```

The generator checks all 130 BCH-64/128 geometries and their status counts,
then writes native PGFPlots inputs under `paper/figures/`. Each input records
the normalized SHA-256 of `CURRENT_RESULTS_TABLES.md`. Full margins are
retained to the table's six decimals; Q1 overlays are reconstructed only
where both the full margin and its loss are present. Missing/weak points
are not interpolated. These plots are diagnostic, unlike the separate
BCH-256 exact-ledger certificate curve.

The first checks table and plot values, generator transcription, the exact
inner spectra, and measured timings. The second authenticates the retained
proof evidence and exact final sums. Neither is a fresh interval-arithmetic
replay. See the bridge's `PAPER_HANDOFF.md` for full replay instructions and
the external artifact packaging boundary.

The LLNCS revision was also compiled with TeX Live 2023 in WSL Ubuntu and writes its PDF to
`../output/pdf/spin_codes_draft.pdf`. The September 7 revision adds the finite
BCH-256 theorem, engineering curve, implementation results, and finite proof
appendix. The archived sources above remain unchanged.
