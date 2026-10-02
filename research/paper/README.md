# SPIN manuscript (LLNCS)

The updated half-rate encoder tables use the optimized precomputed campaign:
0.340, 1.460, and 9.289 ms for transposed encoding at K=2^16, 2^18, and 2^20.
See [the performance ledger](PRECOMPUTED_PERFORMANCE.md) for source bindings,
timing policy, and the standalone-library regression kept separate from this
update. The K=2^16 row uses the additional (64,12,2) inner, described informally;
the original one-round certificate ladder is unchanged. Quarter-rate and standalone PCS timings are unchanged by that encoder update. The lead regular-noise Silent OT result is
2.944 ms per sender batch at K=2^18 (89.0 million hashed OTs/s), using the full
precomputed code. See [the regular-noise timing ledger](REGULAR_OT_PERFORMANCE.md).
The stationary heuristic-refresh comparison remains separate at 84.5 million
OTs/s; see [its timing ledger](STATIONARY_OT_PERFORMANCE.md).
`precomputed_results.py` authenticates
the six new forward/transpose cells from retained local measurements.

Unless stated otherwise, run commands from the repository root. The current
Windows checkout is `C:\Users\peter\repo\permute_conv-github-bch`.

The PCS application revision adds ordinary-encoder measurements, a standalone
SPIN–Brakedown section, and the optimized Flock comparison. See
[the application evidence map](../artifact/APPLICATION_RESULTS.md) for retained
measurements, timing definitions, and the conditional security scope.
Run `python -B research/paper/build_application_tables.py --check` from the repository
root to check its tables without benchmarking.

The asymptotic construction now uses [Independent-Map Transvection (IMT)](../workstreams/inner_design/IMT.md)
at 11% relative distance, with the same shared Golay--BA-3 outer and 39/4
growth constant. Its reviewed argument is in `structured_appendix.tex` and
`structured_imt_appendix.tex`; the latter includes the exact two maps.
The selected finite results now use verified IMT certificates and matching
timings: five half-rate lengths, plus two quarter-rate distance thresholds
at K=2^20. See the [integration ledger](../workstreams/inner_design/finite_migration/PAPER_RESULTS.md).
Figure 1 uses the authenticated 110-cell one-round IMT length study, with
BCH-256 central and smaller steps at short lengths. The two state-size plots
retain the earlier 130-cell Q1 grid. The selected BCH-256 curve compares Q1 and full certificates at
five exactly matching configurations. The diagnostic maps differ from those
selected maps; no full-grid certificate is claimed.
The manuscript presents IMT without a preceding-inner comparison. Historical
source and receipt names remain unchanged in the artifact; Reed--Muller
expansion mathematics and the external RM-based BAA baseline remain in scope.
No library default changes with this manuscript update.

The engineering explainer keeps one transvection per update. Figure 1 compares
fixed steps, length-dependent steps, and fixed-map full-refresh references.
The [length study](../workstreams/inner_design/finite_migration/ADAPTIVE_LENGTH.md)
derives a local cancellation probability, separates proof slack from encoder
changes, and gives numerical replay commands. Reproduce its plot with
`python -B research/paper/build_imt_length_figure.py`; add `--check` to authenticate it.
The earlier multi-round study is retained as supporting research, not an
active figure or an additional construction parameter in this section.

Run `python -B research/paper/check_imt_integration.py` from the repository root to
check all 38 shared map words and the retained 11% evidence bindings.
The finite checker authenticates the seven original certificate targets and
the additional two-round small-length certificate, seven timing cells, and
the selected maps. These checks need
the local generated evidence; a clean source checkout does not include it.
The asymptotic check is additional to the finite checker and is not
included in the historical `research/artifact/reproduce.py quick` command.
Use `python -B research/paper/build_imt_comparison.py --check` for the current
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

The September 29 Flock refresh uses the integrated compact-layout implementation:
94/351 ms for SPIN and 215/571 ms for Ligerito at 16k/65k compressions
(2.29x/1.63x throughput). The two-process measurements, proof sizes and
provenance are pinned in `data/flock_integrated_20260929.json` and imported
into `data/application_results.json`. Its `historical_flock` entry preserves
the preceding campaign. Regenerate with:

```powershell
python -B research/paper/build_application_tables.py --flock-summary research/paper/data/flock_integrated_20260929.json
```

The blog builder reads the same application data. The Bolt--Flock projection
uses the updated surrounding-prover residual while retaining its original
commitment measurements and opening calibration.

## Compile the paper

TeX Live 2026 is installed natively on this Windows machine. From PowerShell:

```powershell
Set-Location 'C:\Users\peter\repo\permute_conv-github-bch\research\paper'
latexmk -pdf -interaction=nonstopmode -halt-on-error '-outdir=../../output/pdf' '-jobname=spin_codes_draft' main.tex
```

This builds the author version for ePrint. To build the anonymous version from
that same directory, use:

```powershell
latexmk -pdf -interaction=nonstopmode -halt-on-error '-outdir=../../output/pdf' '-jobname=spin_codes_submission' submission.tex
```

For another checkout or Linux, start at the repository root and run:

```sh
cd research/paper
latexmk -pdf -interaction=nonstopmode -halt-on-error '-outdir=../../output/pdf' '-jobname=spin_codes_draft' main.tex
```

Both builds write to the repository's `output/pdf/` directory. Quote the complete
`-outdir=...` and `-jobname=...` arguments in PowerShell. LaTeX dependencies are
BibTeX, latexmk, PGFPlots, and placeins; building the paper does not require the
numerical evidence. Repeating the command rebuilds changed inputs as needed.

`main.tex` defaults to `\submissionfalse`, showing authors, affiliations,
repository links, the repository citation, and the contents list.
`submission.tex` selects anonymous mode without editing that default. It clears
PDF author metadata and suppresses identifying repository links and the citation.
Prior papers remain cited in the third person.

Return to the repository root before running the checks below:

```sh
cd ../..
```

## Reproduction checks

Start at the repository's [artifact guide](../artifact/README.md) for the
one-command checks and [paper-to-code map](../artifact/PAPER_MAP.md).
The historical wrapper `python -B research/artifact/reproduce.py paper` first
checks retained evidence, then builds under `research/output/pdf/`. It requires
local research inputs. Use the direct build above for the ePrint PDF under
`output/pdf/`; that build has no numerical-evidence dependency.

The current revision plan is [REVISION_PLAN.md](REVISION_PLAN.md). It covers
the asymptotic/finite structured-SPIN narrative, BCH-256 certificates, and
measured implementation. The root-level restructure plan is historical.

The pre-restart root-level TeX draft is preserved in
`../old/paper_draft_pre_spin_2026-08-31/`. Its `ARCHIVE_MANIFEST.md` records
the moved files, byte lengths, and SHA-256 hashes.

Check finite integration from the repository root:

```text
python -B research/paper/check_finite_integration.py
```

Reproduce or check the IMT Q1 slices and matched certificate curve from the retained
pinned Q1 grid and full-certificate receipts (no grid search or numerical proof replay):

```text
python -B research/paper/build_imt_parameter_figures.py
python -B research/paper/build_imt_parameter_figures.py --check
```

The generator authenticates both the no-constant state grid (130 geometries)
and the adaptive-length grid (110 cells), with their replay receipts and maps,
then writes native PGFPlots
inputs under `research/paper/figures/`. Each input records the Q1 producer's SHA-256.
The selected BCH-256 plot compares Q1 with the full certificate for the exact
same maps at five lengths. It does not certify the different smaller-outer
diagnostic maps. See the [reproduction guide](../artifact/REPRODUCING.md) for
the evidence-package boundary and separate replay instructions.

Earlier revisions were also compiled with TeX Live in WSL Ubuntu. Their old
checkout and output paths are historical; use the commands above for this tree.
The archived manuscript sources remain unchanged.
