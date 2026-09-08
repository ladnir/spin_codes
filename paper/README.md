# New SPIN manuscript

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

The revised draft builds with TeX Live 2026 and writes its PDF to
`../output/pdf/spin_codes_draft.pdf`. The September 7 revision adds the finite
BCH-256 theorem, engineering curve, implementation results, and finite proof
appendix. The archived sources above remain unchanged.
