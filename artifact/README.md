# SPIN paper artifact

This is the reader-facing entry point for the paper's code and numerical
evidence. Use the [paper-to-code map](PAPER_MAP.md) to locate a specific result,
or the [reproduction guide](REPRODUCING.md) to run the checks.
The [validation record](VALIDATION.md) lists tests performed and the local
evidence archive's checksum.

## Reproduction levels

The selected finite results now use IMT. Their current entry point is
`python -B paper/check_finite_integration.py`: seven certificate targets,
four timing cells, and exact map transcription. It needs the pinned local
evidence plus Python and Git. See the [current results ledger](../workstreams/inner_design/finite_migration/PAPER_RESULTS.md).
The three parameter plots now use IMT Q1 diagnostics. Five matched BCH-256
Q1/full comparisons provide certified operating points, not certificates
for the different diagnostic maps.

| Command, from the repository root | What it does | Requirements |
|---|---|---|
| `python -B artifact/imt_reproduce.py check` | Checks the selected finite IMT certificate, map, and timing bindings. | Python 3.11+, Git, pinned local evidence. |
| `python -B artifact/imt_reproduce.py figures` | Regenerates the adaptive-step length comparison, two state-size slices, and the matched Q1/full curve. | Same, including both Q1 studies and their replay receipts. |
| `python -B artifact/imt_reproduce.py inventory --include-q1` | Includes the state-grid and adaptive-length Q1 inputs in the selected finite evidence inventory; also supported by `pack`. | Complete matching local evidence. |
| `python -B artifact/imt_reproduce.py inventory` | Authenticates the selected finite IMT dependencies and reports missing files. Fails on absent or mismatched evidence. | Python and Git for tracked-file status. |
| `python -B artifact/imt_reproduce.py pack --output output/artifact/imt-evidence.zip` | Packages only those IMT pins and verifies every archived hash. Requires a new output path. | Complete matching local evidence. |
| `python -B artifact/reproduce.py quick` | Checks selected finite IMT values and bindings, historical parameter figures, and 31 imported asymptotic entries. | Python 3.11+, Git, pinned local evidence. |
| `python -B artifact/reproduce.py figures` | Regenerates historical RM2Sub figures, no longer used by the manuscript. | Same. |
| `python -B artifact/reproduce.py paper` | Runs quick checks and builds the PDF. | Same, plus TeX Live/latexmk. |
| `python -B artifact/reproduce.py inventory` | Lists historical BCH/RM2Sub audit dependencies, not a full IMT inventory. | Python and Git. |
| `python -B artifact/reproduce.py evidence` | Requires every retained BCH audit dependency locally and authenticates exact retained sums. | Complete local evidence plus Git; no numerical replay. |
| `python -B artifact/reproduce.py pack-evidence --output output/artifact/bch-evidence.zip` | Packages the matching BCH evidence and rechecks every archived SHA-256. Refuses to overwrite an existing ZIP. | Complete local evidence plus Git. |

The `evidence` and `pack-evidence` commands also retain their historical
BCH/RM2Sub scope; they do not package the current IMT proofs.
Use `python -B paper/build_imt_comparison.py --check` for the current
external comparison, whose SPIN row comes from the new IMT series.
These commands stop on failure. A hash check establishes that the retained
file matches its manifest, not that every mathematical reduction is correct.

The manuscript's new **IMT asymptotic theorem at 11%** has an additional check:
`python -B paper/check_imt_integration.py`. It checks the exact map table and
new evidence bindings and needs the local IMT records plus NumPy, SciPy,
mpmath, and python-flint. The historical `quick` command above does not cover
this addition. The [IMT guide](../workstreams/inner_design/imt_asymptotic/README.md)
gives the separate interval and exact-polynomial replay commands. Those local
records are not yet part of the external artifact release.

## Artifact layout

```text
artifact/
  README.md              start here
  PAPER_MAP.md           paper claim -> source, data, checker
  REPRODUCING.md         setup, commands, scope, and release gaps
  reproduce.py           small command-line entry point
  imt_reproduce.py       selected finite IMT check, inventory, and packaging
  test_reproduce.py      regression tests for artifact checks
  data/asymptotic/       compact relocated frozen dependency
```

Manuscript sources, numerical producers, and encoder kernels remain in their
documented directories. The file under `data/asymptotic/` preserves the bytes
of a small selected-map input whose historical path was inside an ignored
receipt tree. The original manifest still supplies its hash.

## What is ready, and what remains

The sources support checking current selected IMT values when the pinned
local receipts are present. They do not bundle those receipts. The three
historical small-BCH plots still regenerate from tracked rounded RM2Sub
tables, not from the new IMT grid. These are separate evidence sources.

A selected finite IMT inventory now authenticates 753 files (640,779,033
uncompressed bytes), with no missing or mismatched files locally. This is
not an inventory of the asymptotic IMT proof or unfinished parameter study.
A complete current release still needs separately supplied numerical receipts,
a fresh-environment replay test, and the finished full-occupancy parameter
figures. Publication must pin the final source revision and
record the external evidence archive checksum. No data files are added to
the source commit. See [release preparation](REPRODUCING.md#release-preparation).
