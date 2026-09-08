# SPIN paper artifact

This is the reader-facing entry point for the paper's code and numerical
evidence. Use the [paper-to-code map](PAPER_MAP.md) to locate a specific result,
or the [reproduction guide](REPRODUCING.md) to run the checks.
The [validation record](VALIDATION.md) lists tests performed and the local
evidence archive's checksum.

## Reproduction levels

| Command, from the repository root | What it does | Requirements |
|---|---|---|
| `python -B artifact/reproduce.py quick` | Checks paper values, inner spectra, all 130 parameter-study geometries, and 31 asymptotic manifest entries. | Python 3.11+, standard library only. |
| `python -B artifact/reproduce.py figures` | Regenerates the three parameter-slice PGFPlots sources, then checks them. | Same. |
| `python -B artifact/reproduce.py paper` | Runs quick checks and builds the PDF. | Same, plus TeX Live/latexmk. |
| `python -B artifact/reproduce.py inventory` | Lists missing and hash-mismatched BCH audit dependencies, with byte totals. | Python and Git. |
| `python -B artifact/reproduce.py evidence` | Requires every retained BCH audit dependency locally and authenticates exact retained sums. | Complete local evidence plus Git; no numerical replay. |
| `python -B artifact/reproduce.py pack-evidence --output output/artifact/bch-evidence.zip` | Packages the matching BCH evidence and rechecks every archived SHA-256. Refuses to overwrite an existing ZIP. | Complete local evidence plus Git. |

These commands stop on failure. A hash check establishes that the retained
file matches its manifest, not that every mathematical reduction is correct.

## Artifact layout

```text
artifact/
  README.md              start here
  PAPER_MAP.md           paper claim -> source, data, checker
  REPRODUCING.md         setup, commands, scope, and release gaps
  reproduce.py           small command-line entry point
  test_reproduce.py      regression tests for artifact checks
  data/asymptotic/       compact relocated frozen dependency
```

Manuscript sources, numerical producers, and encoder kernels remain in their
documented directories. The file under `data/asymptotic/` preserves the bytes
of a small selected-map input whose historical path was inside an ignored
receipt tree. The original manifest still supplies its hash.

## What is ready, and what remains

The compact artifact supports checking reported values and rebuilding plots
and the paper. The small-BCH plots are regenerated from the tracked, rounded
tables, not from fresh full-grid computations. The BCH-256 exact ledgers are
separate from those diagnostic data.

A complete external numerical replay release still needs publication of the larger BCH
worker evidence and restoration of the original small-BCH grid inputs. The inventory command
makes missing BCH dependencies explicit. It does not fetch or manufacture
them. Publication must also pin the final source revision and record the
archive checksum. See [release preparation](REPRODUCING.md#release-preparation).
