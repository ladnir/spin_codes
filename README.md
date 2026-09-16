# SPIN codes

Research code and artifacts for **SPIN Codes: Single-Permutation INterleaved
Codes**. GitHub is the authoritative workspace; Overleaf is an optional
paper-only mirror.

Start with the [artifact guide](artifact/README.md). It gives reproduction
commands, expected results, and the boundary between the compact repository
and the larger numerical evidence bundle.

## Quick start

From the repository root, with Python 3.11 or newer:

```sh
python -B artifact/reproduce.py quick
```

This checks manuscript numbers, selected inner maps and spectra, plot data,
and the historical asymptotic manifest. It needs Python's standard library,
Git history, and the pinned local evidence; it does not launch a search or
benchmark. The selected finite results now use IMT. Its separate 11%
asymptotic check is `python -B paper/check_imt_integration.py` (see the
[artifact guide](artifact/README.md) for dependencies). The broader
parameter plots now show IMT Q1 diagnostics, with a separate matched Q1/full
comparison at the five certified BCH-256 lengths. They do not certify the full grid.
To build the paper, additionally
install TeX Live with latexmk, BibTeX, PGFPlots, and placeins:

```sh
python -B artifact/reproduce.py paper
```

The output is `output/pdf/spin_codes_draft.pdf`.

## Where to look

| Directory | Purpose |
|---|---|
| [artifact/](artifact/README.md) | Reader-facing reproduction guide, checks, and [paper-to-code map](artifact/PAPER_MAP.md). |
| [paper/](paper/README.md) | Current LaTeX manuscript and reproducible vector figures. |
| [Finite IMT results](workstreams/inner_design/finite_migration/PAPER_RESULTS.md) | Current selected certificates, matching timings, implementation paths, and migration status. |
| [workstreams/bare_bch_rm2sub/](workstreams/bare_bch_rm2sub/README.md) | Shared implementation machinery and historical RM2Sub encoder. |
| [IMT inner](workstreams/inner_design/IMT.md) | New inner: terminology, certified instances, and source-name mapping; the supported default remains RM2Sub. |
| [workstreams/transposed_comparison/](workstreams/transposed_comparison/README.md) | Same-host transposed SPIN, chosen BAA, Expand--Convolute, and binary RAA comparison. |
| [workstreams/bch_rm2sub_bridge/](workstreams/bch_rm2sub_bridge/PAPER_HANDOFF.md) | Finite BCH-256 proof sources, compact ledgers, and numerical replay notes. |
| [workstreams/paper_architecture/certificates/single_sampled_ba_rm2sub/](workstreams/paper_architecture/certificates/single_sampled_ba_rm2sub/README.md) | Frozen asymptotic Structured SPIN certificate snapshot. |
| [workstreams/finite_asymptotic_theory/landscape_db/](workstreams/finite_asymptotic_theory/landscape_db/CURRENT_ENGINEERING_RESULTS.md) | Small-constituent parameter study and retained numerical tables. |
| `constructions/`, other `workstreams/`, `explorations/` | Supporting implementations and research history; not the starting point for paper reproduction. |
| `BA_paper/`, `enumerator_paper/`, `expander_codes/` | Related research material. |
| `output/`, `out/`, `tmp/` | Local generated PDFs, builds, and scratch results; not release inputs. |

Frozen source paths are retained because manifests bind their contents and
scripts import them directly. The artifact guide supplies a stable navigation
layer without relocating or duplicating the research implementations.

## Reproducibility and publication

The quick checks verify retained results; they do not independently rerun
every numerical proof. Full BCH evidence is inventoried by
`python -B artifact/reproduce.py inventory`. The [reproduction guide](artifact/REPRODUCING.md)
separates figure regeneration, numerical replay, encoder correctness, and
benchmarking. Run benchmarks serially, never concurrently.

See [the publication policy](GITHUB_PUBLISH_POLICY.md) before adding data.
Commit source, compact selected results, and manifests, not experiment trees.
Before submission, publish a reviewed evidence archive and record an immutable
Git revision in the paper; no release archive has yet been published by this
artifact pass.
