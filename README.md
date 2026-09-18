# SPIN codes

Research code and artifacts for **SPIN Codes: Single-Permutation INterleaved
Codes**. GitHub is the authoritative workspace; Overleaf is an optional
paper-only mirror.

Start with the [core implementation guide](artifact/README.md). The artifact
covers encoder code, build instructions, and correctness tests, not the
research experiments or numerical proof archives.

## Quick start

Follow the [core build and test instructions](artifact/README.md#build-and-test-the-half-rate-encoder)
for the selected BCH/IMT implementation.

To build the paper, install TeX Live with latexmk, BibTeX, PGFPlots, and placeins,
then run from the repository root:

```sh
cd paper
latexmk -pdf -outdir=../output/pdf -jobname=spin_codes_draft main.tex
```

The output is `output/pdf/spin_codes_draft.pdf`. Compiling the committed
manuscript does not require regenerating its tables or replaying numerical proofs.

## Where to look

| Directory | Purpose |
|---|---|
| [artifact/](artifact/README.md) | Core implementation guide; also retains historical author-side scripts and notes. |
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
The artifact is limited to the core implementation. A complete experimental
or numerical-evidence archive is not a release requirement for that artifact.
Historical reproduction commands above remain author-side tools.
