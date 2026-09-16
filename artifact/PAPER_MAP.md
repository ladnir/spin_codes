# Paper-to-code map

Section numbers refer to the current draft. LaTeX labels are included so this
map can survive renumbering. Repository paths are authoritative; older notes
may contain historical worktree paths.

| Paper result | Source and retained evidence | Reproduction entry point |
|---|---|---|
| Common first-moment framework, Section 3 | [framework.tex](../paper/framework.tex) | Analytic argument; not a numerical benchmark. |
| Accumulator and Random SPIN, Sections 4--5 | [accumulator_spin.tex](../paper/accumulator_spin.tex), [random_spin.tex](../paper/random_spin.tex), [proof_appendix.tex](../paper/proof_appendix.tex) | Analytic proofs in the manuscript. |
| Scalable Structured SPIN, Theorem 6.3 (`thm:structured-spin-scalable`), Appendix A | [certificate snapshot](../workstreams/paper_architecture/certificates/single_sampled_ba_rm2sub/README.md), [manifest](../workstreams/paper_architecture/certificates/single_sampled_ba_rm2sub/SINGLE_SAMPLED_BA_RM2SUB_CERTIFICATE_MANIFEST.json) | `reproduce.py quick` checks the 31 pinned entries. The snapshot README lists interval producers; see the replay boundary below. |
| Parameter slices, Section 7.1, Figures 1--3 (`fig:finite-k-b`, `fig:finite-s-t`, `fig:finite-k-s`) | [retained tables](../workstreams/finite_asymptotic_theory/landscape_db/CURRENT_RESULTS_TABLES.md), [interpretation](../workstreams/finite_asymptotic_theory/landscape_db/BCH_DOMINANCE_ANALYSIS.md), [figure generator](../paper/build_parameter_figures.py) | `reproduce.py figures`; 130 geometries, 77 useful full bounds, 7 weak bounds, 46 first-moment obstructions. No fresh grid evaluation. |
| BCH-256, Theorem 7.1, Table 2, Figure 4 (`thm:finite-bch-spin`, `tab:finite-bch-margins`, `fig:finite-bch-curve`), Appendix B | [finite handoff](../workstreams/bch_rm2sub_bridge/PAPER_HANDOFF.md), [ledger directory](../workstreams/bch_rm2sub_bridge/generated/), [paper integration checker](../paper/check_finite_integration.py) | `reproduce.py quick` checks exact-ledger margins and the selected inner. `reproduce.py evidence` authenticates the larger retained evidence set. Neither reruns all interval calculations. |
| Performance, Section 9, Table 3 (`tab:finite-bch-performance`) | [encoder and build guide](../workstreams/bare_bch_rm2sub/README.md), [measurements](../workstreams/bare_bch_rm2sub/PERFORMANCE.json), [methodology](../workstreams/bare_bch_rm2sub/PERFORMANCE.md) | Build and run CTest for correctness; run the benchmark separately and serially for new timings. |

The additional performance comparison (`sec:transposed-comparison`,
`tab:transposed-comparison`) is reproduced by the
[comparison guide](../workstreams/transposed_comparison/README.md).
Its [no-reset observations](../workstreams/transposed_comparison/results_no_reset_20260911.json)
contain 54 serial runs. Run `python -B workstreams/transposed_comparison/report.py --check`
to check the 18 table cells without benchmarking; this separate check is
not currently part of `reproduce.py quick`.

## The five finite ledgers

All paths below are within `workstreams/bch_rm2sub_bridge/generated/`.

| log2 K | Ledger | Full margin (display only) |
|---:|---|---:|
| 16 | [t128_s19_m16_full_split_coverage_v1.json](../workstreams/bch_rm2sub_bridge/generated/t128_s19_m16_full_split_coverage_v1.json) | 53.944367 |
| 18 | [t128_s19_m18_full_split_coverage_v1.json](../workstreams/bch_rm2sub_bridge/generated/t128_s19_m18_full_split_coverage_v1.json) | 52.346388 |
| 20 | [t128_s19_m20_ladder_full_v1.json](../workstreams/bch_rm2sub_bridge/generated/t128_s19_m20_ladder_full_v1.json) | 50.448203 |
| 22 | [t128_s19_m22_ladder_full_v1.json](../workstreams/bch_rm2sub_bridge/generated/t128_s19_m22_ladder_full_v1.json) | 48.470684 |
| 24 | [t128_s19_m24_ladder_full_v1.json](../workstreams/bch_rm2sub_bridge/generated/t128_s19_m24_ladder_full_v1.json) | 46.476222 |

The theorem uses rational upper bounds, not these rounded decimals. These
five lengths share the selected (t,s)=(128,19) map. The small-BCH plots use
different, nested calibration maps; they are not additional matched points
on the BCH-256 curve.

## Replay boundary

The asymptotic manifest preserves 31 accepted sources and receipts. A fresh
numerical replay invokes the producers and checks their output, in addition
to reading the analytic proof. The finite BCH source checkout retains compact
ledgers and the completion audit, but not every worker input or receipt.
The [reproduction guide](REPRODUCING.md) gives the separate procedures.
