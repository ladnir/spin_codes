# Paper-to-code map

LaTeX labels identify results across renumbering. GitHub is authoritative;
older notes may contain historical worktree paths.

| Paper result | Current sources | Check or reproduction |
|---|---|---|
| Common first-moment framework | [framework.tex](../paper/framework.tex) | Analytic argument. |
| Accumulator and Random SPIN | [accumulator_spin.tex](../paper/accumulator_spin.tex), [random_spin.tex](../paper/random_spin.tex), [proof_appendix.tex](../paper/proof_appendix.tex) | Analytic proofs. |
| Scalable Structured SPIN, IMT at 11% (`thm:structured-spin-scalable`) | [IMT proof guide](../workstreams/inner_design/imt_asymptotic/README.md), [11% refinement](../workstreams/inner_design/imt_asymptotic/d11/PROOF_UPDATE.md) | `python -B paper/check_imt_integration.py`; guide supplies numerical replay commands. |
| BCH-256 finite theorem and engineering curve (`thm:finite-bch-spin`, `fig:finite-bch-curve`) | [selected finite IMT results](../workstreams/inner_design/finite_migration/PAPER_RESULTS.md), [finite appendix](../paper/finite_appendix.tex) | `python -B paper/check_finite_integration.py` authenticates accepted evidence and checks exact unions, maps, and transcription. No interval replay. |
| Quarter-rate operating points (`thm:finite-quarter-spin`) | Same ledger; [independent-map transfer argument](../workstreams/inner_design/asymmetric/TRANSFER_ARGUMENT.md), [outer construction](../workstreams/rate_quarter_bch/SMALLER_OUTER.md) | Same checker; current margins come from IMT receipts, not the older RM2Sub outer study. |
| Selected transpose performance (`tab:finite-bch-performance`) | [timing and proof bindings](../workstreams/inner_design/finite_migration/PAPER_RESULTS.md) | Same checker; 3 half-rate sizes and 1 certified quarter-rate size. New benchmarks must run serially. |
| External comparison (`tab:transposed-comparison`) | [external campaign](../workstreams/transposed_comparison/README.md), [combined generator](../paper/build_imt_comparison.py) | `python -B paper/build_imt_comparison.py --check`; IMT replaces only SPIN, retaining external measurements. |
| IMT Q1 parameter slices (`fig:finite-k-b`, `fig:finite-s-t`, `fig:finite-k-s`) | [study and scope](../workstreams/inner_design/finite_migration/PARAMETER_SLICES.md), [generator](../paper/build_imt_parameter_figures.py) | `python -B paper/build_imt_parameter_figures.py --check` authenticates 130 Q1 cells and five selected matched Q1/full anchors. No full-grid claim. |

## Selected finite IMT certificates

| Outer / rate | log2 K | Relative distance | Full margin (display only) |
|---|---:|---:|---:|
| BCH-256 / 1/2 | 16 | .10 | 41.818327 |
| BCH-256 / 1/2 | 18 | .10 | 50.189076 |
| BCH-256 / 1/2 | 20 | .10 | 50.062088 |
| BCH-256 / 1/2 | 22 | .10 | 48.390492 |
| BCH-256 / 1/2 | 24 | .10 | 46.457255 |
| BCH-128 / 1/4 | 20 | .165 | 41.048168 |
| BCH-128 / 1/4 | 20 | .19 | 30.033491 |

The half-rate lengths share both IMT maps with the asymptotic construction.
The quarter-rate code shares the expansion but uses weight-three feedback.
Theorems use exact rational upper bounds, not rounded display values.
The certified curve is not an estimate of the true BCH-256 spectrum.

## Evidence and release boundary

The IMT checker authenticates accepted receipt pins and their source inputs,
checks exact sums (including outward-rounded quarter-rate ceilings), and
compares the selected maps to the manuscript and implementation headers.
It is not a fresh numerical replay or an independent analytic proof review.

The external comparison retains its 2026-09-11 campaign. Two shared SPIN
source files have since changed; the combined generator authenticates their
measured bytes from commit `aafb3f59e7c3541e19b8623520c8a407ee2219e8`.
It does not ignore their hashes or substitute today's source. The new IMT
series was measured separately under the same host/compiler/protocol.

Generated IMT evidence remains local and uncommitted. A source-only checkout
must obtain or regenerate the pinned inputs before these checks can pass.
The older `reproduce.py inventory/evidence/pack-evidence` commands cover
the historical BCH/RM2Sub evidence set, not a complete current IMT release.
The [reproduction guide](REPRODUCING.md) retains those historical procedures.

Use `python -B artifact/imt_reproduce.py inventory` for the selected finite
IMT evidence instead. Its 753-file inventory authenticates the eight accepted
root receipts and their declared source pins. `pack --output <fresh.zip>`
packages that set and verifies its archived bytes. Neither command runs
interval arithmetic. Add `--include-q1` to include the diagnostic grid (766
files total rather than 753). The local
inventory is complete; no current IMT archive has been published.

The broader plots now follow the agreed Q1-focused scope. Their 130-cell grid
does not by itself establish full margins. The selected BCH-256 curve supplies
five matched Q1/full comparisons, not endpoint certificates for the different
diagnostic BCH-64/128 maps. Full-grid proof search is no longer a prerequisite
for this engineering explanation.
