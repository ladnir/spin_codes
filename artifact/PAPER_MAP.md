# Paper-to-code map

Author-side research map, not the scope of the distributed artifact.
The [core-code guide](README.md) defines the artifact; numerical evidence and
comparison experiments are outside that scope.

LaTeX labels identify results across renumbering. GitHub is authoritative;
older notes may contain historical worktree paths.

## Manuscript layout

The main paper keeps the construction, central proof, finite guarantees,
transposed-encoder comparison, and standalone PCS comparison.
The main [Applications section](../paper/applications.tex) summarizes
PCG, PCS, and Flock results; detailed Flock and Bolt tables are supplementary.

- [Central structured proof](../paper/structured_proof.tex): main argument;
  the existing structured appendices retain its complete supporting estimates.
- [Finite appendix](../paper/finite_appendix.tex): spectrum bounds, proof
  details, and the complete half- and quarter-rate margin tables.
- [Engineering appendix](../paper/engineering_appendix.tex): all four
  parameter/length figures, cancellation explanation, and parameter-selection
  guidance. Diagnostic Q1 curves remain distinct from full certificates.
- [Implementation appendix](../paper/implementation_appendix.tex): complete
  timing tables, ordinary-encoding comparison, memory costs, and protocols.
- [Application appendix](../paper/applications_appendix.tex): measured Silent
  OT timings, full PCS and Flock methods, and labeled comparison projections.
  [OT performance](OT_PERFORMANCE.md) pins the per-party libOTe measurements
  and their reproduction contract.

All are included in both PDF builds. Figure/table generators and data
bindings are unchanged; moving a result does not change its evidence scope.

## Result-to-code index

| Paper result | Current sources | Check or reproduction |
|---|---|---|
| Common first-moment framework | [framework.tex](../paper/framework.tex) | Analytic argument. |
| Accumulator and Random SPIN | [accumulator_spin.tex](../paper/accumulator_spin.tex), [random_spin.tex](../paper/random_spin.tex), [proof_appendix.tex](../paper/proof_appendix.tex) | Analytic proofs. |
| Scalable Structured SPIN, IMT at 11% (`thm:structured-spin-scalable`) | [IMT proof guide](../workstreams/inner_design/imt_asymptotic/README.md), [11% refinement](../workstreams/inner_design/imt_asymptotic/d11/PROOF_UPDATE.md) | `python -B paper/check_imt_integration.py`; guide supplies numerical replay commands. |
| Shared BA outer tails (`lem:structured-ba-tails`) | [Appendix proof](../paper/structured_appendix.tex), [original outer certificate](../workstreams/paper_architecture/certificates/single_sampled_ba_rm2sub/linear_time_audit/LINEAR_OUTER_CERTIFICATE.md) | See the outer-tail replay commands below; exact finite regression checks are supplementary. |
| BCH-256 finite theorem and engineering curve (`thm:finite-bch-spin`, `fig:finite-bch-curve`) | [selected finite IMT results](../workstreams/inner_design/finite_migration/PAPER_RESULTS.md), [finite appendix](../paper/finite_appendix.tex) | `python -B paper/check_finite_integration.py` authenticates accepted evidence and checks exact unions, maps, and transcription. No interval replay. |
| Quarter-rate operating points (`thm:finite-quarter-spin`) | Same ledger; [independent-map transfer argument](../workstreams/inner_design/asymmetric/TRANSFER_ARGUMENT.md), [outer construction](../workstreams/rate_quarter_bch/SMALLER_OUTER.md) | Same checker; current margins come from IMT receipts, not the older RM2Sub outer study. |
| Selected transpose performance (`tab:finite-bch-performance`) | [timing and proof bindings](../workstreams/inner_design/finite_migration/PAPER_RESULTS.md) | Same checker; 3 half-rate sizes and 1 certified quarter-rate size. New benchmarks must run serially. |
| External comparison (`tab:transposed-comparison`) | [external campaign](../workstreams/transposed_comparison/README.md), [combined generator](../paper/build_imt_comparison.py) | `python -B paper/build_imt_comparison.py --check`; IMT replaces only SPIN, retaining external measurements. |
| IMT Q1 parameter slices (`fig:finite-k-b`, `fig:finite-s-t`, `fig:finite-k-s`) | [study and scope](../workstreams/inner_design/finite_migration/PARAMETER_SLICES.md), [generator](../paper/build_imt_parameter_figures.py) | `python -B paper/build_imt_parameter_figures.py --check` authenticates 130 Q1 cells and five selected matched Q1/full anchors. No full-grid claim. |

## IMT notation and adjoint check

The IMT step inputs, outputs, and state vectors use columns; the full
message/codeword generator convention remains row-oriented. Thus the local
map matrices have shapes `A: t x s` and `C: s x t`. The reverse recurrence
uses `C^T` for output feedback and `A^T` for state feedback.
Weighted transfer measures propagate as rows, with `e_Z` the zero-state unit
row in both the finite and asymptotic arguments. Collatz witnesses and the
terminal all-ones vector are columns. The occupation envelope is introduced
as the binomial mixture of transfer bounds for uniform fixed-weight inputs.

Run `python -B -m unittest discover -s paper -p test_imt_adjoint.py` to compare
the reverse recurrence against the transpose of an independently assembled
full forward matrix. The exact binary checks cover independent small maps
and the selected manuscript maps over one and three steps. They test the
algebra, not the native implementation or the distance certificate.
The 2026-09-17 check passed; the finite and asymptotic integration checks
also passed after these notation-only changes.

## Shared outer-tail replay

The appendix now contains the full sparse/dense tail argument, rather than
referring back to the main body's assertion. In particular, the upper-tail
reduction costs a factor of the constituent length `b`; the sparse sum is
`O(b^(-3/2))`, so this factor still leaves a vanishing bound. Dense upper
tails use reflection of the limiting accumulator exponent, with the finite
parity shifts handled by continuity. This does not assert symmetry of a
realized constituent's spectrum.

From the repository root, run the exact finite regression checks:

```text
python -B -m unittest discover -s paper -p test_ba_outer_tails.py
```

For the two numerical proof obligations, use Python with `mpmath` installed:

```text
cd workstreams/paper_architecture/certificates/single_sampled_ba_rm2sub/linear_time_audit
python -B verify_outer_sparse_constants.py
python -B verify_outer_interval.py --mode low-tail --max-boxes 100000 --max-depth 100 --progress-every 2000
```

The 2026-09-17 replay passed both checks at 60-digit interval precision.
It reproduced the sparse geometric base upper endpoint `0.6565638131135586`
and the dense upper endpoint `-7.686951035770771e-08`, with 6,747 processed
boxes and no unresolved boxes. These checks concern the BA outer only; they
are not a fresh replay of the complete IMT distance certificate.
The historical directory name does not change that scope.

## Application measurements

The ordinary-encoding, standalone PCS, and Flock tables are mapped in
[APPLICATION_RESULTS.md](APPLICATION_RESULTS.md). Their compact source data are
in `paper/data/application_results.json`; run
`python -B paper/build_application_tables.py --check` to check the tables.

The asterisked Bolt-max row in `tab:spin-pcs-standalone` is an amortized-limit
projection, not a complete measured prover. [BOLT_PCS_ESTIMATE.md](BOLT_PCS_ESTIMATE.md)
records its timing and communication assumptions. The arithmetic-only
[communication calculator](../paper/bolt_communication.py) supplies the proof-size
estimate; `python -B -m unittest discover -s paper -p 'test_bolt_communication.py'`
checks the Merkle expectation against exhaustive small-tree enumeration.

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

## Finite probability reductions

The finite appendix states the output-tilt domain, the triangular index range
of the sparse recurrence, the contribution covered by each density interval,
and the independence and Bessel integral used in the tilted density comparison.
These clarifications do not change a certificate or parameter.

Run `python -B -m unittest discover -s paper -p test_finite_probability.py`
for four exact small-instance checks: ordered region coefficients, shuffled
Bernoulli domination, the sparse envelope including the all-one label, and
the exponential-tilt likelihood ratio. These use rational arithmetic and
synthetic noncommuting transfers; they do not replay the selected IMT bounds.
Use `python -B paper/check_finite_integration.py` for the separate check of
seven selected certificate targets, map transcription, and exact unions.

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
