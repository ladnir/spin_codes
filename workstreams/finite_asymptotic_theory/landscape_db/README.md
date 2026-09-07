# Finite-SPIN landscape database

Git contains the analysis code, documentation, and catalog configuration.
Generated grids, receipts, databases, exports, plots, and compressed snapshots
are local research outputs and are ignored. No experiment history or data
bundle is published. The current worktree retains its existing data unchanged.
On a fresh checkout, run the producers described in `COMPLETE_GRID_PLAN.md`
before data-dependent reports, integration tests, or audits. The catalog
describes the recorded studies; it does not itself supply their result files.
Run `python check_source_only_git.py` before pushing to check both the tip
and every new object reachable from the branch.

The active BCH-64/128 follow-up is in
[`BCH_DOMINANCE_ANALYSIS.md`](BCH_DOMINANCE_ANALYSIS.md). It tests whether
the Q1 surface describes the full bound, with explicit higher-occupation
coverage and aggregation costs. Run `study_bch_dominance_v1.py` for the
130-point Q2..4 pass, then `verify_bch_dominance_v1.py` and
`report_bch_dominance_v1.py`. Q5 and higher require separate tail evidence.
Run numerical producers sequentially.

`bch_zero_state_lower_v2.py` supplies a fast lower-bound screen on all
130 geometries. `bch_zero_state_exact_grid_v3.py` strengthens the 78
t/s settings at K=2^20 using exact integer kernel coefficients; replay
positive witnesses with `verify_bch_zero_state_v3.py`. The report
`report_bch_evidence_v3.py` joins the sparse grid, verified lower bounds,
and available full-reference replays. Its evidence labels distinguish
useful full bounds, weak complete upper bounds, first-moment obstructions,
and unresolved full tails. A coverage map displays continuous margin
losses at every useful full-bound setting, without interpolating missing
values. The earlier v2 report is retained for its historical snapshot.
Generated CSV, JSON, and figures remain ignored.

Full-reference replay uses `verify_bch_full_reference_v2.py`: BCH-128
uses the default dense cutoff Q257; BCH-64 uses `--block 64
--dense-minimum 1025`. The required interval receipts must exist first.
The sparse interval producers and dense-cover refiners retain their
own versioned checkpoints. Do not overwrite an existing seed or change
a producer whose hash is bound by a retained receipt.

The t128/s20 references use a smaller type cover from the exact BCH
spectrum. For each block in {64,128}, run these steps sequentially:

```powershell
python seed_bch_dense_v3.py --block 128 --step 128 --state 20 --minimum 257 --nodes 127
python close_bch_dense_v11.py --block 128 --step 128 --state 20 --minimum 257 --target-bits 55 --maximum-refinements 1000
python close_bch_sparse_tail_v2.py --block 128 --step 128 --state 20 --minimum 5 --maximum 256
python verify_bch_full_reference_v3.py --block 128 --step 128 --state 20
```

Replace `--block 128` by `--block 64` for BCH-64. These commands require
the original sparse grid and input maps. The seed refuses to overwrite
an existing file; the dense producer resumes its own saved checkpoint.
A complete dense cover may still give a weak bound. Only the full
replay establishes the margin after all occupation intervals are added.
The v3 replay checks the bulk cover and the character-based activation
transfer. It remains a binary64 diagnostic, with selected high-precision
checks, rather than an outward arithmetic certificate.

`complete_bch_reference_batch_v1.py` runs the same stages strictly
sequentially for a list such as `--geometry 64:128:19:20 --geometry
128:128:19:20`, where entries are B:T:S:log2(K). It reuses existing
intervals and always runs the complete verifier. This version requires
L>=257 and complete epochs. Its refinement limit controls search effort;
a weak complete upper bound is retained and labeled separately by the
report. Do not run another numerical producer while the batch is live.

The current engineering comparison is in
[`CONSTITUENT_ENGINEERING_SURFACES.md`](CONSTITUENT_ENGINEERING_SURFACES.md).
It treats exact BCH and RM, then the random ensemble, with the same
four-state Q1 transfer and measured K/s interaction slices. It explains
the different block-size trends and tests a long-region model with explicit
state cancellation. Its tables and figures remain local and ignored.

The preceding BCH-focused analysis is in
[`BCH_GROWTH_ANALYSIS.md`](BCH_GROWTH_ANALYSIS.md). It separates message-length
counting from state-size effects, explains the exact BCH-64/128 Q1 plateaus,
and documents a tighter four-state transfer. Its 212-row study and figures
are local outputs, separate from the complete-grid ledger. Follow its
sequential reproduction commands after generating the required pilot inputs.

This directory indexes the known finite results used to compare Structured
SPIN variants. The current database contains only RM2Sub inners. A `random`
row denotes a random outer-code reference; it never denotes RandomStepConv.

The database keeps four result classes separate:

- `certified`: an authenticated outward bound;
- `diagnostic`: a reproducible numerical calculation without outward rounding;
- `estimated`: a stated spectrum or performance estimate; and
- `reference`: an ensemble calculation used only for comparison.

The source receipt remains authoritative. Each database row records its source
path, SHA-256 hash, and row locator. The database does not convert a diagnostic
or estimate into a theorem.

The complete-grid study and its finite scope are documented in
[`COMPLETE_GRID_PLAN.md`](COMPLETE_GRID_PLAN.md). The grid has 3,108 native
parameter tuples. `complete_grid_coverage.json` and
`complete_grid_coverage.csv` report the completed and missing occupation
ranges; run `python register_complete_grid.py` to refresh this snapshot.
Q1..64, composition-preserving Q2/Q3/Q4, and a typed Q65..L cover have
been evaluated on every native tuple. `complete_grid_unions.csv` separates
coverage, useful positive diagnostics and certificates. Many dense bounds
remain loose. See [`GRID_FINDINGS.md`](GRID_FINDINGS.md) for measured
parameter comparisons and [`NEXT_GRID_REFINEMENTS.md`](NEXT_GRID_REFINEMENTS.md)
for the next work. Run `python audit_complete_landscape.py` to verify the
coverage, source hashes, union components, reports, and compressed snapshot.

The preceding pilot snapshot contained 1,605 observations from 22 source files:
414 historical observations, 1,012 preferred activation-aware Q1 screens,
36 activation-aware Q2 screens, and 143 superseded coarse-grid screens.
The preferred Q1 pilot covers four exact
BCH spectra, four exact RM spectra, and three random-outer references at
five message lengths and 22 RM2Sub configurations. The nine smaller-state
configurations currently cover message exponents 16,18,20. See
[`ACTIVATION_PARAMETER_STUDY.md`](ACTIVATION_PARAMETER_STUDY.md) and
[`SMALL_STATE_AND_Q2_STUDY.md`](SMALL_STATE_AND_Q2_STUDY.md).

Schema version 4 preserves the historical `result_class` and a separate
`transfer_review_status`. The seven historical RM certificate entries are
`activation_under_reaudit`. Other historical entries are
`historical_pending_review`. The new screens are `activation_aware` but remain
binary64 diagnostics or ensemble references. Exact outer spectra do not
resolve an issue in the inner transfer.

Conditional random bounds have model kind `random_setup_spectrum_caps`,
with explicit `setup_event_id` and `setup_failure_bits` fields. The shared
failure term must be added once when combining occupations. They are not
unconditional ensemble moments and are not products of expected spectra.
`witness_shift` records the auxiliary density parameter used in band bounds.
The balanced producers instead record `witness_probability_scale` in their
authenticated CSV and name the selected scale in each database row's notes.
The coverage and union reports check containment before reusing an older
conditional result under a stronger spectrum event.

`certified_results` excludes entries under review. It is currently empty:
the completed BCH-256 certificate is still in the owning task's worktree.
This database has not imported or independently replayed that certificate.

## Build and query

From this directory, run:

```text
python build_landscape_db.py
python query_landscape.py summary
python query_landscape.py curves --family bch
python query_landscape.py curves --family rm
python query_landscape.py curves --family random
```

Export the joined `landscape` view with:

```text
python query_landscape.py export landscape_export.csv
```

Regenerate the diagnostic parameter extrapolation with:

```text
python extrapolate_parameters.py
```

The output `parameter_extrapolation.json` fits only the exact-spectrum BCH
occupation-one rows under the matched-persistence RM2Sub schedule. It also
constructs RM2Sub capacity staircases. The output is an experiment-planning
artifact, not a distance certificate.
The old fit is now explicitly marked historical and pending activation review;
do not use its fit or its certificate-anchored schedule for current selection.

The local `spin_landscape.sqlite3.gz` and `landscape_export.csv.gz` are
ignored convenience snapshots. Run `python landscape_snapshot.py restore`
to obtain the ordinary local SQLite and CSV files, or rebuild and export
them using the commands above. The restore command verifies both compressed
and expanded hashes and preserves any different existing local file.
After changing `catalog.json` or an imported source, rebuild, export, and run
`python landscape_snapshot.py pack`. The uncompressed convenience files are
not tracked in Git. The compressed files and their generated manifest are
also local; regenerate them with `pack` after a rebuild and export.

`python parameter_cost_frontiers.py` exports tested state and epoch choices
at Q1, Q1..4, Q1..16, Q1..32, Q1..64, and full coverage. Its XOR counts and
state-update frequency are separate operation proxies. They exclude outer
encoding, field arithmetic, routing, and optimized circuit sharing, so they
do not establish runtime rankings across constituents.

Arbitrary read-only SQL is also available:

```text
python query_landscape.py sql "SELECT * FROM certified_results"
python query_landscape.py sql "SELECT * FROM landscape WHERE comparison_eligible=1"
python summarize_activation_pilot.py
python summarize_sparse_comparison.py
```

Run the tests with:

```text
python -m unittest -v test_landscape_db.py
python -m unittest -v test_extrapolate_parameters.py
python -m unittest -v test_activation_q1.py
powershell -ExecutionPolicy Bypass -File build_activation_q2.ps1
python -m unittest -v test_activation_q2.py
```

## Current coverage

The import catalog includes:

- exact-spectrum BCH and RM occupation-one curves under RM2Sub;
- random outer-spectrum references under the same RM2Sub configurations;
- RM2Sub epoch-length and state-size calibration results;
- the RM(4,9) occupation ladder; and
- the historical outward RM(4,9), `t=64`, `s=14` receipt, under re-audit; and
- the new activation-aware nested-map Q1 pilot.

Overlapping studies remain separate observations. The query layer does not
select the largest margin automatically because different studies can use
different witness grids or proof reductions.
`comparison_eligible=1` selects the preferred Q1 and Q2 observations. The 143 original
`e=24` observations remain indexed but are superseded by their expanded-grid
recomputations. A preferred Q1 observation is not a full-distance certificate.

The completed BCH `[256,128]` RM2Sub work should enter through a new explicit
catalog entry. Its importer must identify whether each spectrum input is
exact, certified, estimated, or an ensemble reference. RandomStepConv results
are intentionally outside this database.
