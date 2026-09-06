# Finite-SPIN landscape database

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

The current snapshot contains 1,272 observations from 20 source files:
414 historical observations, 715 preferred activation-aware Q1 screens, and
143 superseded coarse-grid screens. The preferred pilot covers four exact
BCH spectra, four exact RM spectra, and three random-outer references at
five message lengths and 13 RM2Sub configurations. See
[`ACTIVATION_PARAMETER_STUDY.md`](ACTIVATION_PARAMETER_STUDY.md).

Schema version 2 preserves the historical `result_class` and adds a separate
`transfer_review_status`. The seven historical RM certificate entries are
`activation_under_reaudit`. Other historical entries are
`historical_pending_review`. The new screens are `activation_aware` but remain
binary64 diagnostics or ensemble references. Exact outer spectra do not
resolve an issue in the inner transfer.

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

The checked-in `spin_landscape.sqlite3` and `landscape_export.csv` are
reproducible convenience artifacts. Rebuild them after changing `catalog.json`
or any imported source.

Arbitrary read-only SQL is also available:

```text
python query_landscape.py sql "SELECT * FROM certified_results"
python query_landscape.py sql "SELECT * FROM landscape WHERE comparison_eligible=1"
python summarize_activation_pilot.py
```

Run the tests with:

```text
python -m unittest -v test_landscape_db.py
python -m unittest -v test_extrapolate_parameters.py
python -m unittest -v test_activation_q1.py
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
`comparison_eligible=1` selects the preferred pilot grid. The 143 original
`e=24` observations remain indexed but are superseded by their expanded-grid
recomputations. A preferred Q1 observation is not a full-distance certificate.

The completed BCH `[256,128]` RM2Sub work should enter through a new explicit
catalog entry. Its importer must identify whether each spectrum input is
exact, certified, estimated, or an ensemble reference. RandomStepConv results
are intentionally outside this database.
