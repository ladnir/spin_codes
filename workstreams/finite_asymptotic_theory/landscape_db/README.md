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

The initial snapshot contains 414 observations from 17 authenticated source
files. It indexes 16 outer models and 23 RM2Sub configurations. Seven rows are
components or the full union of the outward RM(4,9) certificate. The remaining
rows are binary64 diagnostics or random-outer references.

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

The checked-in `spin_landscape.sqlite3` and `landscape_export.csv` are
reproducible convenience artifacts. Rebuild them after changing `catalog.json`
or any imported source.

Arbitrary read-only SQL is also available:

```text
python query_landscape.py sql "SELECT * FROM certified_results"
```

Run the tests with:

```text
python -m unittest -v test_landscape_db.py
python -m unittest -v test_extrapolate_parameters.py
```

## Current coverage

The import catalog includes:

- exact-spectrum BCH and RM occupation-one curves under RM2Sub;
- random outer-spectrum references under the same RM2Sub configurations;
- RM2Sub epoch-length and state-size calibration results;
- the RM(4,9) occupation ladder; and
- the complete outward RM(4,9), `t=64`, `s=14` finite certificate.

Overlapping studies remain separate observations. The query layer does not
select the largest margin automatically because different studies can use
different witness grids or proof reductions.

The pending BCH `[256,128]` RM2Sub work should enter through a new explicit
catalog entry. Its importer must identify whether each spectrum input is
exact, certified, estimated, or an ensemble reference. RandomStepConv results
are intentionally outside this database.
