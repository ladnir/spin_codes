# Window screens for two candidate geometries

The bounded searches found no dimension obstruction for either GRS[33,22]
or GRS[64,44] on route seeds zero, one, and two. The largest dimensions found
were 91 and 78, respectively, below their exit-constraint bounds of 132 and
256. These are search outcomes, not distance certificates or upper bounds on
the best possible window objective.

The GRS[64,44] screen used positive-dimensional anchored starts and relaxed
parity continuation. It therefore does not mistake a zero-objective local
plateau for evidence that the geometry is sound. The earlier GRS[32,22]
counterexample fixtures remain unchanged.

## Geometry and exact rejection criterion

The two candidates use four state symbols and packet width sixteen. Their
source rows are independently routed within each outer-coordinate region.
For one aligned contiguous window per region, let \(c_i\) be the number of
selected outer coordinates in source row \(i\).

For an outer GRS[\(b,k\)] code, write \(p=b-k\). Its shortened source dimension
inside the selected windows is

\[
 D_p=\sum_i\max(0,c_i-p).
\]

The [window obstruction](WINDOW_DIMENSION_SEARCH.md) applies when
\(D_p>4b\) and the total window length is at most
\(e=\lceil N/10\rceil-1\). It then gives a nonzero low-output word for every
fixed assignment of the four-state inner coefficients. The word can depend
on those coefficients. The criterion is sufficient for rejection, not
necessary.

| Outer | Rate | Rows \(L\) | Output \(N=bL\) | Initial width | Total window length | Strict limit \(e\) | Exit constraints |
|---|---:|---:|---:|---:|---:|---:|---:|
| GRS[33,22] | \(2/3\) | 47680 | 1573440 | 4752 | 156816 | 157343 | 132 |
| GRS[64,44] | \(11/16\) | 23872 | 1527808 | 2384 | 152576 | 152780 | 256 |

Both region lengths and all window endpoints are packet-aligned. The initial
width is \(16\lfloor e/(16b)\rfloor\). Subsequent length exchanges preserve
the listed total length, keep all windows nonempty, and remain within the
strict output budget.

## Threshold-aware initialization

An arbitrary random start can have \(D_{20}=0\) for GRS[64,44]. If no row has
twenty hits, every one-region true-objective move has zero gain. Retaining
the current start on ties then leaves this plateau unchanged.

The new helper,
[`repaired_geometry_windows.py`](theory/repaired_geometry_windows.py),
constructs eight initial layouts for each route:

- two random aligned layouts;
- three layouts centered on one selected row's position in every region;
- three layouts using two, three, or four selected anchor rows.

A centered layout contains its anchor in all \(b\) regions. That row alone
therefore contributes \(b-p=k\) to the true shortened dimension. In particular,
a GRS[64,44] anchored start has \(D_{20}\ge44\), before any optimization.

For several anchors, each region chooses an aligned window maximizing
weighted anchor coverage. The weights are descending powers of two, so
covering the first anchor outweighs every combination that excludes it.
The first anchor remains present in all regions; the others resolve available
window choices. Ties prefer the first anchor's centered window. This weighted
coverage is only an initialization objective, not a source-dimension bound.

Every layout then receives exact coordinate ascent at the **true** redundancy
\(p\). The helper additionally applies a relaxed-parity ladder to two layouts:
the first random layout and the four-anchor layout. It ascends at each integer
redundancy from ten through \(p\), retaining the best encountered layout by
its independently reconstructed true dimension \(D_p\). A final true-parity
ascent polishes that retained layout.

The distinction matters numerically. For the GRS[64,44] seed-zero random
ladder, the parity-ten surrogate reached 3315 while the true dimension
\(D_{20}\) was zero. That value is not a counterexample. The same continuation
eventually reached true dimension 44. Only true dimensions enter the result
table below.

## Bounded search and results

Each ascent has a limit of twelve sweeps. The search retains the best true
dimension among the eight direct paths and two ladder paths. From that
layout, it makes 512 pairwise length-exchange proposals, using transfers of
16, 64, 256, or 1024 symbols. Each proposal receives at most twelve
coordinate-polishing sweeps and is accepted only when the true dimension
does not decrease. Ties can change the layout.

Routes use PCG64 seeds zero, one, and two. Initialization uses seed 193711;
the length-exchange seed is 9173. Per-path seeds are explicit in the helper.
The replay used NumPy 2.4.4. Each final count vector was independently
reconstructed by concatenating the selected row labels and applying
`bincount`, rather than trusting incremental search state.

| Outer | Route seed | True dimension before exchanges | True dimension after exchanges | After excluding padding rows | Required for rejection |
|---|---:|---:|---:|---:|---:|
| GRS[33,22] | 0 | 71 | 81 | 81 | \(>132\) |
| GRS[33,22] | 1 | 66 | 82 | 82 | \(>132\) |
| GRS[33,22] | 2 | 69 | 91 | 91 | \(>132\) |
| GRS[64,44] | 0 | 48 | 71 | 71 | \(>256\) |
| GRS[64,44] | 1 | 57 | 78 | 78 | \(>256\) |
| GRS[64,44] | 2 | 58 | 78 | 78 | \(>256\) |

For GRS[64,44], five of the six direct random starts remained at dimension
zero. The remaining random start began at one and ascended to 44. Anchored
starts began at 44 or 45. The winning pre-exchange layouts came from parity
ladders, so the reported search was not merely a replay of zero plateaus.

The padding check conservatively excludes the last eighteen complete rows
for GRS[33,22] and the last forty-one for GRS[64,44]. In both cases, the
retained complete rows use 1048564 source symbols, at most the nominal
\(2^{20}\). This check assumes row-major padding; a different restriction
requires its actual source-coordinate constraints. No witness was found
under either the full or restricted source dimensions.

## Verification and limitations

Seven tests pass in
[`repaired_geometry_windows_tests.py`](theory/repaired_geometry_windows_tests.py).
They verify aligned anchor containment, exact weighted anchor choices,
positive true dimension for 64-region anchored starts, strict budgets, and
padding counts. They also verify separation of surrogate and true objectives
and deterministic replay of a complete small bounded screen.

```sh
python -B research/workstreams/char2_field_spin/theory/repaired_geometry_windows_tests.py
python -B research/workstreams/char2_field_spin/theory/repaired_geometry_windows.py --geometries 33/22 64/44 --route-seeds 0 1 2 --exchanges 512 --max-sweeps 12
```

The helper reuses the exact scoring and length-exchange routines from the
earlier window search. It changes neither those routines nor the frozen
GRS[32,22] witnesses. No native implementation, benchmark, remote run, or
setup-failure probability claim is included.

These two repaired geometries survive this bounded necessary-condition
screen only. Neither has a whole-code 10% distance certificate. The next
useful step is a global rank/count argument for a chosen geometry, with
threshold-aware window searches retained as an independent obstruction check.
