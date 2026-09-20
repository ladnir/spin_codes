# Practical BCH-256 parameter search

Current proof-search policy: a weak pooled endpoint screen must be followed
by the explicit all-one-row split before spending more effort on that endpoint.
Use `search_split_certificates.py` and `fine_split_search.py` for the selected
t128_s19 candidate. See [T128_S19_M16_CLOSURE.md](T128_S19_M16_CLOSURE.md).
The catalog screener below remains a proposal tool; its frozen source and
historical receipts are unchanged.

The purpose of this workflow is to choose which inner parameters to investigate
or benchmark next. A successful run produces a ranked grid with evidence labels.
It does not require every configuration to receive a distance certificate.

`parameter_grid.py` sweeps the locally audited map catalog over message sizes.
The initial catalog is t64_s16, t64_s20, t128_s15, t128_s19, and t256_s14.
Each name identifies a snapshotted map, not just nominal t and s values.
The fixed BCH-256 outer, cutoff, setup randomness, and input hashes are retained.
This version does not generate new maps or enumerate arbitrary state dimensions.

## Run the search

From the bridge directory:

```powershell
Set-Location C:/Users/peter/repo/permute_conv-github-bch/workstreams/bch_rm2sub_bridge
C:\Python314\python.exe -B parameter_grid.py --directory generated/my_grid --m 16 18 20 --target-bits 40
```

The command runs the five maps at K=2^16, 2^18, and 2^20. It always includes
t64_s20 as a matched reference for each size. To choose a subset, pass
`--configurations t64_s16 t128_s19`. Use `--list-configurations` to list the
available selected maps. Message exponents 13 through 20 are supported when
the row count contains complete epochs; unsupported geometry is reported.

The required numerical dependencies and ignored audit inputs must already be
installed/restored, as for the existing bridge scripts. The search authenticates
the selected maps and frozen BCH inputs. It does not fetch another worktree's
files or reinterpret a missing map as a different selection.

## Two stages, with a stopping rule

First, screen every requested cell. The worker computes:

- The existing Q1 upper bound at four fixed tilts, with a 256-bit Arb producer.
  This grid screen does not spend time replaying Q1 separately for every cell.
- The existing four-state mixed-density diagnostic at four occupancies:
  approximately L/16, L/4, L/2, and L, where L=K/128 is the number of outer rows.
  Defaults are 17 density points, 13 tilts, and three existing row-probability banks.
- The exact map's minimum A and kernel weights, and the number of epoch updates.
- The established full certificate for t64_s20 at K=2^20, if local legacy reuse
  is enabled. No other cell inherits that certificate.

Second, spend a limited refinement budget on selected candidates. By default,
the driver chooses the top two shortlisted cells, starting at the largest K.
It probes Q2 and approximately L/3 with the existing fixed-weight engine.
Each occupancy gets at most two automatically predicted neighboring tilts.
A useful bound is replayed at 512 bits. A passing anchor does not cover its
neighbors, and unresolved probes do not trigger an unlimited closure search.

The endpoint is deliberately not the default expensive refinement target.
The preceding investigation showed that repeated endpoint tuning can consume
the budget without helping parameter selection. Full coverage remains available
through the separate `search_certificates.py` workflow when a selected finalist
needs it. The grid does not invoke that full-coverage workflow automatically.

## How ranking works

At each K, the t64_s20 screen supplies a reference computed with the same settings.
Its sampled bound may be weak even when a full certificate is known. Therefore
the driver compares screen values with that reference rather than requiring the
sampled mixed bound itself to reach 40 bits.

A candidate enters `SHORTLIST_HEURISTIC` when both default search preferences hold:

1. Its Q1 producer margin is at least the target plus two bits.
2. Its worst sampled mixed margin is no more than 0.05 per output bit below
   the reference's worst sampled mixed margin.

The second quantity is `(reference margin - candidate margin)/(2K)`.
The threshold 0.05 is an exposed exploration tolerance, not a theorem or a fitted
error bar. Change it with `--reference-tolerance`; change the Q1 preference
with `--q1-headroom`. A candidate outside either preference is `WATCH_HEURISTIC`,
not rejected as insecure or mathematically impossible.

Within each class, fewer epoch updates rank first, then fewer state bits.
The update count is 2K/t. It is an operation-count proxy, not measured runtime:
different state arithmetic and implementations can change the actual ranking.
Remaining ties use Q1 and sampled-screen values. The JSON also marks the
Pareto set for these search measurements; that flag has no security meaning.

Cells with incomplete screens or a missing matched reference receive no rank.
Every completed candidate keeps its numerical observations and evidence level,
including candidates outside the shortlist. Optional anchor results are reported
alongside the screening rank; a weak upper bound does not establish an obstruction.

## Reading the output

Each run writes a numbered JSON summary and a Markdown table under its output
directory. Important fields are:

| Field | Meaning |
|---|---|
| `priority`, `rank` | Heuristic exploration order at this K |
| `q1_margin_bits` | One-occupancy producer bound, not a full margin |
| `reference_deficit_per_output_bit` | Matched-reference screen comparison |
| `epoch_updates` | Operation-count proxy, not timing |
| `replayed_occupancies` | Individually checked anchors only |
| `refinement_signal` | Whether tested anchors passed, remained unresolved, or have not been tried |
| `full_margin_bits` | Full margin only when supported by the reused complete certificate; otherwise null |
| `screen_state` | Completed, pending, unsupported, timed out, or process error |

`screening_complete` means every requested cell has a terminal screening state.
Check `successful_screens` as well: timeouts and errors are terminal states,
but are not successful screens. No report field converts missing evidence into
a margin estimate.

The standalone report check reconstructs the latest table from cached jobs,
reauthenticates their dependencies, and checks replay receipts for accepted anchors:

```powershell
C:\Python314\python.exe -B parameter_grid.py --verify-report generated/my_grid/summary_0000.json
```

Use the latest numbered summary. Earlier snapshots describe an earlier cache
state and may no longer match after more jobs finish. This checker verifies
the report and provenance, not a security certificate for the grid.

## Budget and resume

Defaults are a 300-second work budget, at most 24 jobs, 30 seconds per screen,
and 45 seconds per refinement job. Input audit and reporting add some overhead.
Workers run serially and do not start grandchildren. A timed-out worker is
terminated and reaped before another starts. Do not launch concurrent searches
or benchmarks on the same resources.

Examples:

```powershell
# Cheap grid only; no anchor refinement.
C:\Python314\python.exe -B parameter_grid.py --directory generated/my_grid --refine-top 0

# Resume the same grid, adding two shortlisted cells and one uncertain cell.
C:\Python314\python.exe -B parameter_grid.py --directory generated/my_grid --refine-top 2 --refine-uncertain 1 --seconds 240
```

Repeating a command reuses completed jobs. Timeouts, errors, and interrupted jobs
are recorded and skipped rather than retried indefinitely. Use a fresh directory
to retry them with a larger per-job budget. Work budgets and refinement quotas
can change on resume; the map grid, screen resolution, and ranking preferences
cannot. A changed policy or selected-map fingerprint requires a new directory.

This is intentionally a small, explicit workflow. It has no automatic map
generation, runtime benchmark harness, adaptive global optimizer, or guarantee
that its shortlist contains the best configuration. Those are not prerequisites
for using it to direct the next experiments.

Only source, tests, and interpretation belong in Git. The per-job numerical
outputs and logs remain in ignored directories.

## Demonstrated grid

The local run is `generated/parameter_grid_v1/`. It screened all five selected
maps at message exponents 16, 18, and 20: 15 successful cells. The first invocation
was deliberately stopped after three jobs; the next invocation reused them.

With the default headroom and reference tolerance, all three K sizes have the
same exploration order:

| Class | First | Second |
|---|---|---|
| Shortlist | t128_s19 | t64_s16 |
| Watch | t256_s14 | t128_s15 |

t64_s20 remains the matched reference, outside the candidate ranking.
At K=2^20, the normalized sampled deficits are approximately 0.03505 for
t128_s19 and 0.01520 for t64_s16. Both meet the default 0.05 exploration tolerance;
t128_s19 ranks first because it uses half as many epoch updates. This is not
a measured speedup or a conclusion that either candidate meets the full margin.

Selective refinement at K=2^20 produced replayed Q2 and Q2730 anchors for both
shortlisted maps. These are new individually verified points, not coverage
between them. An additional watch-list cell, t256_s14, was probed separately.
Its positive Q2 result does not erase its weak screening evidence or establish
dense coverage.

For t128_s19 at Q2730, the first predicted tilt (-10 in the integer CLI
convention) gave a weak bound. The automatic neighboring tilt (-5) produced
a passing bound and completed 512-bit replay. For t256_s14, the first Q2730
bound was weak and its second attempt timed out at the remaining run budget.
Neither result was converted into a failure claim about that code.

The resumed search executed 20 jobs in approximately 240 seconds after the
initial three screening jobs. A further invocation executed zero jobs. Its
final table and data are `summary_0002.md` and `summary_0002.json`; the standalone
report checker verified the latter against the cached evidence.

The only reused full margin in this grid is t64_s20 at K=2^20: 50.487298 bits.
All other full-margin fields remain null. Thirty-eight tests passed, covering
ranking, reference handling, budgeted resume, report tampering, evidence labels,
and the existing arithmetic/receipt regressions.

The next practical use is to benchmark the shortlisted implementations serially
or add another audited map near their parameter frontier. Pursue full coverage
only for a finalist that needs it. The grid workflow itself is not waiting on
the mixed-weight endpoint research.
