# Full selected (128,19) certificate at K=2^20

The selected BCH-256 / RM2Sub (t,s)=(128,19) construction now has a full
distance/setup certificate at K=2^20 with margin **50.448203312946134 bits**.
Every occupancy Q=1,...,8192 is covered. The exact rational full union is
less than 2^-40.

The code, selected maps, random setup, and event are those in
`CURRENT_UNDERSTANDING.md`. There are L=8192 outer rows, N=2097152 output
bits, and bad-output cutoff H=209715. State persists across regions, output
precedes update, and there is no final flush. This certificate bounds the
probability, over the shared setup, that any nonzero message has output
weight at most H. It is not the complete SPIN protocol's security level.

## Exact ledger and numerical evidence

Full ledger: `generated/t128_s19_m20_ladder_full_v1.json`.

SHA-256:
`657198cb67f805ec5fda4c9caa48e95063b67090badbb1eea5a38d6343ff03cd`.

| Contribution | Covered occupancies | Margin of its upper bound |
|---|---:|---:|
| Q1 | 1 | 50.44820424785166 bits |
| Sparse remainder | 2 through 511 | 71.00564656314114 bits |
| Dense tail | 512 through 8192 | 577.1425238434481 bits |
| Full union | 1 through 8192 | 50.448203312946134 bits |

The decimals summarize exact rational inequalities; they are not the
arithmetic used for acceptance. The three component bounds were separately
added back to the stored union, and the result was checked against 2^-40.
Q1 dominates this certificate. Further dense-tail tightening would not
materially improve the full margin.

The sparse report is
`generated/t128_s19_m20_ladder_sparse_v1/report_0001.json`. Every accepted
producer has a matching 512-bit numerical replay. The first bounded search
closed Q1..250; its last replay hit the run budget. A resumed run retried
that task and closed the remaining occupancies through Q511.

The dense certificate is
`generated/t128_s19_m20_ladder_dense_refined_v1.json`, with its exact replay
in the adjacent `_replay.json` file. Its source is
`generated/t128_s19_m20_ladder_fixed_reference_v2/cover_0000.json`.
All 1,011 rectangles passed a 512-bit replay recorded in
`replay_cover_0000.json` in that directory. The partition tree covers all
integer occupancies 512..8192 and every feasible label coordinate.

The source rectangles used the loose region-density factor (L+1)^256.
Their finite upper bounds cover the whole domain even though some fail
the source search's allocation. The final dense certificate multiplies
their complete union by the exact rational (364/8193)^256, using the
stronger comparison in `POISSON_DENSITY_REFINEMENT.md`. No source leaves
are omitted. The input/output tilt argument and the required composition
count are derived in `DENSE_TWO_TILT.md`.

`combine_ladder_dense.py` verifies sparse receipts, dense receipt provenance,
the exact partition, the density-factor conversion, complete occupancy
coverage, and the final rational union. It does not rerun the numerical
workers; those replays were executed before the ledger was constructed.

## Reproduction and next rung

Run commands from `workstreams/bch_rm2sub_bridge` in the authoritative
GitHub worktree. Existing receipts are write-once. Use new output paths
for another production run.

```text
python -B search_ladder_dense_fixed.py --directory generated/NEW_dense --m 20 --minimum 512 --seconds 180 --nodes 2000
python -B search_ladder_dense_fixed.py --directory generated/NEW_dense --m 20 --minimum 512 --verify
python -B refine_ladder_density_cover.py --source generated/NEW_dense/cover_0000.json --output generated/NEW_dense_refined.json
python -B refine_ladder_density_cover.py --output generated/NEW_dense_refined.json --verify
python -B search_sparse_only.py --directory generated/NEW_sparse --m 20 --last 511 --seconds 240 --max-jobs 40 --job-seconds 70
python -B combine_ladder_dense.py --sparse-report generated/NEW_sparse/report_0000.json --dense generated/NEW_dense_refined.json --output generated/NEW_full.json
```

These are bounded searches. Resume them if a budget expires, use the newest
numbered cover/report, and replay that cover before refining its density
factor. The full ledger refuses an occupancy gap or a union above target.

The combined regression run passed 63 tests. It includes exact toy changes
of measure, scalar interval bounds, partition validation, occupancy-gap
rejection, selected-map arithmetic, and the inherited search tests. Source
pins for the K16/K18 ledgers and all three ladder Q1 certificates remained
intact.

The subsequent K22 and K24 rungs are now complete; see
`T128_S19_M22_M24_CLOSURE.md` and `T128_S19_LADDER_TO_M24.md`.
The K20 result pins the new numerical sources and proof note; extend them
with new versions or wrappers rather than editing frozen inputs in place.
