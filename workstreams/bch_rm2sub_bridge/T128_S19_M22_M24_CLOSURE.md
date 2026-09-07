# Full selected (128,19) certificates at K=2^22 and K=2^24

Both remaining rungs have full distance/setup certificates above the 40-bit
target. They use freshly computed bounds for the same selected BCH-256 and
RM2Sub maps as the completed K20 rung.

| Message bits K | Outer rows L | Output bits N | Bad-weight cutoff H | Full certified margin |
|---:|---:|---:|---:|---:|
| 2^22 | 32768 | 8388608 | 838860 | 48.47068387175386 bits |
| 2^24 | 131072 | 33554432 | 3355443 | 46.476222182258425 bits |

The outer is the fixed p37-syndrome BCH [256,128] subcode with minimum
distance at least 38, not an arbitrary BCH code with these dimensions.
The selected inner has (t,s)=(128,19). The setup independently permutes
each outer row and each transposed region. It also supplies fresh independent
nonzero GF(2^19) multipliers. The state starts at zero, persists across
regions, and is updated after output. There is no final flush.
`CURRENT_UNDERSTANDING.md` gives the precise encoder recurrence.

Let E_theta be the resulting encoder for the shared setup theta. For each
row of the table, the ledger records an exact rational U satisfying

\[
\Pr_\theta[\exists x\ne0:\operatorname{wt}(E_\theta(x))\le H]
\le U < 2^{-40}.
\]

The displayed margin is -log2(U). This is a distance/setup statement, not
the complete SPIN protocol's security level. The proof does not assume
independence between codewords, an exact BCH-256 spectrum, or a heuristic
growth curve. It bounds the bad-message expectation at every occupancy
and sums those bounds.

## Contribution bounds and receipts

| Contribution | Occupancies | K=2^22 margin | K=2^24 margin |
|---|---:|---:|---:|
| Q1 | 1 | 48.4706841091504 | 46.476222245808 |
| Sparse remainder | 2 through 511 | 71.00564656314114 | 70.91253715874966 |
| Dense tail | 512 through L | 1540.0158429921617 | 1855.0867908870882 |
| Full union | 1 through L | 48.47068387175386 | 46.476222182258425 |

Each number is the margin of an upper bound, not an estimate of the true
failure probability. Q1 dominates both full unions. The much larger dense
margins do not add to the full margin.

Full ledgers:

- `generated/t128_s19_m22_ladder_full_v1.json`, SHA-256
  `f7586170c31160c0c1ca1d56999ebd5e4a00ae6cbbbbc1bfd208097592417620`.
- `generated/t128_s19_m24_ladder_full_v1.json`, SHA-256
  `a1d112ee11148d9077362eac871faeac51e42a62a3dcb5d9b3755887e752e8ea`.

For m in {22,24}, the sparse report is
`generated/t128_s19_m{m}_ladder_sparse_v1/report_0000.json`.
Each report covers Q1..511 with matching 512-bit numerical replays.
Q1 uses the size-specific coefficients and exact BCH weighted inequality.
The sparse range workers use the existing positive polynomial recurrence,
with the new row count and cutoff supplied explicitly.

The dense source is
`generated/t128_s19_m{m}_ladder_fixed_reference_v1/cover_0000.json`.
Its `replay_cover_0000.json` receipt checks every retained box at 512 bits.
The K22 cover has 241 boxes; the K24 cover has 245. Both source covers
already pass their dense allocations using the loose density factor.
The final receipts also apply the sharper factor proved in
`POISSON_DENSITY_REFINEMENT.md` by exact rational multiplication.
They are `generated/t128_s19_m{m}_ladder_dense_refined_v1.json` and their
adjacent `_replay.json` files.

The dense proof sums the all-one label explicitly. Its fixed-reference
variant pays the stated composition count before varying the input tilt.
`DENSE_TWO_TILT.md` gives the change of measure and scalar interval bound.
No sampled density substitutes for complete coverage.

## Size extension, resource bounds, and reproduction

`ladder_sparse_backend.py` preserves the frozen numerical recurrence and
extends instance validation through m24 using `ladder_instance.py`.
`verify_ladder_sparse.py` replays the new workers, and
`search_ladder_sparse.py` schedules them serially with bounded run and job
budgets. Tests check agreement with the old backend at a supported size,
producer/replay agreement at both new sizes, and rejection of a mismatched
instance. No old-size probability is substituted for a new-size bound.

The sparse polynomial degree never exceeds 511. The dense calculation uses
fixed four-state matrices, bounded caches, and a finite rectangle cover.
It does not allocate or enumerate the quadratic table of (Q, all-one count)
cases at L=131072. The search schedulers run producer/replay jobs serially;
no encoder benchmarks were run for this goal.

From `workstreams/bch_rm2sub_bridge`, use m=22 and then m=24 with fresh
output directories and names:

```text
python -B search_ladder_dense_fixed.py --directory generated/NEW_dense --m 24 --minimum 512 --seconds 180 --nodes 2000
python -B search_ladder_dense_fixed.py --directory generated/NEW_dense --m 24 --minimum 512 --verify
python -B refine_ladder_density_cover.py --source generated/NEW_dense/cover_0000.json --output generated/NEW_refined.json
python -B refine_ladder_density_cover.py --output generated/NEW_refined.json --verify
python -B search_ladder_sparse.py --directory generated/NEW_sparse --m 24 --last 511 --seconds 240 --max-jobs 40 --job-seconds 70
python -B combine_extended_ladder.py --sparse-report generated/NEW_sparse/report_0000.json --dense generated/NEW_refined.json --output generated/NEW_full.json
```

If a bounded run stops early, resume it and use the newest numbered report
or cover. A full ledger rejects missing occupancies, mismatched instances,
missing replays, or an exact union above 2^-40.

`audit_ladder_completion.py` audits all three rungs together. It freshly
evaluates every dense leaf at 768 bits, authenticates the sparse 512-bit
receipts, reconstructs Q1's BCH aggregation and all component unions, and
checks that the selected maps and setup agree across sizes. It also checks
preservation of the earlier K16, K18, and K20 ledgers. Its result is recorded
separately as `generated/t128_s19_ladder_completion_audit_v1.json`.

The combined regression suite passed 71 tests, including the extended sparse
backend, scalar bounds, partition checks, and rejection of incomplete or
cross-size full unions.

The ladder establishes the three specified message lengths, not every
intermediate message length. Next, use these proved choices in the
performance comparison; keep heuristic spectrum/curve estimation separate
from certified margin lower bounds.
