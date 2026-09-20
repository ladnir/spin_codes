# BCH-256 / selected RM2Sub (128,19): full closure at K=2^18

The selected t128_s19 construction has a full certified distance/setup margin
of **52.3463883689 bits** at K=2^18. All 2048 occupancies are covered, and
the exact union clears the 40-bit target. This result recomputes the bounds
at K=2^18; it does not transfer numerical bounds from K=2^16.

## Statement and scope

Fix the BCH [256,128,d>=38] constituent and selected t128_s19 maps identified
by the instance manifest. Set K=262144, L=K/128=2048, N=2K=524288, and
H=floor(N/10)=52428.

Use the setup and encoder from
[CURRENT_UNDERSTANDING.md](CURRENT_UNDERSTANDING.md): independent row and
region permutations, fresh independent uniform nonzero GF(2^19) multipliers,
zero initial state, output before update, persistent state across regions,
and no final flush. One sampled setup is shared by all messages.

For a sampled setup S, let Z_Q(S) count nonzero messages with Q nonzero outer
rows and encoded output weight at most H. The authenticated exact ledger gives

\[
\sum_{Q=1}^{2048}\mathbb E_S[Z_Q(S)]\le U<2^{-52},
\qquad -\log_2 U\approx52.3463883689.
\]

Here U is approximately 1.7464933497e-16. Acceptance uses its exact rational
value, not this decimal approximation. If any such message exists, the total
count is at least one. Thus Markov's inequality bounds that setup event by U.
Except with probability at most U, every nonzero message has output weight
at least 52429.

This is a distance/setup guarantee, not the total security level of SPIN.
It uses certified BCH spectrum constraints rather than an assumed spectrum.
The selected map's K=2^20 instance remains uncertified.

## Coverage and contribution

| Contribution | Covered occupancies | Retained margin |
|---|---|---:|
| Q1 | 1 | 52.3600177467 bits |
| Sparse certificates | 1 through 479 | Included in exact full union |
| All-one-split certificate | 480 through 2048 | 59.0794486935 bits |
| Full exact union | 1 through 2048 | 52.3463883689 bits |

The dense certificate includes all h=0,...,Q all-one-row counts at each Q.
It contains 1984785 cases. Every retained case bound is 2^-80, so its exact
contribution is 1984785/2^80. The 80-bit cap limits retained precision; it
does not estimate the actual per-case margin. Q1 dominates the full union.

For comparison, the separate [K=2^16 certificate](T128_S19_M16_CLOSURE.md)
has a full margin of 53.9443672720 bits. These are conservative proved points,
not a fitted estimate of the unknown BCH spectrum.

## Search and verification

The proof uses the same explicit all-one-row split as K=2^16. Zero rows are
outside Q. The all-one-row count is explicit, and only the 12 nonconstant
bands enter the ordinary-row envelope. State persists through all regions.

`search_scaled_split.py` scales witness locations from the K=2^16 bank,
then recomputes region transfers, Bernoulli probabilities, and every retained
bound for the requested size. Scaling proposes witnesses; it supplies no bound.
Fifteen witnesses close the dense interval. Quarter-index tilts fill the gaps
left by the initial coarse witnesses, with lambda=exp(index/10).

`scalable_split_grid.py` reuses the frozen directed recurrence and fixed-width
four-state terminal arithmetic. It skips a case only after that case has
reached the retained 80-bit cap. Consecutive binomial coefficients are computed
by exact integer recurrence. These changes reduce repeated work without changing
the inequality.

Each retained case records the witness that supplies its bound. The 512-bit
replay recomputes that assigned witness and checks the stored upper bound.
It verifies every case exactly once and rejects missing or invalid assignments.
This is a numerical replay of the same mathematical bound, not an independent
formalization of the transfer lemmas.

The dense producer took about 211 seconds of numerical search locally.
Run-length encoding reduces its receipt to 512461 bytes. Its replay checked
all 1984785 cases. The original K=2^16 sources and receipts remain unchanged.

`search_sparse_only.py` handles Q1 and the sparse ranges using the frozen
producer and replay backend. It excludes pooled dense-endpoint trials.
The first bounded run reached Q=456; its final trial exhausted the time budget.
A fresh retry of that trial certified Q=457,...,479. Every accepted sparse
certificate has a 512-bit replay. Q1 replay also uses linear epoch iteration.

The combined ledger authenticates both inputs, checks instance identities,
and adds all accepted upper bounds as exact rationals. No heuristic spectrum,
interpolated coverage, or bound from another message length enters the union.

## Records and reproduction

Accepted local records:

- `generated/t128_s19_m18_sparse_v1/report_0001.json` and its referenced certificates/replays.
- `generated/t128_s19_m18_scaled_split_v1.json`.
- `generated/t128_s19_m18_scaled_split_v1_replay.json`.
- `generated/t128_s19_m18_full_split_coverage_v1.json`.

The full ledger SHA-256 is
`acac24114eb25d658ffcd86ee3a7fe2c012f0c156e0fa73fb50cbfa4f5819550`.
The dense producer SHA-256 is
`39824756f5963381fadd10a7b2409426ccf5d968333a8bb5900341f79d56acfe`.

From this directory, run serially with the retained audit inputs available.
Example output paths must be fresh:

```powershell
C:/Python314/python.exe -B search_scaled_split.py --output generated/example_m18_dense.json --m 18 --first 480 --last 2048 --max-jobs 28 --seconds 600
C:/Python314/python.exe -B search_scaled_split.py --output generated/example_m18_dense.json --verify
C:/Python314/python.exe -B search_sparse_only.py --directory generated/example_m18_sparse --m 18 --last 479 --seconds 240 --max-jobs 40 --job-seconds 40
```

Repeat the sparse command if its report has gaps. The controller retries an
interrupted logical task once with fresh output paths. It does not treat a
timeout as a weak bound. Use the newest numbered sparse report when combining:

```powershell
C:/Python314/python.exe -B combine_scaled_split.py --output generated/example_m18_full.json --sparse-report generated/example_m18_sparse/report_0001.json --grid generated/example_m18_dense.json
```

The report number depends on how many runs were needed; do not assume it is
always 0001. Dense continuation requires a replayed `--resume` bank with the
same instance and interval. Time limits are scheduling limits, not acceptance
conditions. Full closure requires complete coverage and the exact union test.

The targeted bridge suite passes 43 tests, including the existing all-one
endpoint and fractional-tilt regressions. New tests compare selected evaluation
with the full directed grid, verify exact binomial counts, reject malformed
case encodings, and exercise sparse-only scheduling and timeout retry.

No encoder code or performance measurement changed. The existing
[performance report](../bare_bch_rm2sub/PERFORMANCE.md) records 2.322 ms at
K=2^18 for this implementation. Numerical artifacts remain ignored; preserve
source, tests, and the proof notes in Git.

## Next step

Close K=2^20 for the same selected map. Use occupancy intervals to bound memory:
the scaled controller intentionally rejects tables larger than ten million cells.
Recompute all bounds and the full union at that size. Further tightening the
completed K=2^18 dense bound has low payoff because Q1 already dominates.
