# Bounded certificate search for BCH-256

## Current dense-search entry point

Use the explicit all-one-row split for new dense BCH-256 searches. Start with
`search_split_certificates.py`, then continue gaps with `fine_split_search.py`.
The latter uses exact rational tilt indices. The pooled-band search described
below is retained as a historical interface and for its sparse certificates;
its weak endpoint bounds are not the current closure status.

For K=2^18, use `search_scaled_split.py` with `search_sparse_only.py` and
`combine_scaled_split.py`. This version skips capped cases, replays assigned
witnesses, and uses compact receipts. See the
[K18 closure and reproduction](T128_S19_M18_CLOSURE.md).

The current workflow and numerical interpretation are documented in
[T128_S19_M16_CLOSURE.md](T128_S19_M16_CLOSURE.md). Zero rows stay outside the
nonzero occupancy. The all-one-row count is explicit, and only the 12
nonconstant bands enter the ordinary-row envelope. The regression suite rejects
zero or all-one weights in that envelope.

## Earlier controller and pilot

For routine parameter selection, start with
[PARAMETER_GRID_SEARCH.md](PARAMETER_GRID_SEARCH.md). That workflow screens and
ranks a grid, then checks only a few candidates. This document describes the
optional full-coverage search, not a prerequisite for every grid point.

Follow-up: [ENDPOINT_BAND_DIAGNOSIS.md](ENDPOINT_BAND_DIAGNOSIS.md) isolates the
remaining endpoint problem to controlling mixed weight types for the tested
witness. Single-band optimization and removal of never-activated paths do not
close it; changing recurrence scaling did not improve its numerical value.

The search automates witness selection and coverage bookkeeping for the fixed
BCH-256 outer and snapshotted RM2Sub inner maps. Its default target is a
40-bit upper bound on the total expected bad-message count. It does not
optimize for 80 bits at each occupancy, and it does not assume a BCH spectrum.
The setup randomness and distance event are those in
[CURRENT_UNDERSTANDING.md](CURRENT_UNDERSTANDING.md).

## Interface and acceptance

`search_certificates.py` runs one instance at a time. The instance records the
selected map hash, outer-input manifest hash, message exponent, row count,
cutoff, and setup convention. The initial interface supports message exponents
13 through 20 when the row count contains complete inner epochs.
Supporting an exponent in the interface does not certify it.
The checker also authenticates the actual frozen BCH files against that
manifest. A manifest hash alone would not detect altered numerical inputs.

Each numerical trial has three separate stages:

1. Produce an outward bound with a saved witness at 256-bit Arb precision.
2. Replay that witness at 512 bits without optimizing it. Q1 replay also uses
   linear epoch iteration instead of binary powering.
3. Let the coverage checker include only the explicitly recorded occupancies
   from the replayed trial. Add their upper bounds as exact rational numbers.

The two enabled numerical backends are the four-state Q1 bound and the
activation-aware three-state fixed-weight range bound. Higher precision is
a check of numerical evaluation, not an independent proof of their shared
mathematical reductions. Binary64 screens and fixed-type box diagnostics
are not certificate backends in this version.

`CERTIFIED` means every integer occupancy from 1 through the number of outer
rows is covered and the exact sum is at most the requested target.
`UNRESOLVED` includes coverage gaps, weak bounds, numerical errors, and exhausted
search budgets. It does not assert that the code or a first-moment argument fails.
The reported margin for a partial sum is explicitly marked as not a full margin.

The scheduler initially reserves at most half the target for Q1. It divides the
remaining budget equally among uncovered occupancies and rounds that allowance
down to a dyadic number. This allocation is a conservative search policy,
not a necessary condition for a certificate. The verifier checks the exact
final union, not this particular allocation policy. The reused evaluator may
retain a bound of 2^-80 when it computed a smaller value; this storage convention
is not the acceptance threshold.

## Scheduling and resumption

The scheduler first tries Q1 and useful saved anchor witnesses. Imported
witnesses are recomputed on the requested instance; their old numerical bounds
are not transferred. It then alternates sparse-gap work with a bounded local
tilt search at the dense endpoint. Sparse proposals reuse a predicted tilt,
then try neighboring tilts and smaller intervals when needed.

The search uses serial child processes. Each trial's budget includes producer
and replay, and a run has both a wall-clock budget and a maximum trial count.
A timed-out child is terminated and reaped before another is started. Startup,
input authentication, and final reporting add overhead outside numerical work.
There is no background service and no benchmark is started by this driver.
Do not run multiple copies concurrently.

Each search directory has an immutable instance/target manifest, write-once
trial inputs, logs, optional numerical certificates and replay receipts, and
numbered report snapshots. Repeating a command resumes that directory: it
reuses accepted bounds and skips already attempted task identities, including
weak or timed-out trials. A producer without a completed replay does not count.
This version does not retry an interrupted replay automatically; use a separate
search directory with a larger budget if that is needed. A changed instance or
target requires a new directory. The scheduler itself may evolve without
changing the instance or invalidating its mathematical witnesses.

The selected t64_s20 instance at K=2^20 has an optional authenticated legacy
baseline. This route reconstructs the existing complete ledger and improved
Q1 bound. It does not rerun every historical producer, and its success is not
evidence that the new scheduler can rediscover the full certificate from scratch.
The historical cutoff is one larger than the new floor cutoff, so its bound
also covers the new event.

## Running and checking

From this document's directory, with the existing ignored audit inputs restored:

```powershell
C:\Python314\python.exe -B search_certificates.py --configuration t64_s16 --directory generated/search_t64_s16_example --m 20 --target-bits 40 --seconds 180 --job-seconds 120 --max-jobs 4
C:\Python314\python.exe -B verify_certificate_search.py --report generated/search_t64_s16_example/report_0000.json --fresh
```

Run the same search command again to continue, then verify the newly numbered
report. `--fresh` reruns the accepted new numerical certificates; it still does
not rerun legacy producers. Omit `--fresh` to authenticate replay receipts and
reconstruct coverage without another numerical evaluation. For the established
baseline, use `--configuration t64_s20 --reuse-legacy` in its own directory.

Only code, tests, and this interpretation belong in Git. Search directories and
large prerequisite numerical artifacts remain ignored. A clean checkout needs
the separately preserved audit inputs; the driver does not fabricate missing
receipts or silently fetch inputs from another worktree.

## Scope of the first implementation

This is an initial bounded search controller, not a complete parameter optimizer.
It has no outward fixed-type fallback, no shared polynomial cache across child
processes, and no measured encoder cost objective. In particular, automation
does not repair the weak dense-endpoint inequality identified in
[INNER_PARAMETER_OPTIMIZATION.md](INNER_PARAMETER_OPTIMIZATION.md).

The next useful extension is a diagnostic that identifies the dominant weight
band and state contribution at the endpoint, followed by a replayable stronger
witness backend if that diagnosis supports one. Expanding a large grid before
addressing that bottleneck would mostly automate inconclusive computations.

## Pilot evidence

The authoritative first pilot outputs use the local directories
`generated/search_<map>_v2/`. The earlier `_v1` directories are development
outputs superseded by the additional frozen-input authentication; do not use
them as current replay receipts.

All pilots use K=2^20 and a 40-bit total target, with at most four trials,
180 seconds of numerical-search budget per run, and 120 seconds per trial.
The baseline uses its explicit legacy-reuse option. The candidates use fresh
numerical producers and higher-precision replay, including for imported anchors.
No candidate receives the baseline's coverage.

| Selected inner | Replayed/accepted occupancies | Full 40-bit result |
|---|---|---|
| t64_s20 | 1 through 8192, authenticated legacy reuse | Certified: 50.487298 bits |
| t64_s16 | 1 through 3, and 512, after one resumed trial | Unresolved |
| t128_s19 | 1, 2, and 2620 | Unresolved |

The initial candidate runs took approximately 112 and 138 seconds, respectively.
A one-trial resume for t64_s16 automatically selected [3,10] at tilt -6.5,
added Q3, and did not repeat an earlier trial. In the integer CLI convention,
this tilt is written -65. Its final report is `report_0001.json`; the other
two pilots use `report_0000.json`.

Both candidate reports passed an additional fresh 512-bit numerical replay of
every accepted certificate followed by exact union reconstruction. Twenty-three
tests passed, including global-budget accounting, false-closure rejection,
cross-instance rejection, dependency checks, process timeout, and resume behavior.
Repository hygiene checks passed. These checks do not establish full coverage
for either cheaper map.

The automatically selected new endpoint tilts were 0.2 for t64_s16 and -0.2
for t128_s19. Both gave weak producer bounds and were excluded without spending
another numerical trial on replay. They did not improve the earlier endpoint
witnesses. The next investment should be diagnosis or a stronger backend,
not an unbounded expansion of this local tilt search.

The sparse trial [2,9] illustrates why subset accounting matters. At its single
chosen tilt it certifies Q2, but its bounds at Q3 through Q9 are weak. Only Q2
enters the report. A different tilt can be tried automatically at the next gap;
the existence of a range result does not make its entire range pass.

The endpoint search is likewise an upper-bound search. A negative diagnostic
margin means that this witness gives a weak upper bound, not that a lower bound
establishes insecurity. Weak endpoint trials are logged but do not consume
the certificate's union budget or enter coverage.
