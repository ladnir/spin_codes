# BCH-256 / selected RM2Sub (128,19): full closure at K=2^16

Later result: the same selected map is also certified at K=2^18, with a
52.3463883689-bit margin. See the [separate K18 proof](T128_S19_M18_CLOSURE.md).
The statement and receipts below concern K=2^16 only.

The selected t128_s19 construction now has a full certified distance/setup
margin of **53.9443672720 bits** at K=2^16. All 512 occupancies are covered.
The target was 40 bits. This is an outward first-moment result using the
certified BCH spectrum constraints, not the heuristic spectrum model.

The encoder is unchanged. Its previously measured online time at this size
is 0.560 ms on Peach; see [the performance report](../bare_bch_rm2sub/PERFORMANCE.md).
No new encoding benchmark was run for this proof work.

## Certified statement

Fix the BCH [256,128,d>=38] constituent and the selected t128_s19 maps
identified by the instance manifest. Set K=65536, L=K/128=512, N=2K=131072,
and H=floor(N/10)=13107.

The setup independently permutes each outer row and each transposed region.
Each inner update uses a fresh independent uniform nonzero GF(2^19) multiplier.
The state starts at zero, persists between regions, and updates after output.
There is no final flush. One sampled setup is shared by all messages.
The selected maps and encoder interface are specified in
[CURRENT_UNDERSTANDING.md](CURRENT_UNDERSTANDING.md).

For a sampled setup S, let Z_Q(S) count nonzero messages with Q nonzero
outer rows whose encoded output has weight at most H. The authenticated
ledger gives the exact rational bound

\[
\sum_{Q=1}^{512}\mathbb E_S[Z_Q(S)]\le U<2^{-53},
\qquad -\log_2 U\approx53.9443672720.
\]

Numerically, U is about 5.7693562117e-17. The exact numerator and denominator,
not that decimal approximation, determine acceptance.

If any nonzero message has output weight at most H, then the total count
is at least one. Markov's inequality therefore bounds that setup event by U.
Except with probability at most U, every nonzero message has output weight
at least 13108. This is the distance/setup margin, not the total security
level of SPIN. The result applies to this selected map and K=2^16.
K=2^18 now has its own certificate; K=2^20 remains a separate proof task.

## What closed the proof

The [all-one-row split](CONSTANT_ROW_CLOSURE.md) is now the default dense
proof-search path. Zero rows remain outside Q. At each Q, the bound sums
every all-one-row count h=0,...,Q. Only the other Q-h nonzero rows use the
12 nonconstant weight bands.

`constant_split_grid.py` computes the same matrices V_d(h) as the existing
single-occupancy evaluator. It retains every offset h during the recurrence,
so one witness serves every Q in the requested interval. Batching shares
calculations; it does not infer untested bounds from nearby occupancies.

The first nine witnesses covered most of Q=150,...,512. A bounded search
at integer tilt indices left Q=210,...,219, Q=231,...,251, and Q=284 unresolved.
These were narrow witness windows, not another failure at the all-one endpoint.
For example, the tested indices -3 and -2 do not close Q=241, whereas -5/2
does. Here an index a means lambda=exp(a/10).

`fine_split_search.py` selects the midpoint of a current gap and its weakest
all-one-row case. It proposes quarter-index tilts near interpolated anchor
tilts, with bounded neighboring trials. Each proposal receives new rational
Bernoulli probabilities and a complete outward evaluation. Interpolation
selects a witness; it never supplies a bound or coverage.

Starting from the replayed nine-witness bank, ten additional witnesses closed
the full dense interval. The accepted bank has 19 witnesses and includes
120516 distinct (Q,h) cases. Every case was replayed at 512-bit Arb precision.

## Numerical evaluation and exact aggregation

The region transfers and row-cost roots use Arb at 256 bits in the producer
and 512 bits in replay. The scaled positive recurrence uses directed binary64
upper bounds. The new terminal evaluator squares batches of 4-by-4 matrices
eight times, rounding each nonnegative product and addition upward separately.
It rescales by powers of two after each square and retains the exponent exactly.
No BLAS reduction or floating logarithm contributes to a retained bound.

Integer arithmetic converts each resulting positive upper bound into a dyadic
upper bound. The evaluator retains at most 80 bits per case. It then sums
all h cases exactly. Thus the stored margins are conservative guarantees,
not estimates of the true spectrum-weighted margin.

The final ledger authenticates the existing sparse report and its replay
receipts, then merges the replayed dense bounds by occupancy. It uses the
smaller valid upper bound where coverage overlaps. All integers Q=1,...,512
must be present before the ledger can report `CERTIFIED`.

| Quantity | Retained margin |
|---|---:|
| Q=1 contribution | 54.0067373509 bits |
| Q=150,...,512 contribution from the new bank | 58.5073428910 bits |
| Full exact union, Q=1,...,512 | 53.9443672720 bits |

The earlier combined partial ledger reported 48.7642 bits. Some of its loose
sparse bounds were in the new dense interval. Replacing those bounds explains
why the full union is stronger than that partial ledger. No code change or
assumed spectrum improvement accounts for this gain.

Replay uses higher-precision coefficients with the same mathematical bound
and directed terminal algorithm. It is a numerical cross-check, not an
independent formalization of the transfer lemmas. The new terminal arithmetic
also passed comparisons against scalar Arb matrix powering.

## Records and reproduction

The accepted local records are:

- Sparse report: `generated/closure_t128_s19_m16_v1/report_0002.json`.
- Initial bank and replay: `generated/t128_s19_m16_split_grid_v1.json` and
  `generated/t128_s19_m16_split_grid_v1_replay.json`.
- Final bank and replay: `generated/t128_s19_m16_fine_split_grid_v1.json` and
  `generated/t128_s19_m16_fine_split_grid_v1_replay.json`.
- Full ledger: `generated/t128_s19_m16_full_split_coverage_v1.json`.

The full ledger SHA-256 is
`3684d355322ee59d9cffdeb3ac5edd29ff3bdc22ca37b73fb723a4cbe95b0b18`.

The integer-only continuation `t128_s19_m16_split_grid_v2.json` is a partial
producer, has no replay, and does not enter the full ledger.
Its controller, `refine_split_search.py`, is retained for diagnosis;
use the rational-tilt controller for new continuations.

Run the following serially from this directory, with the retained audit inputs
available. These example output paths must be fresh:

```powershell
C:/Python314/python.exe -B search_split_certificates.py --output generated/example_split_seed.json
C:/Python314/python.exe -B search_split_certificates.py --output generated/example_split_seed.json --verify
C:/Python314/python.exe -B fine_split_search.py --seed generated/example_split_seed.json --output generated/example_split_fine.json --max-jobs 16 --seconds 180
C:/Python314/python.exe -B fine_split_search.py --output generated/example_split_fine.json --verify
C:/Python314/python.exe -B fine_split_search.py --output generated/example_split_full.json --sparse-report generated/closure_t128_s19_m16_v1/report_0002.json --grid generated/example_split_fine.json
```

A continuation stops on complete per-occupancy coverage or its proposal budget.
The time budget is checked between witnesses; the final witness can exceed it.
Authentication, replay, and final aggregation are outside that search budget.
To resume an incomplete run, replay its bank and use it as `--seed` with a
fresh output path. The automatic proposal policy currently supports only
t128_s19 at m=16. Supporting another size numerically does not certify it.

The targeted bridge suite passed 32 tests. It covers the actual pooled-failure
versus split-success endpoint, the fractional-tilt window at Q=241, every offset
of small batched recurrences, exact dyadic rounding, and incomplete-case rejection.
The existing seven activation-density tests also pass.

Generated numerical banks remain ignored; each new accepted bank is about
1.6 MB. Commit source, tests, and these notes, not the experimental directories.
Previously frozen producers and their receipts remain unchanged.

## Next step

The K=2^18 extension is complete. Move next to K=2^20, reusing the structural
split and rational-tilt search but recomputing coefficients, witnesses, and the
full union. Q1 already dominates the completed K=2^16 bound; spending more
effort on its dense occupancies is not the next priority.
