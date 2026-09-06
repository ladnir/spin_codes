# First activation-aware parameter tranche

This study asks how constituent size, message length, epoch length, and state
size affect the finite Structured SPIN bound. Exact-spectrum BCH and RM
constituents supply the primary data. The bounded-spectrum BCH-256 certificate
remains a secondary comparison and does not enter this tranche or a fit.

The first completed tranche contains 715 preferred occupation-one (Q1)
observations. All use the activation-aware transfer and nearest binary64
arithmetic. These are diagnostic bounds, not complete distance certificates.
The database also retains 143 superseded coarse-grid observations and all 414
historical observations, giving 1,272 indexed results.

## Fixed construction and comparison grid

Fix one rate-one-half outer constituent of dimension D and length B. Encode
each of the L=k/D message rows with this same constituent. Independently
permute coordinates within each encoded row, transpose into B regions, and
independently permute positions within each region. The output length is N=2k.

Within each t-bit epoch, the RM2Sub recurrence is

    q_0 = 0;  Y_i = X_i + A q_i;  q_(i+1) = alpha_i q_i + A^T X_i.

The state has s bits. Each alpha_i is independently uniform on the nonzero
elements of GF(2^s). State continues across regions and there is no final
flush. One setup fixes all permutations and multipliers for every message.
The Q1 quantity counts bad messages with exactly one nonzero outer row.
The bad event is output weight at most floor(N/10).

The pilot uses:

- message exponents log2(k) in {16,18,20,22,24};
- t in {64,128,256};
- s in {14,16,18,20} at each t, plus s=19 at t=128;
- exact BCH [8,4,4], [32,16,8], [64,32,12], and [128,64,22];
- exact RM(1,3), RM(2,5), RM(3,7), and RM(4,9); and
- random full-rank outer references [64,32], [128,64], and [512,256].

This is 13 inner configurations times five message lengths times eleven
outer models. Every tuple satisfies D | k and t | L. No wrapper is used.
The current inventory has four exact BCH constituents; a fifth has not been
identified in this checkout.

For each t, one unselected ordered RM(2,log2(t)) basis supplies all tested
states by prefix restriction. Its seed is 3390173185+t. This controls map
variation when s changes. A comparison across t still includes variation
between these three chains. It does not establish the typical behavior of
sampled maps or an optimal map at any t.

Every map was checked for full rank, A^T A=0, and distinct nonzero columns
of A^T. The A spectrum was exhaustively enumerated. The complete kernel
spectrum was obtained by an integer MacWilliams transform, with an independent
pair-collision check of its weight-four coefficient. These algebraic audits
do not constitute a distance proof.

Random references use the exact ensemble expectation for one uniformly
sampled full-rank outer, reused in every row. Linearity of the Q1 sum permits
this average. It must not be substituted into products of spectra for higher
occupations.

## Why the inner transfer changed

The old family driver uses a two-state live-law invariant. After activation
from zero, the next state is A^T e_J, supported on only t possibilities.
The BCH bridge showed that the asserted near-uniform live law need not hold
at that point. An exact outer spectrum does not repair this invariant.

`activation_q1.py` uses the BCH bridge's three weighted-measure classes:
zero, arbitrary nonzero, and nonzero with bounded density. Activation enters
the arbitrary class. A later fresh multiplier supplies the live-law bound.
The transfer uses the exact nonzero A spectrum for the density-bounded class
and minimum A weight for the arbitrary class.

The kernel spectrum is recorded for subsequent higher-occupation work. Q1
only needs nonzero columns of A^T; it does not test the kernel behavior that
can dominate middle or dense occupations.

The implementation batches tilt witnesses, keeps all coefficients in the
log domain, and uses degree-one matrix-polynomial squaring for long regions.
Its matrix contractions retain three explicit states. It reuses the complete
coefficient array when different outer spectra share the same geometry.

Exact rational tests independently enumerate outer support patterns and
compare sequential positive region products with the batched log calculation.
A separate GF(4) encoder calculation checks actual first moments against the
envelope on a small instance. These checks validate the implementation on
those instances; binary64 results remain diagnostics.

## Observed epoch/state trade-off

At k=2^20, the following Q1 margins use exact outer spectra:

| Outer | s | t=64 | t=128 | t=256 |
|---|---:|---:|---:|---:|
| BCH [128,64,22] | 14 | 31.243 | 31.529 | 31.658 |
| BCH [128,64,22] | 18 | 32.669 | 32.676 | 32.675 |
| BCH [128,64,22] | 20 | 32.758 | 32.752 | 32.744 |
| RM(4,9) [512,256,32] | 14 | 38.164 | 38.444 | 38.572 |
| RM(4,9) [512,256,32] | 18 | 39.137 | 39.146 | 39.147 |
| RM(4,9) [512,256,32] | 20 | 39.188 | 39.184 | 39.177 |

Within these chains, increasing s from 18 to 20 buys less than 0.1 Q1 bit
at k=2^20. This does not say that those state bits are unnecessary for
higher occupations. All three epoch sizes should remain implementation
candidates until those occupations have been compared.

The preferred epoch also depends on k. For BCH [128,64,22] at s=14,
t=64 gives 35.558 bits and t=256 gives 35.254 bits at k=2^16.
At k=2^24, the respective margins are 18.592 and 25.518 bits.
Thus a fixed ranking of epoch sizes would miss a substantial scaling effect.
At k=2^24 and s=20, the margins return to a narrow range, 28.635--28.743 bits.

The plotted comparisons are `activation_ts_tradeoff.png` and
`activation_k_scaling.png`. They show Q1 screens, not measured failure rates.

## Witness-grid correction

The first runs used log-surprisal witnesses {-12,-11.9,...,0}. At k=2^24,
some dominant shells selected the lower grid endpoint. Every configuration
at that message size was recomputed on {-16,-15.9,...,0}, which contains
the entire original grid. The largest aggregate improvement was 6.240339 bits.
No preferred observation now has a dominant witness at a grid endpoint.

An interior dominant witness does not prove continuous optimization. It
removes an observed truncation of this exploratory grid. The coarse receipts
remain unchanged and indexed with comparison_eligible=0.

## Evidence, cost, and database outputs

The new manifests bind the CSV files, maps, spectra, and producer sources by
SHA-256. Imports reject a changed dependency or CSV. Git attributes preserve
the line endings used by these byte hashes. The original 414 result IDs are
preserved.

`transfer_review_status` distinguishes activation-aware screens from historical
results pending review. Historical certificate labels remain in `result_class`
for provenance. `certified_results` excludes entries under review, so the old
RM(4,9) receipt cannot silently supply a current certificate anchor.
`comparison_eligible=1` selects the 715 preferred Q1 observations.

`activation_q1_frontier.csv` records the smallest tested s clearing each of
0, 20, 40, and 60 Q1 bits, separately for every outer, t, and k. A missing
entry means that no tested state passed that Q1 screen. It does not prove
failure of the code, and a passing entry does not prove full distance.

`activation_pilot_audit.json` records direct XOR reduction counts for A and
A^T, output addition, and epoch frequency. These are preliminary operation
counts. They exclude field multiplication, state-update addition, routing,
outer encoding, circuit sharing, and SIMD. They do not rank runtime.

## Reproduction and next experiments

From this directory, rebuild the index and reports with:

```text
python build_landscape_db.py
python query_landscape.py export landscape_export.csv
python summarize_activation_pilot.py
python -m unittest -v test_landscape_db.py test_extrapolate_parameters.py test_activation_q1.py
```

The sequential screen commands used three output directories:

```text
python run_activation_pilot.py --output-dir activation_pilot_v1
python run_activation_pilot.py --output-dir activation_pilot_scaling_v1 --message-exponents 18 22 24
python refine_activation_grid.py --source activation_pilot_scaling_v1/q1.csv --output-dir activation_pilot_e24_refined_v1
```

Use new directories when rerunning numerical jobs; existing screen CSVs are
write-once. Run only one numerical job or benchmark at a time.

The next smallest useful tranche is:

1. Extend the same chains to s=12,13,15 near the small-k Q1 frontier, keeping
   the existing higher-state members fixed. Compare matched persistence as
   well as fixed s.
2. Implement activation-aware Q2 and selected sparse, middle, and dense
   screens for the exact-spectrum constituents. Include both an economical
   state choice and a larger-state control at each t.
3. Fit BCH and RM scaling only after the inner screens are comparable. Hold
   out a constituent size and a message exponent. Keep Q1 and full-occupation
   models separate. Retain prediction errors, signed margins, and map identity.
4. Compare the completed BCH-256 certificate as a bounded-spectrum secondary
   anchor. Import compact evidence after review; do not treat a shell cap as
   an exact spectrum measurement.
5. Benchmark the competitive choices only after their certificate landscape
   is useful. Keep t=64,128,256 until map and multiplication costs are measured.

No new full-distance claim, extrapolation fit, or implementation benchmark was
produced in this tranche. The older Q1 fit is retained explicitly as historical.
