# Finite IMT migration

Goal opened 2026-09-15: replace the finite RM2Sub certificates and SPIN
performance numbers in the paper with results for IMT. Keep the historical
sources and evidence intact. Do not commit generated data.

Scope refined 2026-09-16: use the complete 130-cell IMT Q1 grid for the
engineering slices, alongside five selected, exactly map-matched Q1/full
certificate comparisons. These four plots are now integrated. Full proofs
for every diagnostic grid cell are no longer a completion requirement.
The selected BCH-256 points do not certify endpoints of the different
BCH-64/128 diagnostic maps. See [current scope](PARAMETER_SLICES.md).

## Replacement inventory

| Paper or artifact item | Required IMT replacement | Starting status |
|---|---|---|
| BCH-256 theorem and engineering curve | Rate 1/2, distance 10%, K=2^16,2^18,2^20,2^22,2^24; full union below 2^-40 at each size | All five lengths closed with the same weight-five feedback |
| Finite proof appendix | Exact independent maps, transvection law, all occupancy reductions, authenticated unions | Draft now uses shared IMT transfers, short-length refinements, and quarter-rate reductions; final review remains |
| BCH-64/128 parameter plots | Q1 sensitivity slices with explicit fixed IMT maps | All 130 cells authenticated and integrated; five selected BCH-256 Q1/full comparisons shown separately |
| Selected half-rate performance table | Complete transpose at K=2^16,2^18,2^20, matching certified maps | Matched serial series authenticated and bound to full certificates at all three sizes |
| Reference inner performance row | Measured IMT alternative, or an explicitly historical comparison outside the current-family table | Old (64,20) RM2Sub row is not an IMT result |
| External-encoder comparison | Rerun the SPIN row with IMT under the same in-place timing policy; reconcile comparison text and ratios | Other families are external baselines, not replacement targets |
| Quarter-rate operating points | BCH [128,32,32], K20, 16.5%/40 bits and 19%/30 bits, with matched transpose measurements | Both IMT certificates and an optimized implementation exist; authenticate and integrate |
| Abstract, introduction, artifact map, checks | Values and links drawn from the accepted replacement records | Update after the relevant cells pass |

The parameter plots are diagnostic evaluations, not outward certificates.
They must retain that distinction after migration. Existing external BAA,
Expand-Convolute, and RAA results do not become IMT results and remain scoped
to their own ensembles and measurement conditions.

## Execution order

Current milestone: all selected half-rate certificates and both quarter-rate
operating points are closed and integrated into the manuscript with matched
timings. [`PAPER_RESULTS.md`](PAPER_RESULTS.md) is the compact integration
ledger. The Q1-focused parameter-slice replacement is integrated; final
manuscript/release audit remains. The historical progress notes below are not the
current certificate inventory.

1. Extend the weight-five BCH-256 engine to the upper ladder without changing
   frozen producers. Recompute Q1, sparse, and dense contributions at each size.
2. Close K16 and K18 using tighter finite bounds. The existing pooled dense
   expression fails there; its failure does not establish poor code distance.
   Adapt the explicit all-one-row split before changing the selected maps.
3. Recompute the smaller-BCH parameter slices and authenticate quarter-rate
   proof and implementation records. Preserve useful slice coverage rather
   than erasing plots to make the migration appear complete.
4. Measure matching implementations serially, with correctness and map binding.
   Reuse fixed-width optimized kernels and routing improvements. No concurrent
   benchmarks. Complete forward implementation is a separate follow-up; this
   migration replaces the paper's transposed-encoder measurements.
5. Update the finite construction, proofs, tables, plots, abstract, and artifact
   checks together. Compile and visually inspect the manuscript.

The goal is complete only when every current-family finite claim and SPIN
timing in the paper has an IMT source. A missing certificate is not replaced
by interpolation, a partial-Q bound, or a timing from another configuration.
If a required target cannot close, record the obstruction and request a
scope decision rather than silently weakening the theorem.

## New sources

`ladder.py` supplies a length adapter for the unchanged weight-five maps and
transfers. Its Q1 command produces only a one-active-row bound, never a full
distance certificate. Replay authenticates sources and uses 512-bit arithmetic
with linear epoch iteration. New records use fresh output paths.

```text
python -B workstreams/inner_design/finite_migration/ladder.py --mode q1 --m 22 24 --output workstreams/inner_design/finite_migration/Q1_UPPER.json
python -B workstreams/inner_design/finite_migration/ladder.py --mode q1 --output workstreams/inner_design/finite_migration/Q1_UPPER.json --verify
python -B -m unittest discover -s workstreams/inner_design/finite_migration -p "test_*.py"
```

Next: adapt the small-length composition split, then replace the diagnostic
plots and measurements. No manuscript number changes merely because a job
was launched.

## Initial progress

The length adapter passes five tests, including agreement with the old
finite producer at K16. Geometry and mixed-witness tests add three checks.
The four existing exact-union rejection tests also pass.

At K22 and K24, the one-active-row bounds have margins 48.3904952900 and
46.4572568186 bits. Both passed 512-bit replay with linear epoch iteration.
The sparse ranges Q=2..511 have margins approximately 67 and 66 bits and
also passed 512-bit replay. These components alone are not full certificates.

K22 now has a complete certificate: all Q=1..32768 are covered, with an exact
union margin of **48.3904916829 bits** at distance 10%. Its 339-leaf dense
cover passed 512-bit replay. The numerical inputs are `Q1_UPPER.json`,
`SPARSE_M22.json`, and `DENSE_M22_mixed_v2.json`, each with a replay receipt.
`FULL_M22_VERIFIED.json` authenticates their source bindings, the unchanged
outer caps, their common instance, complete coverage, and the exact sum.
Implementation binding at this length remains separate; no K22 timing is
claimed.

K24 is also complete: all Q=1..131072 are covered, with an exact union
margin of **46.4572549297 bits** at distance 10%. Its 359-leaf dense cover
passed 512-bit replay. The inputs are `Q1_UPPER.json`, `SPARSE_M24.json`,
and `DENSE_M24_mixed_v2.json`, with their replays. The aggregate is
`FULL_M24_VERIFIED.json`. No K24 timing is claimed.

| log2(K) | Q1 margin | Sparse Q=2..511 margin | Dense Q=512..L margin | Full margin |
|---:|---:|---:|---:|---:|
| 22 | 48.3904952900 | 66.9999999993 | 84.9999972482 | 48.3904916829 |
| 24 | 46.4572568186 | 65.9999999987 | 191.5121599662 | 46.4572549297 |

These are conservative upper bounds on the shared-setup failure event,
not estimates of actual failure probabilities. Different retained witness
quality explains why the dense component margins need not follow a smooth
curve. Q1 dominates each full union. The calculations use the fixed BCH-256
outer envelopes and the exact weight-five IMT maps already certified at K20.

The initial dense search exposed a discovery limitation: it optimized only
the Fourier transfer, although the verifier accepts either complete moment
bound. At K22, Q=512, density coordinate 1/4, an occupation-transfer proposal
changes a failing bound to a margin above 11,000 bits. This is a point check,
not a certificate for an interval. `mixed_dense.py` now proposes from both
families and evaluates each proposal with the unchanged outward bound.
It never mixes entries between transfer matrices.

The K16/Q218 and K18/Q870 bottleneck points at coordinate 21/128 failed
with both original proposal families. Subsequent refinements now pass both
points; see the short-length progress below. Full coverage remains separate.

`dense_ladder.py` reshapes an authenticated subdivision tree to a new root
and recomputes every leaf bound. `sparse_ladder.py` scales only witness
proposals, then recomputes every sparse contribution. `verify_ladder.py`
requires authenticated replays, exact contiguous coverage, and an exact
union below 2^-40 before it reports a full certificate.

To replay and aggregate K22, use the local records above:

```text
python -B workstreams/inner_design/finite_migration/verify_ladder.py --m 22 --q1 workstreams/inner_design/finite_migration/Q1_UPPER.json --sparse workstreams/inner_design/finite_migration/SPARSE_M22.json --dense workstreams/inner_design/finite_migration/DENSE_M22_mixed_v2.json
python -B workstreams/inner_design/finite_migration/verify_ladder.py --m 24 --q1 workstreams/inner_design/finite_migration/Q1_UPPER.json --sparse workstreams/inner_design/finite_migration/SPARSE_M24.json --dense workstreams/inner_design/finite_migration/DENSE_M24_mixed_v2.json
```

That command authenticates the retained numerical replays; it does not
rerun them. The producer commands accept `--verify` for numerical replay
and require a fresh replay output path. No numerical artifact is committed.

## Short-length progress

[`SHORT_LENGTH_BOUNDS.md`](SHORT_LENGTH_BOUNDS.md) derives the new activation,
band-sum, routing-density, and fixed-input-tilt bounds. The all-one split alone
did not close the short instances. The useful additional gains came from
reducing the routing comparison penalty and avoiding a composition-count
penalty in a distinct fixed-input-tilt calculation. Neither change alters
the encoder or its setup distribution.

`FIXED_INPUT_POINTS_v1.json` and its 512-bit replay certify only these points:

| log2(K) | Q | Density coordinate | Point margin |
|---:|---:|---:|---:|
| 16 | 218 | 21/128 | 53.0405001895 |
| 18 | 870 | 21/128 | 654.2693726381 |

The one-active-row calculations in `Q1_LOWER.json` give 43.5928335021 bits
at K16 and 50.1891811139 bits at K18; both passed 512-bit linear replay.
In particular, the old RM2Sub K16
margin cannot be retained for IMT. The new K16 bound has less headroom above
the requested 40 bits; these numbers are conservative proof bounds, not
predictions of the actual failure probability.

Length-scaling the K20 sparse witness bank was insufficient at K18. Fresh
reference probabilities and tilts close every Q=2..127, as recorded in
`SPARSE_M18_v3.json` with a 512-bit replay. The failed scaled bank is retained
as `SPARSE_M18_v1.json`; its partial coverage must not be used as a certificate.
`short_dense.py` now searches a full dense partition beginning at Q=128.
An unresolved leaf remains unresolved even when a nearby point passes.
The bounded `DENSE_M18_short_v2.json` search retained 131 leaves with 40
unresolved regions. Continuing that tree produced `DENSE_M18_short_v3.json`:
181 leaves cover every Q=128..2048, with no unresolved regions. Its 512-bit
replay passed. `FULL_M18_VERIFIED.json` authenticates all three components
and gives a full union margin of **50.1890763816 bits**. Component margins
are 50.1891811139 (Q1), 63.9389235201 (Q=2..127), and 79.9032795381 (dense).

K16 remains open. Fresh discovery closes Q=2..63 in `SPARSE_M16_v2.json`,
with a 512-bit replay. The dense cover starts at Q=64. `budget_dense.py`
now checks a total dense budget of 2^-42 rather than requiring 80 bits
from every leaf. `verify_budget.py` independently enforces the final sum
below 2^-40. The first budgeted search, `DENSE_M16_budget_v1.json`, retained
341 leaves after 160 refinements but did not meet its budget. It is not
a certificate. Its center diagnostics exposed a Q=206, v=49/256 point
with only 29.7814336035 bits, so subdivision alone was insufficient.

`retune_fixed_input.py` selects witnesses against the actual new expression.
At that point it gives **31.4197780178 bits**, recorded in
`M16_DIRECT_POINT_v1.json` with a 512-bit replay. Separating assignments
with at most 1/32 all-one rows from those with at least 1/32 gives
**47.5100238751 bits**. `M16_CONSTANT_POINT_v1.json` and its 512-bit replay
retain this result. The rare, all-one-rich subset contributes a bound with
1958.19 bits. This split changes the proof bound, not the encoder or setup.

`constant_thresholds.py` repeats this calculation for cutoffs from 1/32
through 1/2048. `M16_CONSTANT_THRESHOLDS_v1.json` and its 512-bit replay
give the strongest tested point result at cutoff 1/1024:
**50.3999770231 bits**. The exceptional subset then has 82.85 bits.
This is a sharper point witness, not a completed dense cover.

`constant_density.py` implements the constrained variance calculation;
`constant_point.py` checks the point; `constant_dense.py` applies the split
to rectangles and uses the total-budget driver. The derivation is in
[`SHORT_LENGTH_BOUNDS.md`](SHORT_LENGTH_BOUNDS.md). A new dense-cover search
is required before promoting the point result to a full K16 certificate.

An exact joint-map audit in `joint_overlap.py` also proves rank 38 and
minimum distance at least seven for the combined A/B code. It tightens
the overlap cap but improves this point by only about 0.00004 bits.
It is not used in the all-one-count split or its cover search.

The first split cover, `DENSE_M16_constant_v1.json`, recomputed the prior
partition and used its 300-second refinement budget for 90 refinements.
It retains 431 leaves, with `budget_met: false` and dense union margin
-79.7584094638 bits. This improves the earlier incomplete cover's
-153.7114662197 bits but remains uninformative as a distance guarantee.
`zero_constant_dense.py` subsequently implemented the integer split h=0
versus h>=1. `DENSE_M16_zero_v1.json` retains 702 leaves after 271
refinements, with dense union margin -52.4435354395 bits and an unmet budget.
The five worst-leaf centers include new points with 42.58 and 44.39 bits;
these are retained in `M16_ZERO_CENTERS_v1.json`.

The useful next refinement splits out assignments containing any band
whose reference probability exceeds 3/4, not just the all-one band.
`high_band_dense.py` implements this split with a fixed auxiliary tilt.
`M16_HIGH_BAND_POINTS_v1.json` gives **67.04, 90.08, 158.82, 81.94, and
68.62 bits** at those five centers; every point passed 512-bit replay.
A full-cover search is separate. A failed trial using cutoff 1/2 is retained
but not promoted: its exceptional subset made the bound worse.

`DENSE_M16_high_band_v1.json` completed its bounded search with 887 leaves
after 185 refinements. Its dense union margin is -41.1776627023 bits and
`budget_met` remains false. Thus the stronger center results have not yet
closed the K16 cover. No background certificate or benchmark process remains.

The fast combined continuation completed with a dense union margin of
**-30.721989 bits** in `DENSE_M16_fast_combined_v1.json`.
Its total budget remains unmet; it is not a new certificate. A separate
endpoint-convexity backend, `convex_input_dense.py`, is tested but has not
completed a full cover. It must not be used to relabel this failed receipt.

The tests cover seven-state positive powering, the band-norm sum,
activation densities, routing ratios against exact toy distributions,
constrained variance vertices, joint-map hypotheses, exact union acceptance,
and uniform interval moments. They supplement the derivations; numerical
replay is not an independent mathematical proof implementation.
The current suite has 62 passing tests, including causal moments, exact
type-box coverage, no-constant map selection, and logit interpolation.

`combined_dense.py` retains both complete split bounds; it does not mix
their matrix entries. `fast_combined_dense.py` uses the same bounds with
`fast_density.py`, which skips a variance program only when a point-mass
upper bound of one already makes that weight irrelevant to the maximum.
Both constrained and unrestricted factors are checked against the original
calculations. Use this backend for the next cover continuation, retaining
the completed high-band search as its input geometry.

## Parameter-slice replacement

[`PARAMETER_SLICES.md`](PARAMETER_SLICES.md) specifies the new fixed map
chains, all 130 geometries, scope, and reproduction commands. The complete
Q1 grid is in `PARAMETER_Q1_v1.json`; `PARAMETER_Q1_VERIFIED_v1.json`
reconstructs all 39 map pairs and recomputes every cell. Seven cells also
passed a separate log-domain calculation. No old RM2Sub value was relabeled.
Higher-occupancy diagnostics remain under development, so the manuscript's
full-margin curves are not replaced by these Q1-only values.
For the BCH-128/t64/s20/K20 pilot, balanced reference probabilities and
fixed band compositions give Q2/Q3/Q4 margins 34.153247, 53.990145, and
72.967441 bits. The sparse verifier replayed all 990 composition witnesses.
The higher-occupancy tail still needs a useful bound.

The pilot's Q5..64 aggregate is now 83.534893 bits, with every witness
replayed through a separate log-domain product. Dense-search verifiers
add exact integer type-box coverage checks. The original chain's remaining
worst types are mostly all-one outer rows, and their singleton bounds also
fail. A new chain excludes the all-one expansion word, matching the restriction
already used in the selected BCH-256 implementation. Initial point checks
remove that gap while leaving Q1 essentially unchanged. The two chains have
separate receipts; neither pilot yet replaces a full parameter curve.

## Quarter-rate evidence check

The retained quarter-rate certificate and its implementation-verification
receipt passed source authentication, and the exact dyadic unions passed
again: 41.0481676058 bits at 16.5% distance and 30.0334910634 bits at 19%.
This check reused the authenticated numerical replay; it did not rerun it.

The optimized deployment's 15 compiled-source bindings also matched. Its
three K20 process medians are 15.254337, 15.312936, and 15.204465 ms, giving
a median of 15.254337 ms. These are existing serial in-place measurements,
not a new benchmark. `quarter_binding.py` now checks the explicit bridge from
the certified base implementation to this routing-optimized build. It
reconstructs the exact map header, repeats exact unions and coverage checks,
reconstructs the reviewed source transformation, authenticates compiled-source
and correctness receipts, and checks timing identities and policy.
`QUARTER_DEPLOYMENT_VERIFIED.json` retains the result. No smaller-K
quarter-rate distance claim is inferred.

## Fresh half-rate transpose measurements

`measure_half.sh` ran on the retained isolated Peach build, after checking
its complete source and binary hash lists. All 12 correctness tests passed
again. It measured three variants sequentially at K16, K18, and K20, with
three processes per cell and 101 in-place calls after three warmups. CPU 15,
GCC 15.2, znver4, packed-24 routing, and the existing tile choices were retained.
Setup and workspace preparation are outside the interval; input is neither
copied nor reset between calls.

| log2(K) | Supported RM2Sub (ms) | Routing-tuned RM2Sub (ms) | IMT (ms) | Reduction vs supported |
|---:|---:|---:|---:|---:|
| 16 | 0.562380 | 0.569613 | 0.523367 | 6.94% |
| 18 | 2.305955 | 2.323608 | 2.153740 | 6.60% |
| 20 | 11.138372 | 10.755226 | 10.109862 | 9.23% |

Each entry is the median of three process medians, not the fastest call.
The K20 result is consistent with the earlier six-process IMT result of
10.155 ms; these are separate measurement series and are not pooled here.
Routing tuning alone did not improve the two shorter baseline cells.

`half_binding.py` reconstructs the accepted K18/K20 unions, checks the exact
IMT map columns, authenticates measured sources and the retained binaries,
and validates all timing rows. `HALF_PERFORMANCE_VERIFIED.json` records the
result. The K16 row is an implementation measurement without a full distance
certificate. Generated receipts live under ignored `measurements/`.

```text
python -B workstreams/inner_design/finite_migration/quarter_binding.py
python -B workstreams/inner_design/finite_migration/half_binding.py
python -B workstreams/inner_design/finite_migration/verify_ladder.py --m 18 --q1 workstreams/inner_design/finite_migration/Q1_LOWER.json --sparse workstreams/inner_design/finite_migration/SPARSE_M18_v3.json --dense workstreams/inner_design/finite_migration/DENSE_M18_short_v3.json
```

The first two commands check retained evidence; they do not launch benchmarks.
The measurement runner requires the isolated build and its authenticated
source/binary manifests. Never run it alongside another benchmark.

```text
python -B workstreams/inner_design/finite_migration/retune_fixed_input.py --seed workstreams/inner_design/finite_migration/M16_CENTERS_v1.json --output workstreams/inner_design/finite_migration/M16_DIRECT_POINT_v1.json --verify
python -B workstreams/inner_design/finite_migration/constant_point.py --seed workstreams/inner_design/finite_migration/M16_DIRECT_POINT_v1.json --output workstreams/inner_design/finite_migration/M16_CONSTANT_POINT_v1.json --verify
```

## K16 closure and timing binding

`DENSE_M16_coupled_v1.json` covers Q=64,...,512 with 1,096 leaves. It uses
the coupled input-normalizer bound and passed the full 512-bit replay.
Its dense-only margin is 42.3171086801 bits. The independently reconstructed
Q1 and Q2..63 contributions have margins 43.5928335021 and 60.9887506349 bits.
Their exact rational sum gives **41.8183267960 bits** in
`FULL_M16_VERIFIED.json`. No occupancy is omitted.

`half_ladder_binding.py` repeats the existing source, correctness, binary,
map, and timing checks and adds the K16 proof. `HALF_LADDER_PERFORMANCE_VERIFIED.json`
binds the 0.523367, 2.153740, and 10.109862 ms IMT measurements to full
certificates at K16, K18, and K20. This reuses the authenticated serial series;
it does not claim that a new benchmark was run.

Next: complete the joint-tuned parameter curves, then integrate the finite
IMT construction, proof appendix, tables, and measurements into the manuscript.
The selected-certificate work no longer blocks that integration. The current
curve pilot has a bounded joint-mixed search; its predecessor's full numerical
union remained uninformative despite passing coverage and witness replay.
