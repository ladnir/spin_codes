# Inner calibration for quarter-rate SPIN

Status: 2026-09-11. Keep the BCH-derived [128,32,32] outer fixed, with
K=2^20 message bits and N=2^22 output bits. Compare RM2Sub maps for two
targets: 16.5% distance with 40-bit margin, and 19% distance with 30-bit margin.
The probability model is the same setup randomness used by the parent
[smaller-outer analysis](../SMALLER_OUTER.md). Every constituent here is fixed;
its spectrum is enumerated exactly, rather than replaced by an ensemble average.

## Current decision

The existing selected **t64_s16** map supports the first target with full
numerical occupation coverage. Keep **t128_s19** as the baseline supporting
both targets. No encoder implementation was changed, and no timing was measured.
All margins below use nearest binary64, not outward-certified arithmetic.

The calibration decision is now settled: retain t128_s19 and park t256.
The [subsequent outward certificates](../SMALLER_OUTWARD_CERTIFICATE.md)
prove both baseline targets and passed 512-bit replay. The table below
retains the numerical calibration values; alternatives were not promoted.

Update: the [dedicated t256 run](T256_DEDICATED_RUN.md) found a deterministic
distance ceiling. At t=256, s<=21 cannot reach 16.5%, and s<=24 cannot reach
19%, for any map or setup with this interface. The tested t256_s18 candidate
is ruled out, not merely waiting for a tighter numerical bound.

| Fixed inner | Margin at 16.5%, Q=1..128 | Margin at 19%, Q=1..128 | Full-coverage conclusion |
|---|---:|---:|---|
| Existing t128_s19 | 41.08356* | 30.05299 | Both targets supported numerically |
| Existing selected t64_s16 | 40.47606 | 29.44828 | 16.5%/40-bit target supported numerically |
| Prefix t128_s18 | 41.02992 | 30.00518 | Dense bound unresolved at both targets |
| Extension t256_s18 | 41.05405 | 30.03484 | Both targets impossible: deterministic ceiling 14.06403% |

*The original t128_s19 receipt uses Q=1..64 at 16.5%, followed by its dense
cover from Q=65. The table reports that retained sparse margin for this entry.

For t64_s16 at 16.5%, the dense contribution from Q=129 through all 32,768
outer rows has margin 1454.006 bits. Adding it to the sparse contribution
leaves **40.47605871881398 bits**. At 19%, this map's Q1 screen already misses
30 bits, and the transported dense witnesses also fail to establish the target.
A failed upper bound does not establish failure of the construction.

The two s18 candidates illustrate why Q1 screening is insufficient. Both
survive Q=1..128, but reusing the baseline dense witnesses does not close them.
For t128_s18, continuous witness optimization with three nodes per failing
root still leaves a useless aggregate dense bound. We have not exhausted
alternative maps, finer partitions, or stronger transfer bounds for t128_s18.
For t256_s18, the subsequent deterministic ceiling rules out both targets.

## What was varied

`calibrate.py` constructs nested state-dimension chains for t=64,128,256 and
s=14,...,20. They use the selected t64_s20, t128_s19, and t256_s14 maps as
anchors. Smaller dimensions take quadratic-row prefixes; larger dimensions
append independent quadratic rows. Two additional existing selected maps,
t64_s16 and t128_s15, provide controls. This is a 23-map screen, not an
optimization over every map at each parameter pair.

For each map, the producer checks rank, distinct nonzero columns, and BA=0.
It exhaustively enumerates A, computes the kernel spectrum by an exact integer
MacWilliams transform, and independently counts kernel weight-four words.
The retained selected t128_s19 map is checked against the existing baseline.

Map selection matters beyond s and t. The selected t64_s16 control has kernel
distance 6; the nested t64_s16 map has 16 kernel weight-four words. The nested
t128_s18 map has 64, and t256_s18 has 1280. These differences motivate more
deliberate map construction, but do not alone explain or prove the dense failures.

The sparse checker reuses exact-support Q1/Q2 bounds, composition bounds at
Q=3,4, and kernel-aware adaptive bounds through Q=128. The dense checker
transports the baseline's row-count boxes, but reevaluates every witness with
the candidate's t, s, and spectra. It checks the full integer coverage of the
transported boxes. It does not assume that an old numerical bound transfers.

## Engineering interpretation

t256_s18 uses 16,384 epochs and 294,912 field-row words at this output length.
t128_s19 uses 32,768 epochs and 622,592 field-row words. This motivated the
t256_s18 experiment; the subsequent ceiling rules it out for both targets.
Larger-state t256 choices retain the lower epoch count, but need new proofs
and timings. These counts are not speedup measurements.

t64_s16 reduces the state dimension but doubles the baseline's epoch count.
It uses 1,048,576 field-row words, so its smaller state need not make it faster.
The screen also records lookup and pruned-transform work indicators. They
exclude several costs and are not an encoder runtime model. Existing rate-half
whole-encoder timings must not be relabeled as quarter-rate measurements.

## Reproduction

Use the parent directory's Python dependencies and native Q2 helper. From the
repository root, run these commands in order:

```sh
python workstreams/rate_quarter_bch/inner_calibration/calibrate.py
python workstreams/rate_quarter_bch/inner_calibration/check_candidates.py --dense-only --tags t128_s19_nested t128_s18_nested t256_s18_nested t256_s16_nested t64_s16_selected
python workstreams/rate_quarter_bch/inner_calibration/check_candidates.py --tags t256_s18_nested t64_s16_selected t128_s18_nested
python workstreams/rate_quarter_bch/inner_calibration/refine_candidates.py --tag t128_s18_nested
python -m unittest discover -s workstreams/rate_quarter_bch/inner_calibration -p 'test_*.py'
```

`maps/` retains exact map data and constituent spectra. `Q1_SCREEN.json`
retains all 23 screens and cost indicators. Candidate `_full.json` receipts
contain all occupation ranges, including unsuccessful bounds; `full` denotes
coverage, not target attainment. `_dense.json` contains transported witnesses,
and `_dense_refined.json` contains the bounded refinement. Receipts bind their
inputs and producer code by SHA-256. The sparse search is reproduced by rerunning
its producer; regression tests check retained sparse union arithmetic rather
than independently repeating the entire search.

## Next step

Implement and benchmark the certified t128_s19 quarter-rate encoder.
The t256 investigation is parked because its larger required state makes
the performance benefit uncertain. If resumed, start with s>=22 for 16.5%
and s>=25 for 19%, applying obstruction tests before upper-bound searches.
