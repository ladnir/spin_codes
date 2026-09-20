# IMT inner: BCH-256 migration

IMT means [Independent-Map Transvection](../../IMT.md). The name refers only
to the inner; the BCH outer and randomized bit-transpose permutation are
separate components of SPIN.

This work reuses the fixed BCH [256,128,d>=38] outer and replaces its inner
with the quarter-rate candidate's balanced expansion A, sparse feedback B,
and one transvection per epoch. The target is 10% relative distance with
40 margin bits at K=2^16, 2^18, and 2^20. The supported encoder is unchanged.

## Current result

Follow-up: [weight-five feedback](weight5/README.md) now has a full certificate
at K=2^20, 10% relative distance, and 50.0620882264 margin bits. Its exact
union and all 512-bit replays pass. Its isolated optimized transpose is now
tested and benchmarked; the two smaller sizes remain unresolved. The results below
refer to the original weight-three feedback map.

The Q=1 calculation passes at all three sizes. Q counts nonzero outer rows.
The new producer uses 256-bit outward arithmetic. Its 512-bit replay checks
the retained bounds and uses linear epoch iteration for Q=1.

| K | Q=1 margin | Contiguous sparse coverage | Full distance certified? |
|---|---:|---:|---|
| 2^16 | 43.59281 bits | Q=1..91 | No |
| 2^18 | 50.18916 bits | Q=1..255 | No |
| 2^20 | 50.06226 bits | Q=1..511 | No |

These are not full-code margins. Each partial union bounds only messages
whose number of nonzero outer rows lies in its covered interval. The event
is output weight at most floor(N/10), with N=2K. Probability is over the
shared setup of independent row permutations, region permutations, and
per-epoch transvections. State starts at zero, persists across regions, and
updates after output; there is no final flush.

The Q=1 calculation reuses the existing rational joint-shell objective by
coefficient domination. Higher occupancies use the existing certified shell
caps, not an assumed BCH spectrum. `verify_progress.py` authenticates those
receipts and reconstructs exact partial unions. It does not rerun the old LP
solvers or constitute an independent formalization of their proofs.

`SPARSE_Q64.json` evaluates every Q=2..64 separately with Arb. The batched
`SPARSE_M*.json` banks extend that coverage using directed positive folds
and Arb terminal powers. All accepted banks have 512-bit replay receipts.
The K16 range search stops at Q=92; it does not infer that Q=92 is impossible.

## What remains difficult

The first dense attempt reused the old K20 partition and replaced every
inner moment. After 24 witness retunings, 895 of 1,011 leaves remained weak.
Those records retain a complete partition, not a certificate of the target.

Searching with the independent-map Fourier moment fixed all weak sampled
points at K18 and K20 in the initial coarse grid. However, that grid missed
an intermediate-density bottleneck. A new adaptive K20 search produced 186
leaves, of which 88 remain unresolved. Every leaf was replayed at 512 bits,
including the weak leaves; numerical replay does not make a weak bound pass.

At Q=3482 and normalized density coordinate v=21/128, a retuned point bound
still has diagnostic margin about -1537 bits. Here the actual coordinate is
nu=p_min+(1-p_min)v, using the fixed outer-label probabilities. A negative
margin means this upper bound is useless there. It is not a lower bound on
failure probability or a counterexample to the code.
`diagnose_dense_gap.py` reproduces this witness and checks it at 512 bits;
`DENSE_GAP_POINT.json` records the point and its outward bound.

Further rectangle subdivision alone cannot repair this particular witness's
point value. A better witness family, tighter transfer, or tighter outer
label bound is needed. Do not keep refining the same cover without testing
the difficult points first.

## Completed overlap experiment and next step

The exhaustive audit is now complete. It improves the difficult point by
only about 1.68 bits. The exact zero-state contribution is the main limit
at that witness. Two weight-five feedback candidates pass the same point,
but neither has a full certificate. See `overlap/README.md` for the retained
results, scope of the obstruction, and reproduction limitations. The next
priority is bounded certificate searches for those candidate feedback maps.

The following describes the overlap tightening that was tested.

The current Fourier cap treats an overlap h between Aq and B^T a using only
their two weights. That allows the most unfavorable overlap separately for
every nonzero dual state a, even when the fixed maps exclude it.

For a fixed nonzero a, let D be the support of B^T a. Restrict the columns
of A to D. If that restriction has rank s, every nonzero q has at least one
one in D, so h>=1. A certified minimum distance d_D of that restricted image
would strengthen this to h>=d_D. A corresponding certificate on the
complement gives h<=wt(Aq)-d_complement. These are exact finite-map claims;
sampling restrictions is not sufficient for a universal Fourier bound.

The audit groups every nonzero dual state by weight and disjoint-basis
packing bounds on the two restrictions. The resulting point evaluations
are not a complete distance certificate. An explicit all-one-row split
remains a sparse-bound option for the short K16 instance.

## Reproduction and layout

Run from the repository root with the inner-design Python dependencies.
Do not use Python `-O`: assertions enforce certificate checks.

```text
python -m unittest discover -s workstreams/inner_design/asymmetric/bch256 -p "test_*.py" -v
python workstreams/inner_design/asymmetric/bch256/verify_progress.py
```

The fast verifier checks saved receipts and exact unions. To recompute the
numerics, replay the retained producers:

```text
python workstreams/inner_design/asymmetric/bch256/certify_bch_q1.py --verify
python workstreams/inner_design/asymmetric/bch256/sparse_bch.py --output workstreams/inner_design/asymmetric/bch256/SPARSE_Q64.json --verify
python workstreams/inner_design/asymmetric/bch256/sparse_ranges.py --output workstreams/inner_design/asymmetric/bch256/SPARSE_M16.json --verify
python workstreams/inner_design/asymmetric/bch256/sparse_ranges.py --output workstreams/inner_design/asymmetric/bch256/SPARSE_M18.json --verify
python workstreams/inner_design/asymmetric/bch256/sparse_ranges.py --output workstreams/inner_design/asymmetric/bch256/SPARSE_M20.json --verify
python workstreams/inner_design/asymmetric/bch256/adaptive_dense.py --output workstreams/inner_design/asymmetric/bch256/DENSE_ADAPTIVE_M20.json --m 20 --minimum 512 --verify
```

Replay receipt paths are write-once. For a fresh run without touching retained
artifacts, first produce to a new `--output` path, then replay that path.
The Q64 producer accepts `--occupations` followed by every integer 2 through
64. The range producer accepts `--m 16 --last 255`, `--m 18 --last 255`, or
`--m 20 --last 511`. The K16 producer is expected to save partial coverage.

`bch_model.py` parameterizes the independent-map transfer without modifying
the frozen quarter-rate producer. `dense_bch.py` inserts its moments into
the existing outer-only two-tilt reduction. `dense_search.py` searches new
witnesses. `adaptive_dense.py` retains the complete partition and explicit
unresolved leaves after a bounded run. `test_bch_adapters.py` tests region
coefficients, the alternate Q1 recurrence, directed folds, and cover geometry.

No performance benchmark or default-encoder change belongs to this result.
The remaining migration gates are listed in `../MIGRATION_PLAN.md`.
