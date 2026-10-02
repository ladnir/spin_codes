# Packed GL32: First Whole-Code Certificate

2026-09-30. The canonical GL32 construction closes **>9.5% relative distance
with >36.68 bits of setup-failure margin** at K=2^20 and rate one-half.
The precomputed transposed encoder remains at **6.203224 ms** for 128-bit
XOR elements. No implementation change was needed to finish the proof.

## Construction

The stable ensemble identifier is `canonical-gl32-width8-shared4-r2`.
The forward encoder maps K=2^20 bits to N=2^21 bits as follows.

1. Encode 8,192 message rows with the fixed BCH[256,128] constituent.
   Collect each four adjacent rows into one of 2,048 groups.
2. Partition each group into 32 canonical blocks of eight adjacent columns.
   Apply an independent uniform GL(32,2) map to each four-row/eight-column block.
3. Independently shuffle the 256 columns of each group, sharing that shuffle
   across its four rows. Each column becomes one four-bit packet in its region.
4. Independently shuffle the 2,048 packet positions within each of the
   256 regions, then concatenate the regions.
5. Apply IMT with t=128, s=19, two independent transvections per step,
   zero initial state, and no terminal flush.

The fixed expansion and feedback maps are `Map128S19` in
[the selected maps](../../../../spin/src/kernels/generated/SelectedMaps.h).
For input block X_i and entering state Q_i, the forward recurrence is

    Y_i = X_i + A Q_i,
    Q_(i+1) = M_(i,2) M_(i,1) Q_i + C X_i,     Q_0 = 0.

Here A maps 19 state bits to 128 output bits, and C maps 128 input bits
to 19 state bits. Each M is I+u v^T, with u uniform nonzero and v uniform
in u's orthogonal complement. The transvections and all preceding setup
choices are independent. Setup is sampled once and fixed for every message.
There is **no additional GF16 packet randomizer**.

For every fixed setup, the outer is injective and the local maps and
permutations are invertible. The inner is also invertible: recover X_i
from Y_i and Q_i, then compute Q_(i+1). Thus every realization has rate
exactly one-half. The transposed implementation applies the adjoint map;
the distance statement below concerns the forward code.

## Certified Claim

Let Z count nonzero messages whose encoded word has weight at most 199229.
The independent 384-bit replay gives the following first-moment bounds.
Occupancy q counts nonzero four-row groups.

| Scope | Coverage | Negative log2 of outward upper bound |
|---|---|---:|
| Sparse | Every q=1..32 | 37.959690941543 |
| Dense | Every q=33..2048; complete 240-interval partition | 37.450099737024 |
| Whole code | Exact sum of all fresh contributions | 36.682511508003 |

The whole-code outward endpoint U is approximately 9.0670e-12 and satisfies
U < 2^-36.68. This last comparison was also checked with exact rationals:
U^25 < 2^-917. Therefore, over the stated ideal setup distribution,

    Pr[d_min <= 199229] <= E[Z] <= U < 2^-36.68.

Outside this bad-setup event, the minimum weight is at least 199230, giving
relative distance strictly greater than 9.5%. The result exceeds the
requested 20-bit whole-code margin. It covers every nonzero message, not
only selected occupancies or comparison means.

## Why the Local Mixing Helps

For a fixed nonzero 32-bit input, uniform GL32 produces a uniform nonzero
32-bit output. Viewed as eight packets, its support W satisfies

    Pr[W=w] = binom(8,w) 15^w / (2^32-1),     1 <= w <= 8.

Conditional on the support, the nonzero packet labels are independent and
uniform. This removes the row-rank dependence of the earlier shared route.
It does not activate an empty canonical block. The proof therefore uses
authenticated BCH shortening bounds to count tuples supported on few blocks.
The [local proof note](LOCAL_KERNELS.md) gives the counting and transport
argument, including why independent group setup permits multiplying the
averaged counting measures.

The numerical closure uses these transported counts and a checked positive
comparison mixture. The sparse checker includes packet supports 5 through
256. The dense checker covers its entire auxiliary mean domain. Six final
depth-ten intervals needed subdivision; all twelve children passed. Neither
the construction nor the distance target changed to close those intervals.

## Replay and Retained Evidence

All receipts below are local, ignored numerical artifacts; they are not raw
data to commit. The [progress ledger](PROGRESS.md) records their search
provenance, earlier checkpoints, and how to regenerate the witnesses.

The complete dense search is
`tmp/packed-closure/r2-d095-complete-search-p256.json`, SHA256
`d99341614abba8f7f35074e2dcfaa71062ec0201cc9f5dd66dca1d5ee3fb48cc`.

The final receipts are under `tmp/packed-closure/gl32-d095-whole-p384/`:

| Receipt | SHA256 |
|---|---|
| `whole.json` | `a10e72479ce480839a461aff5b096ee4f0a3e828c705c6aabfc4657c6ad4b67e` |
| `sparse.json` | `d5c939d1b89ff4410f67a9e25fcec0b88b27ea3864b61b4d920c52ec5fde6a7f` |
| `dense.json` | `9896b5d2db5c65f8dd8a9e6ae1c5eb8ebe3a507248dc4b9672c6ca246a7c63b5` |

Run from the repository root, choosing a new output directory:

```text
python -B research/workstreams/permutation_locality/packed_mixing/whole_replay.py --dense tmp/packed-closure/r2-d095-complete-search-p256.json --sparse tmp/packed-closure/sparse-q1-d096.json tmp/packed-closure/sparse-q2-d096-smalltilts.json tmp/packed-closure/sparse-q2-4-d096.json tmp/packed-closure/sparse-q5-7-d096-smalltilts.json tmp/packed-closure/sparse-q5-32-d096.json --output-dir tmp/packed-closure/gl32-d095-whole-p384-rerun --precision 384 --target-bits 20 --dense-workers 4
```

The sparse and dense stages run sequentially in fresh processes. Four
uncached workers evaluate the dense witnesses. The checker rebuilds the
BCH premises and inner operators; saved search endpoints are not proof
inputs. Exact dyadic summation precedes one upward rounding of the final
endpoint. All 154 local tests pass.
A separate receipt audit verified all eight referenced file hashes, exact
scope and partition coverage, and the sum against the final outward endpoint.

## Performance and Next Step

The [matched performance report](../packed_driver_REPORT.md) retains the
6.203224 ms result. The same-build optimized shared-GF16 control takes
5.805222 ms, so GL32 costs 0.398002 ms for the stronger certified distance.
Setup and allocation are outside these timings. No new benchmark was run
during closure, and production defaults and paper claims remain unchanged.

Keep this as the proved baseline. The next proof experiment should reoptimize
near comparison means .104, .114, .120, and .208 before attempting 9.6%.
Simply retargeting the present witnesses fails; the search stopped when
each interval met its budget, so that failure does not establish a distance
ceiling for the encoder. Further implementation tuning can separately target
the remaining roughly 0.2 ms above 6 ms.
