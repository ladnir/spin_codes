# First full-range proof with two-bit packets

Lowering the distance target closes the proof without changing the encoder.
At K=2^20 and N=2^21, the current two-bit ensemble has relative distance
greater than 9.25%, except with setup probability less than 2^-49.11.
This is a statement about a sampled code working for every nonzero message,
not a separate fresh-code experiment for each message.

The construction remains the one defined in [README.md](README.md):
BCH[256,128] outer rows, independently shuffled row coordinates, pairs of
adjacent rows, 256 independent regional packet permutations, independent
packet lane swaps, and IMT(128,19) with two independent transvections per
step. The expansion and feedback maps are the authenticated maps loaded
by the proof model. The state starts at zero and is not flushed.
The theorem concerns this ideal independent setup distribution.

No inner parameter or encoder implementation was changed to obtain this result.
Performance was not measured during this proof experiment. The subsequent
[tiled implementation](../TWO_BIT_IMPLEMENTATION.md) preserves the construction;
the theorem itself does not establish a 5 ms implementation.

## How the two parts cover every message

For a nonzero message, let q count pairs containing at least one nonzero
BCH input row. Every message belongs to exactly one q in {1,...,4096}.

For q<=400, use the earlier complete support covers at output-weight
threshold 209715. These bounds include every unequal support configuration
and every set of active pairs. Since 193986<209715, the same bounds apply
to the smaller bad event without rerunning those sparse calculations.

For q>=401, use the complete affine mixture cover at threshold 193986.
The positive row-mixture comparison dominates the averaged BCH counting
measure. Each cell sums all component assignments satisfying its coordinate
constraints and the active-pair threshold. Exact separation checks discard
only cells containing no feasible composition. The split-tree verifier
checks complete coverage; Arb recomputes each accepted cell's bound.

The complete points below all use the same construction. Each dense cover
has no unresolved cell and has been independently replayed at 256-bit
precision. Margins are rounded downward.

| Distance | Method | Visited cells | Terminal leaves | Dense margin (bits) |
|---:|---|---:|---:|---:|
| 5% | Lifted | 205 | 103 | 312.52286603460 |
| 6% | Lifted | 399 | 200 | 294.49211716425 |
| 7% | Lifted | 1451 | 726 | 80.04077802883 |
| 8% | Affine | 1079 | 540 | 86.79627001402 |
| 8.5% | Affine, retargeted from 8% | 854 new | 820 | 74.62992856911 |
| 9% | Affine, retargeted from 8% | 2256 new | 1453 | 73.84716062499 |
| 9.25% | Affine, retargeted from 9% | 2302 new | 2235 | 71.23535577341 |

These margins reflect sufficient witnesses, not optimized estimates of
the true failure probability. In particular, the 8% calculation uses a
stronger bound than the 7% calculation.

Add the sparse and dense expected bad-message counts. The conservative
sparse bounds in `closure.py`, together with each dense bound, give

    E[number of nonzero messages with output weight <= floor(delta*N)]
        < 1.642672754e-15 < 2^-49.11,
    delta in {0.05, 0.06, 0.07, 0.08, 0.085, 0.09, 0.0925}.

The probability that this count is at least one is at most its expectation.
For delta=0.0925, no such message implies minimum distance at least 193987,
which is greater than 0.0925*N. Thus the two occupancy ranges prove the stated
whole-code guarantee. No interpolation between sampled occupancies is used.

## Reproduction and verification scope

The sparse source bounds and their reproduction commands are in
[COVER_STATUS.md](COVER_STATUS.md). The new dense commands are:

```sh
python -B research/workstreams/permutation_locality/two_bit/affine_mixture.py --minimum-groups 401 --theta 2/5 --threshold 167772 --max-cells 2000 --max-depth 48 --output tmp/two-bit-affine-dense-401-d08.json
python -B research/workstreams/permutation_locality/two_bit/affine_mixture.py --retarget tmp/two-bit-affine-dense-401-d08.json --threshold 178257 --max-cells 1500 --max-depth 48 --output tmp/two-bit-affine-dense-401-d085.json
python -B research/workstreams/permutation_locality/two_bit/closure.py --dense tmp/two-bit-affine-dense-401-d085.json --precision 256
python -B research/workstreams/permutation_locality/two_bit/affine_mixture.py --retarget tmp/two-bit-affine-dense-401-d08.json --threshold 188743 --max-cells 2500 --max-depth 48 --output tmp/two-bit-affine-dense-401-d09.json
python -B research/workstreams/permutation_locality/two_bit/closure.py --dense tmp/two-bit-affine-dense-401-d09.json --precision 256
python -B research/workstreams/permutation_locality/two_bit/affine_mixture.py --retarget tmp/two-bit-affine-dense-401-d09.json --threshold 193986 --max-cells 2500 --max-depth 51 --output tmp/two-bit-affine-dense-401-d0925.json
python -B research/workstreams/permutation_locality/two_bit/closure.py --dense tmp/two-bit-affine-dense-401-d0925.json --precision 256
python -B -m unittest discover -s research/workstreams/permutation_locality/two_bit -p 'test_*.py'
```

`closure.py` independently replays the dense witness and recomputes the
final sum at 256-bit precision. It reuses the documented, previously
verified sparse bounds; it does not rerun every sparse cover. It rejects
incomplete dense covers, gaps between the two ranges, and cutoffs above
the sparse lemmas' scope. Cached upper bounds in witness files are ignored.
The unit tests include threshold, range, and safe retargeting checks,
plus the optional geometric and update-count diagnostics. A dense record
with three updates cannot be combined with the two-update sparse lemmas.

Generated witnesses stay in ignored `tmp/`; no production default changes.

## Next step

The 9.5% attempt is paused as of 2026-09-28. Sparse coverage reaches
q=424, but the q>=425 dense checkpoint still has 199 unresolved cells.
Finishing that cover and replaying the aggregate would be required for
a 9.5% claim. These saved calculations do not change the encoder.

The exact proved construction now has a
[tiled implementation and benchmark](../TWO_BIT_IMPLEMENTATION.md).
Use it as the control when comparing packet sizes and inner parameters.
Three updates have repaired a
selected 10% bottleneck and passed the
complete q=1,...,32 checks, but do not yet have a full certificate. See
[HILL_CLIMB.md](HILL_CLIMB.md) for the results and remaining obligations.
The earlier four-bit proof work remains separate and unchanged.
