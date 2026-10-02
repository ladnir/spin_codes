# Packed GL32 Closure Progress

## Target and Status

2026-09-30. Target: a whole-code first-moment certificate for relative
distance greater than 9.5%, with at least 20 bits of setup-failure margin,
near 6 ms of precomputed transposed encoding time. **R2 closed:** fresh 384-bit
whole-code replay certifies **>9.5% distance with >36.68 bits of margin**.
The computed whole-code margin is 36.682511508003 bits. See the
[closure record](FIRST_CLOSURE.md) for the exact claim and reproduction.
Existing certificates in the [proof index](../PROOF_INDEX.md) are unchanged.

The follow-on target is **10% distance / 40 bits**; see the separate
[hill-climb ledger](HILL_CLIMB.md). The R4 variant has passed fresh 384-bit
replay of every sparse occupancy and all 337 dense intervals. Independent
exact-dyadic aggregation gives the following margins:

| R4 fresh replay scope | Audited margin (bits) |
|---|---:|
| Sparse q=1..32 | 58.0241331220138 |
| Dense q=33..2048 | 52.0760654522247 |
| Sum | 52.0528843554777 |

**R4 is now closed at >10% distance with >52.05 bits of setup-failure margin.**
The coordinator in `tmp/packed-hill/r4-d10-whole-parallel-p384/` completed
its exact aggregation and source-provenance checks. Its `whole.json` SHA256 is
`637a1ac049a7aa81441705e5a4fbeebf79e78fa80e63fe54f93cfe6dad05a8f0`.
The outward endpoint U satisfies U^20 < 2^-1041 in exact integer arithmetic.
The cutoff is 209715, so successful setup gives minimum distance at least
209716 for N=2097152, strictly greater than 10%. Failure has probability at
most U under the declared independent ideal setup, sampled once and fixed
for every message. The implementation's PRNG is not part of this claim.
See [the R4 closure record](R4_CLOSURE.md). All **260 packed-mixing tests pass**.

The [fused R4 implementation](../packed_fused_r4_REPORT.md) measures
**6.212096 ms**, versus **6.753111 ms** for sequential R4 in the same build.
Fusion preserves every sampled map; it is not a reduced-randomness variant.
This follow-on does not overwrite the completed R2 result below.

The retained R2 construction is K=2^20, N=2^21, BCH[256,128], four rows per group,
independent uniform GL(32,2) maps on each canonical four-row/eight-column
block, independent shared column permutations per group, independent
regional packet shuffles, and IMT(128,19) with two updates. There are 2,048
groups. GL32 supplies the conditional independent nonzero packet labels;
there is no additional GF16 randomizer.

## Retained R2 Implementation Checkpoint

The [optimized-driver report](../packed_driver_REPORT.md) establishes a
stable **6.203224 ms** median, versus **5.805222 ms** for the same-build
legacy encoder. The difference is +0.398002 ms / +6.86%. The retained old
binary remeasures at 5.829211 ms. Two seeds, both process orderings, and
101 measured calls per process were used, with all benchmarks serial.

The large earlier harness outlined callbacks in the inner loop. The new
fixed-width emitter forces inlining, with disassembly checked before timing.
Independent full-reference, adjoint, and preserved-suffix checks pass at
K=2^14 and K=2^20 for both seeds. The production encoder is unchanged.

## Retained R2 Proof Checkpoints

The BCH premises are regenerated and authenticated, including exact LP and
positive-polynomial witnesses, dual distance, and containment. The original
canonical block counts are transported to rational expected support counts.
The sparse interface now includes supports 5 through 256, including the rare
tail below the old BCH support floor of 38.

1. **Local and interface tests:** 154 tests passed at this R2 checkpoint, including exhaustive toy
   spectra, direct-shell transport, exact rational tails, support-partition
   geometry, replay scope/degree checks, witness reuse, cache invalidation,
   and whole-code aggregation rejection cases. The whole-code tests use toy
   and mocked replay outputs; they are not a certificate for this encoder.
2. **Sparse occupancy:** q=1,2,3,4 are bounded at 9.6% distance with margins
   37.05, 61.31, 86.32, and 184.00 bits, respectively. The initial q=2
   obstruction disappeared after adding smaller output tilts; no construction
   change was needed. The complete q=1..32 search now passes at 9.6%; the
   exact sum of its 32 outward dyadic bounds has **35.19124 bits** of margin.
   Fresh 384-bit replay at 9.5%, cutoff 199229, gives **37.95969094 bits** for
   the entire sparse prefix. All 32 occupancies, source hashes, and dyadic
   aggregation checks pass. This number bounds only the sparse prefix.
3. **Dense selected points:** the table below covers selected comparison
   means for q>=33 at 9.5%, not intervals between them.
4. **Dense exhaustive cover:** the bounded first batch finished with
   direct expected-shell bounds and 42 exact positive comparison components.
   Every accepted interval is checked with outward arithmetic. After 120
   visited cells there are 22 accepted leaves and 77 pending cells. The
   accepted leaves cover 78.125% of the comparison interval, not of the
   message set or its failure probability. The remaining interval is
   approximately [0.0002121, 0.2189157], partitioned into 35 depth-eight and
   42 depth-nine cells. No cell exhausted the configured depth limit.
   Fresh 384-bit replay checks all 22 accepted leaves; their combined
   log2 upper bound is -480.565698. This excludes every pending cell and
   therefore does not bound the full dense range.
5. **Narrow-cell continuation:** the pending cells were split geometrically
   to width at most 1/1024, where the stronger regional-count bound is
   enabled. This preserved all accepted leaves and bypassed 147 broad
   internal cells without evaluating them; splitting itself proves nothing. The next 20 cells
   all passed fresh 256-bit outward checks. The saved 140-visit checkpoint
   has **42 accepted leaves and 204 pending depth-ten cells**, covering
   **80.078125%** of the comparison interval. Its pending range is
   approximately [0.01973918, 0.21891568]. No depth-limit failure occurred.
   At this stage, the independent 384-bit replay required for a final claim
   had not run. The earlier 22-leaf replay remains separately retained.
6. **Seeded continuation:** after freshly rechecking all 42 saved leaves,
   the next 20 new intervals also passed. The retained 160-visit snapshot
   has **62 accepted leaves and 184 pending depth-ten cells**, covering
   **82.03125%** of the auxiliary comparison interval. The remaining range
   is approximately [0.03926629, 0.21891568]. No new cell required splitting
   or exhausted the depth limit. Of these 20 checks, 17 accepted a reused
   proposal, including 16 neighboring-cell witnesses; three used ordinary
   optimization. The inner-operator cache recorded 16 hits. These are search
   reuse counts, not an encoding-speed measurement or a whole-code bound.
7. **Second seeded checkpoint:** the next 20 narrow intervals also passed,
   bringing the saved 180-visit checkpoint to **82 accepted leaves and 164
   pending depth-ten cells**. Coverage is **83.984375%** of the comparison
   interval, and the pending range is approximately [0.05879340, 0.21891568].
   Across the two seeded batches, 34 of 40 intervals accepted reused
   proposals; 33 used neighboring-cell witnesses. No narrow interval in
   these batches required further splitting. The run was stopped at its
   checkpoint to try wider intervals; its completed work is 40 new cells,
   not the 80-cell budget in its filename.
8. **Parallel replay regression:** three fresh, uncached worker processes
   replayed the original 22 accepted intervals at 384-bit precision. Every
   path, interval, exact dyadic upper endpoint, and the aggregate matched
   the retained sequential replay exactly. The comparison used the same
   authenticated mixture and outer premises. Both replays still exclude
   77 unresolved cells; this validates the parallel replay implementation,
   not full dense coverage.
9. **First wider checkpoint:** after rechecking and retaining all 82 leaves,
   the next 20 candidate evaluations accepted four width-1/256 intervals
   and subdivided 16 others. The retained 200-visit snapshot has **86
   accepted intervals and 53 pending cells** (21 at depth eight, 32 at
   depth nine), with **85.546875%** auxiliary comparison-interval coverage.
   Its pending range is approximately [0.07441508, 0.21891568]. No cell hit
   the depth limit, and ordinary fallback optimization was not invoked.
   Three of the four new accepted intervals used warm neighboring witnesses.
   A failed wide hint requests a narrower check; it is not a low-distance
   codeword or a failed narrow-cell certificate.
10. **Return to narrow intervals:** the broad hints eventually requested 111
   subdivisions. They did not improve coverage beyond the four accepted
   broad cells. Once width 1/1024 was reached, the next 25 intervals all
   passed, bringing the saved 320-visit checkpoint to **111 accepted
   intervals and 123 pending depth-ten cells**. Coverage is **87.98828125%**.
   The serial process was stopped with this saved checkpoint preserved;
   unfinished work after it was not retained.
11. **Parallel continuation:** three fresh processes each receive 41
   disjoint pending roots from that immutable checkpoint. They authenticate
   the BCH and mixture premises independently and search only their assigned
   subtrees, without replaying unrelated accepted leaves. The global model
   root and numerical evaluator are unchanged. A separate merger verifies
   exact coverage and provenance, discards worker numeric endpoints, and
   writes search witnesses only. Final independent replay is still required.
   A real cold-start test guards the local/legacy `dense_cover` namespace
   collision found before the workers began their numerical searches.
12. **Complete witness cover:** the three workers finished with 47, 41, and
   41 accepted leaves, respectively. Six of the first worker's depth-ten
   intervals needed one subdivision; all twelve children passed. Exact
   merging with the 111 retained leaves gives **240 accepted intervals and
   no unresolved cells**. The final search has full comparison-domain
   coverage; the following replay establishes the whole-code probability bound.
13. **Whole-code closure:** fresh 384-bit replay passed all 32 sparse
   occupancies and all 240 dense intervals, including 103 regional witnesses.
   The dense stage used four independent uncached workers. The exact fresh
   sum has **36.682511508003 bits of margin**: sparse 37.959690941543 bits
   and dense 37.450099737024 bits. Every sparse row and its aggregate match
   the earlier 384-bit replay exactly. The final outward endpoint also passes
   the exact rational check U^25 < 2^-917, proving margin greater than 36.68.
   The construction, cutoff 199229, and implementation did not change.
   A separate receipt audit checked all eight referenced source/fresh hashes,
   the exact partition, complementary scopes, and the exact aggregate.
14. **Search-cost diagnostic:** the 58 regional witnesses among the 82
   retained leaves contain 928 variance parts and 29,328 MGF candidates.
   The ordinary count-cap loop performs 60,093,072 candidate/count updates
   across those witnesses. The isolated `count_caps_fast.py` instead uses
   floating scores only to select a candidate for each count, then computes
   that candidate with outward arithmetic. Every candidate bounds every
   count; a poor selection can weaken the result but cannot understate it.
   All counts and the original baseline remain. The helper is not wired
   into the search or final replay. On saved leaf `0000111011`, all 16
   variance parts and 2,049 counts agreed with the ordinary evaluator to
   relative difference at most 3.20e-74 at 256-bit precision. The diagnostic
   took 5.30 seconds for the ordinary count caps and 3.22 seconds for the
   helper. These are single-run checker-component timings, not an encoder
   benchmark or an end-to-end checker speedup. Full-suite testing also
   caught and fixed an import-path collision; a cold-import test now guards it.

A fresh read-only review found no blocker in the expected-count bridge:
independent group setup permits multiplying averaged counting measures;
the comparison preserves active-group labels even on artificial zero
draws; and the sparse folding includes support 5 through 256. This review
did not rederive every inherited inner bound or establish the missing
numerical coverage.

| Comparison mean | Outward log2 first-moment upper bound at 9.5% |
|---|---:|
| 0.001 | -1246.57 |
| 0.008 | -277.90 |
| 0.032 | -409.82 |
| 0.064 | -827.62 |
| 0.128 | -929.19 |
| 0.25 | -7938.52 |
| 0.5 | -51491.43 |
| 0.75 | -88468.82 |

The 0.032 row is an outward exact-threshold rescaling of the earlier checked
10% witness. The other seven rows are fresh 256-bit point evaluations using
the expected-CDF shell majorant. The exhaustive cover uses the tighter
direct-shell majorant, not these saved numerical answers. A point bound
cannot establish a whole-code certificate, regardless of its size.

## What Still Constitutes a Proof

The sparse and dense scopes must together cover every q=1..2048. All saved
witnesses must be freshly replayed from authenticated BCH counts and inner
operators. Exact cutoff, construction, and scope metadata must agree.
Finally, the sum of the sparse and dense dyadic upper endpoints must be less
than 2^-20. Per-cell or per-occupancy targets alone do not satisfy that test.

The 9.6% sparse witnesses may be reused at the smaller 9.5% cutoff only after
the stated monotone retargeting/replay. Fresh replay does not treat saved
floating optimizer values or aggregate endpoints as proof inputs.

## Reproduction and Retained Receipts

Run local tests:

```text
python -B -m unittest discover -s research/workstreams/permutation_locality/packed_mixing -p "test_*.py" -v
```

The initial bounded dense cover was run as follows (a reproduction needs a
fresh output path):

```text
python -B research/workstreams/permutation_locality/packed_mixing/dense_cover.py --distance .095 --comparison shell --max-cells 120 --target-bits 36 --precision 256 --output tmp/packed-closure/r2-d095-shell-cover-120-p256.json
```

Use `--resume <prior-receipt>` with a new output path to continue; accepted
leaves are recomputed before retention. Use `--replay <prior-receipt>` for
a fresh bound on its accepted leaves. A replay explicitly reports remaining
unresolved cells and does not turn a partial cover into a complete one.
The first dense checkpoint's SHA256 is
`c56c153ca005b14cd98d9814a3bc9d470a494fa0b584c254a74d64e741bde7a4`.
Its partial fresh replay is
`tmp/packed-closure/r2-d095-shell-cover-120-replay-p384.json`, SHA256
`edb926636ea7694c4648617365ae8e054f200a29e8b6299741b2413dfa268f0f`.
Both record 77 unresolved cells. All jobs for this checkpoint have finished.
The independently rebuilt parallel replay is
`tmp/packed-closure/r2-d095-shell-cover-120-parallel-replay-p384.json`, SHA256
`357d0ebd9199370188f0a8a86972861b61a70ee6190e8c1f63e1fd36b451da1c`.
Its invocation adds `--replay-workers 3` to the ordinary dense replay.
Worker logs are retained beside the receipt in its `-worker-logs` directory.

The narrow-cell continuation is
`tmp/packed-closure/r2-d095-regional-presplit-320-p256.json`, SHA256
`70d1a5273df42cf2479f2ddfa9163794bae3029646afd5c52d2b62d65354daff`.
Despite the filename's 320-cell work budget, this run was stopped after its
first 20-cell checkpoint to switch search strategies. Its actual cumulative
visit count is 140; do not infer completed work from its filename.

The continuation search adds `--point-seeds` for the seven 9.5% point
witnesses and the earlier mean-0.032 witness, plus `--region-cache`.
Successful neighboring-cell witnesses are also reused as proposals. The
cached regional polynomial depends on the actual inner, output tilt,
precision, and local refinements, not on the mean cell or outer counts.
Those cell-dependent terms are recomputed. Reuse retains only one operator
and one neighboring witness in memory; it does not import numerical bounds
from old receipts. The final replay rejects all these search options.

The first seeded continuation was launched with an 80-cell budget:

```text
python -B research/workstreams/permutation_locality/packed_mixing/dense_cover.py --resume tmp/packed-closure/r2-d095-regional-presplit-320-p256.json --point-seeds tmp/packed-closure/r2-d095-selected-p256.json tmp/packed-mixer-results/canonical-r2-d10-m032-refined-p256.json --region-cache --max-cells 80 --target-bits 36 --precision 256 --output tmp/packed-closure/r2-d095-seeded-cached-80-p256.json
```

Use a new output path when rerunning. The output's `cover.visited`, accepted
leaves, pending cells, and reuse counters describe completed work; the
configured work budget and mere existence of an output file do not.

Its first 20-cell checkpoint is retained separately as
`tmp/packed-closure/r2-d095-seeded-checkpoint160-p256.json`, SHA256
`edb4c21c9720f2f5b0d785780ce999525615389578b2b81fa13ada1c6fbd50aa`.
The original 80-cell process subsequently reached visit 180 and was stopped
at that checkpoint. Its final receipt is the `r2-d095-seeded-cached-80-p256.json`
path above, SHA256
`20604f3ad3f4d4b01e88141681ea4aabb2aa34e44bba6e9c8acf8f59608b0b44`.
The separate 160-visit snapshot remains unchanged.

The next continuation uses `--seed-width 1/256`, `--region-cache`, the same
point sources, and a 320-cell work budget. Its output is
`tmp/packed-closure/r2-d095-wide-warm-320-p256.json`. The geometric
preprocessing preserves all 82 accepted cells and turns the 164 unresolved
cells into 41 wider candidates. It performs no numerical evaluations and
does not add coverage. Failed wide hints request subdivision; ordinary
fallback optimization starts only at width 1/1024. Successful warm hints
skip floating screening but still get a fresh exact bound.
The fresh recheck retained all 82 saved leaves before new-cell search began.
The first new checkpoint is retained separately as
`tmp/packed-closure/r2-d095-wide-warm-checkpoint200-p256.json`, SHA256
`62b55c113399e16eba6c9922e88f8043c4fd25324fd70f7b86ab37dd946b0e5a`.
It records 86 accepted intervals and 53 pending cells. The original
320-cell-budget process was later stopped at its saved visit-320 checkpoint.
Its immutable parallel base is
`tmp/packed-closure/r2-d095-parallel-base-p256.json`, SHA256
`2c1cbc4d877a3a2c902cd8dc19d9931b561854fd55c239479402f247b1365f9a`.
This base has 111 accepted intervals and 123 pending cells. Each worker uses:

```text
python -B research/workstreams/permutation_locality/packed_mixing/cell_search.py tmp/packed-closure/r2-d095-parallel-base-p256.json --shard-index I --shard-count 3 --precision 256 --target-bits 36 --max-cells 160 --checkpoint-every 5 --output tmp/packed-closure/r2-d095-shardI-p256.json
```

Here `I` is 0, 1, or 2. These are proof-search processes, not concurrent
encoder benchmarks. Once all workers finish, use `cell_merge.py --source
<base> --partials <three-results> --output <new-dense-path>`. A merged cover
may still contain unresolved descendants; the merger does not certify a bound.

The completed worker receipts have SHA256 values:

| Worker | SHA256 |
|---|---|
| 0 | `a4b033d6d1998a3e747c0d395da8baf70011689b7faff1963e04017bca6dfd1c` |
| 1 | `89942e0d7c9391bca3afe75c8137d7bf3a31621f6bedff513c211415c55629f0` |
| 2 | `f3c68d10a9a2e7a85dee7726519d09ecccee9de3439267f425a54c48a4538a11` |

Their merged search is `tmp/packed-closure/r2-d095-complete-search-p256.json`,
SHA256 `d99341614abba8f7f35074e2dcfaa71062ec0201cc9f5dd66dca1d5ee3fb48cc`.
It has 240 leaves and zero unresolved cells. The merger discards the new
workers' numerical bounds. This file supplies witnesses to final replay,
not a claimed probability endpoint.

The fresh sparse receipt is
`tmp/packed-closure/sparse-q1-32-d095-replay-p384.json`, SHA256
`5c1f8b954d874ef82b3f498ca1652d8b38aaa2fce4d30a94ad6462e4273d43e6`.
It replays the passed, disjoint scopes from `sparse-q1-d096.json`,
`sparse-q2-d096-smalltilts.json`, `sparse-q2-4-d096.json`,
`sparse-q5-7-d096-smalltilts.json`, and `sparse-q5-32-d096.json` in the same
directory. Use `sparse.py --replay <these-five-paths> --distance .095
--precision 384 --target-bits 24 --output <fresh-path>` to reconstruct it.
The per-occupancy replay target is not the aggregate certificate margin.

Local receipts are under ignored `tmp/packed-closure/`: the selected-point
run, the dense checkpoint, and separate sparse searches. Performance logs
are under ignored `tmp/packed-driver-results/`. They are not raw data to
commit. The sources, proof interfaces, and this progress record are the
durable work product.

## Next Decision

Preserve the closed R2 9.5% / 36.68-bit / 6.203224 ms baseline and the new
**R4 >10% / >52.05-bit / 6.212096 ms** result. The
[hill-climb ledger](HILL_CLIMB.md) records the completed whole replay and
matched exact-map fusion. Retain the official receipt and independently
audited component endpoints. Next, target the packed GL32/BCH phase using
the fused R4 baseline, keeping the sampled distribution unchanged.

For a future return to the historical R2 9.6% plan, first reoptimize point
or narrow-cell bounds near comparison means
.104, .114, .120, and .208. Reusing the current witnesses unchanged is
insufficient: several near-budget bounds lose hundreds of bits when the
cutoff increases by 2,097. Search early stopping means this is not evidence
of a distance ceiling. Keep further speed tuning separate from that test.
The count-cap helper remains an isolated prototype; it is not enabled in
the continuation or final replay.

The successful `whole_replay.py` invocation ran fresh 384-bit sparse and
dense replays in `tmp/packed-closure/gl32-d095-whole-p384/`. Its `whole.json`
SHA256 is `a10e72479ce480839a461aff5b096ee4f0a3e828c705c6aabfc4657c6ad4b67e`.
The [closure record](FIRST_CLOSURE.md) gives all three final hashes and the
full invocation. All proof processes for that R2 closure have finished.

For another replay, `whole_replay.py` runs fresh 384-bit sparse and
dense replays in separate processes and sums every freshly computed endpoint.
It checks the complementary scopes q=1..32 and q=33..2048, the exact cutoff,
construction, count premises, and source hashes. Only a sum below 2^-20
produces `whole.json`; a partial dense cover is rejected before replay.
Pass the completed dense receipt to `--dense`, the five original sparse
search receipts listed above to `--sparse`, and a new directory to
`--output-dir`. Do not pass the aggregate sparse replay as though it were an
original search witness file.
The optional `--dense-workers 4` runs only the dense stage in parallel;
the sparse and dense stages still run sequentially. Each dense worker
independently rebuilds the inner and uses the ordinary uncached evaluator.
One worker remains the default. Missing results, scope mismatches, and
worker failures prevent a final certificate.
