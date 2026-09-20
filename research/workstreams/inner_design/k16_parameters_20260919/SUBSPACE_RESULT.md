# A certified t64, s12 inner for K = 2^16

Integration update: the reusable encoder now exposes this map as the explicit
K16-only `Configuration::T64S12` option in both libraries. See the
[integration record](../../spin_optimized/README.md#certified-k16-inner-option)
for forward/transpose tests and matched timings. The isolated measurements
below remain the original map-selection experiment; no default was changed.

## Outcome

The new expansion map gives a full **45.3819529593-bit bound** for half-rate
SPIN with the existing BCH [256,128] outer, K=65536, and bad-word cutoff 13107.
The union covers all 512 possible positive outer-row occupancies. All three
components passed 512-bit replay against their 256-bit outward bounds.

The measured transpose takes 0.332–0.335 ms with shared emission, versus
0.342–0.345 ms for the current (128,19) inner in the matched batch. This is a
2.0–3.7% improvement. The margin improves by about 3.56 bits over the current
full 41.818-bit result. The earlier 5% timing belonged to the different,
uncertified diagnostic (64,12) map; it is not this result.

No selected default, paper result, upstream Hypercat source, or submitted
artifact was changed. This is a verified experimental configuration for K16,
not a claim about the best parameters at other lengths.

## What changed

Relative to the first t64,s12 diagnostic candidate, only the expansion map
changes. Its feedback columns, transvection distribution, state dimension,
outer, and routing distributions stay fixed. Setup still uses independent
row and region permutations and one independent transvection per update.
The state starts at zero, persists across regions, and has no final flush.
Output is taken before the state update.

The old diagnostic expansion had three nonzero states with weight 16 and
one with weight 48. The new map has this exact nonzero spectrum:

| Weight | Multiplicity |
|---:|---:|
| 24 | 286 |
| 28 | 880 |
| 32 | 1743 |
| 36 | 912 |
| 40 | 274 |

There are 4095 nonzero states. Thus the expansion is a binary [64,12,24]
linear map, with no all-one output. Its minimum relative weight is 3/8,
matching the selected t128 expansion's 48/128 minimum.

`subspace.py` constructs the map from the existing t64,s15 diagnostic parent.
It preserves the six linear generators and constrains only the quadratic
coordinates. In the parent's state coordinates, the ordered new basis is:

```text
1, 2, 4, 8, 16, 32, 4160, 24704, 12544, 29184, 29696, 22528
```

The three constraints on the nine quadratic coordinates are 420, 473, and
227. Exhaustive enumeration verifies the embedding and the spectrum above.
The search also found an s13 map with weights in [24,40], but it has not
received the full proof/performance treatment. No s14 map was found in the
bounded search; that is not a nonexistence result.

At K16 a region has eight updates when t=64, versus four when t=128.
The smaller step gives more mixing opportunities. Choosing a better expansion
map also prevents a retained state from repeatedly emitting a weight-16 word.
The experiment shows that map selection materially affects the bound; it
does not isolate every contribution to the true distance-failure probability.

## Complete bound

| Contribution | Margin, bits |
|---|---:|
| One active outer row | 45.3838190316 |
| All occupancies 2 through 63, combined | 61.9776318483 |
| All occupancies 64 through 512, combined | 54.9887454409 |
| Exact sum of all three bounds | **45.3819529593** |

The BCH contribution uses the existing deterministic shell bounds, not a
guessed BCH spectrum. `assemble.py` rechecks their supporting evidence,
authenticates the component sources and replay receipts, verifies identical
instances, checks exact occupancy coverage, and sums rational upper bounds.
The resulting setup-failure bound is below 2^-40 for the stated cutoff.

The sparse cover uses eleven witnesses. Five were found automatically after
the initial point bank left gaps between its sampled occupancies. Every
accepted witness is reevaluated over the full sparse range with directed
positive folds and Arb terminal products. The dense cover has 1112 leaves,
after 33 refinements of the historical search geometry. Historical numerical
bounds were not transferred to the new inner.

The higher-density-band split and routing comparison are reused as outer-only
reductions. Their inner moments are rebuilt with t=64 and the new maps.
Q1 replay uses linear region products rather than repeated squaring. A finer
Q1 diagnostic gives 45.822 bits, but the full certificate uses the smaller
outward value in the table.

## Implementation and timing

The new map retains degree-at-most-two expansion, so it uses the existing
pruned zeta transform. Its finishing circuit has 29 XORs, compared with 30 for
the old diagnostic t64,s12 expansion. The optional shared feedback-emission
circuit has 136 XORs. Both circuits are symbolically checked against the
exact map columns before compilation.

The implementation retains the optimized four-row BCH kernel and direct K16
routing. No extra permutation, mixer round, or runtime polymorphism is added.
The measurements below use isolated generated research builds. The subsequent
public configuration selector and its separate measurements are linked above.

| Route seed | Current (128,19) | New map, plain emission | New map, shared emission |
|---:|---:|---:|---:|
| 1 | 0.342098 ms | 0.336658 ms | 0.335175 ms |
| 17 | 0.345194 ms | 0.333673 ms | 0.332420 ms |

Peach Ryzen 7950X, CPU 15, GCC 15.2.0, O3, x86-64-v3 baseline with znver4
scheduling; the selected BCH kernels are compiled separately for AVX-512.
Each cell is the median of three process medians. Each process times 101
in-place calls after three warmups, with setup and input initialization
excluded. There is no input reset/copy in the timing. The second repetition
reverses variant order. Benchmarks run serially under the shared locks.

Retained setup is 540672 bytes, versus 532480 for the control; workspace is
3 MiB for both. The extra setup stores the additional transvections.

Both emission variants pass the dense-oracle implementation tests, including
two routing layouts, in-place suffix preservation, compaction, alternate
routing seed, and boundary impulses. Both also pass ASan and UBSan with leak
detection enabled. Seven Python tests check the bound adapter, unchanged
historical globals, subspace linear algebra, exhaustive map embedding, and
exact spectrum. `bind_subspace.py` matches the full-bound columns to the
generated map manifest, header, and both measured remote header hashes.

## Reproduction and handoff

This is a research-checkout workflow; it uses the existing BCH bound assets
and historical witness files. Numerical evidence and generated builds remain
under ignored `measurements/`. Source and documentation can be committed
without those data files. No commit or push was performed.

1. `subspace.py` reconstructs the candidates; the selected one has s=12.
2. `subspace_check.py` produces Q1 and selected sparse witnesses.
3. `subspace_cover.py --mode sparse` fills gaps and sums q=2..63.
4. `subspace_cover.py --mode dense` rebuilds and refines the dense cover.
5. Repeat each proof producer with `--verify`, then run `assemble.py`.
6. `generate_map.py --maps <subspace-receipt>` generates the measured circuits.
7. `measure.sh` measures a fresh named batch; `sanitize_subspace.sh` checks both
   emission variants. `summary.py` checks hashes and reports matched timings.
8. `bind_subspace.py` matches the full proof's maps to the measured headers.

The retained successful receipts are `sub12_q1.json`,
`sub12_sparse_final.json`, `sub12_dense_cover.json`, their replay files, and
`sub12_full.json`. The measured generated map is
`measurements/subspace/Map64S12.h`. The timing batch is
`measurements/remote_subspace/subspace_v1`.

The reusable integration and forward/transpose adjoint tests are now complete.
Recommended next step: select the new option in K16 callers, retaining the
current K20 choice. Further shrinking the state is an optional later search,
not needed to obtain this improved certified operating point.
