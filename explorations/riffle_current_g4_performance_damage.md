# Performance of the current proof-oriented g=4 construction

## Question

The `g=4` certificate includes more setup randomness than the intended packet-only deployment. This probe measures the complete transposed map on `GF(2^128)` elements.

The probe implements the proof-oriented construction without a whole-word permutation pass. Every permutation is fused into an address schedule or a local operation.

## Implemented transpose

The parameters are

```text
K = 2^20
N = 2^21
packet width = 4
inner nodes = N/64 = 32768
```

The implementation applies these operations in order.

1. It traverses the recursive inner from right to left.
2. Each node gathers sixteen independently permuted four-element packets.
3. The generated circuit computes `Acc^T(Y + P^T S)` in 930 field XORs.
4. One uniform 64-coordinate permutation maps the resulting adjoint state to the preceding node.
5. The outer transpose reads the three physical bands sequentially.
6. The tile read fuses the per-group lane bijection and each block's within-band coordinate permutations.
7. Three generated partial BCH circuits accumulate each data block directly at its assigned output address.
8. The graph transpose gathers the 128 replacement positions, applies one BCH transpose, and adds `R^T` using three byte tables.

The packet permutation moves no values. The state permutation touches one hot 64-element array per node. The outer permutations introduce no `N`-element workspace.

The benchmark samples a uniform packet permutation, a uniform data-block assignment, and every declared local permutation during setup. It also implements the global-lane puncture rule.

## Correctness gate

The optimized inner was compared with a direct implementation of `P^T`, the accumulator transpose, and every state permutation. The optimized outer was compared with 16,384 independent full BCH transposes and direct address reconstruction.

The comparison covered all `N` input elements and all `K` output elements. It included the 128 graph replacements and the dense 24-by-`K` graph map. The comparison passed.

## Matched measurement

The executable used one pinned Windows thread. It was compiled by MSVC with `/O2`, `/arch:AVX2`, and C++20. The timed variants ran serially.

The matched run used 31 trials for each executable. Times are medians in milliseconds.

| Operation | Time |
| --- | ---: |
| Four-element packet chain, accumulator, no state permutations | 10.861 |
| Current inner, including state permutations | 12.267 |
| Fixed three-band outer without random maps | 7.766 |
| Current randomized three-band data outer | 11.259 |
| Current inner inside the complete path | 12.514 |
| Current data outer inside the complete path | 12.170 |
| Graph BCH and dense `R^T` | 1.291 |
| Complete current transposed encoder | **26.001** |

The 31-trial end-to-end median was 26.001 ms. The executable prints the complete sample list.

The current hot path reads about 13.06 MiB of setup metadata. An additional 4 MiB stores inverse maps used only by the independent correctness check.

## Historical implementation

The preserved historical source was rebuilt against the same `cryptoTools` library. The rebuild used the same compiler flags and thread-affinity code as the current probe.

The historical inner has the following structure.

1. It partitions the input into 32,768 contiguous groups of 64 field elements.
2. It samples one permutation of those groups.
3. It traverses the permuted groups from right to left.
4. It retains one 64-element state in cache.
5. A generated 913-XOR circuit applies the cyclic split BCH transpose.
6. The loop prefetches every cache line of the group used two iterations later.

The historical construction has no four-element packet permutation. It has no per-node state permutation. Its sequential outer avoids the current lane and coordinate maps.

The historical measurements were:

| Historical operation | Time |
| --- | ---: |
| Inner only, random 64-element group chain | 5.103 |
| Sequential outer only | 6.594 |
| Complete historical encoder | **11.085** |
| Historical encoder with the fixed striped outer probe | **12.801** |

The separately measured inner and outer times do not sum exactly to the complete time. The complete path has different cache state and no benchmark boundary between stages.

The 11.085 ms result reproduces the earlier 11.17 ms measurement. The discrepancy is 0.085 ms, or less than one percent.

## Direct comparison

| Construction | Complete time | Relative to historical |
| --- | ---: | ---: |
| Historical group-chain encoder | 11.085 ms | 1.00x |
| Historical fixed striped repair | 12.801 ms | 1.15x |
| Current proof-oriented candidate | 26.001 ms | 2.35x |

The current inner is 2.45 times slower than the historical inner-only kernel. The current complete encoder is 2.03 times slower than the earlier striped repair.

These ratios measure implementations, not isolated mathematical changes. The current candidate changes the inner map, permutation granularity, state transition, outer layout, and graph path simultaneously.

| Feature | Historical implementation | Current candidate |
| --- | --- | --- |
| Randomized inner unit | One contiguous 64-element group | One four-element packet |
| Number of randomized inner units | 32,768 | 524,288 |
| Recursive inner map | Cyclic split BCH transpose | Systematic accumulator transpose |
| Generated inner circuit | 913 field XORs per node | 930 field XORs per node |
| State-coordinate permutation | None | One independent permutation per node |
| Outer data layout | Sequential BCH blocks | Three sloped bands |
| Outer local maps | None in the sequential path | Lane and coordinate permutations |
| Graph work | Included in the historical outer | Separate BCH and dense `R^T` stage |

This table is the starting point for the proof audit. For each added feature, the audit must identify the exact proof obligation that required it. A feature without a necessary proof role should not remain in the deployment candidate.

## Interpretation

The earlier 1.25 ms estimate measured only the change from a contiguous outer layout to a fixed striped outer. It did not include the packet-width-four chain, state permutations, within-band coordinate permutations, data-block assignment, or graph map.

The current construction is therefore more than twice the 11--13 ms frontier that motivated the co-design. The measured costs are distributed across both halves:

- packet-width-four inner plus state permutations: about 12.5 ms;
- randomized outer plus graph map: about 13.5 ms.

The permutations need no standalone full-size pass in this implementation. Their fused cost is still material. The outer random maps add about 3.49 ms over the fixed three-band baseline. The state permutations add about 1.41 ms over the same packet chain.

This probe measures the current proof-oriented candidate. It does not establish that this candidate is the construction we should deploy.

## Reproduction

The source is `libOTe_Tests/RiffleCurrentG4_Bench.cpp` in the `codex/permute-conv-code` libOTe worktree. The generated three-band BCH headers remain experimental build inputs in that worktree's `out` directory.

The measured executable is `out/riffle_current_g4_bench.exe`. Run it with an odd trial count:

```text
out\riffle_current_g4_bench.exe 31
```

The historical source is `out/riffle_profile.cpp` in the same worktree. The rebuild requires `LIBOTE_RIFFLE_PROFILE` because the profiler uses an explicitly gated access hook.

The matched historical commands were:

```text
out\riffle_profile_rebuilt_20260813.exe 1048576 31 full-prefetch2-stride64
out\riffle_profile_rebuilt_20260813.exe 1048576 31 inner-prefetch2-stride64
out\riffle_profile_rebuilt_20260813.exe 1048576 31 outer-sequential
out\riffle_profile_rebuilt_20260813.exe 1048576 31 full-prefetch2-stride64-striped-probe
```
