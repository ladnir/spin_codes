# RM2Sub s=16 inner benchmark

`RiffleRm2SubS16_Bench.cpp` implements the complete transposed
SplitState-PreAddMul-RM2Sub inner at `t=128`, `s=16`.  One call processes
2^21 128-bit transposed blocks (32 MiB), or 16,384 epochs.  Input reset and
setup are outside the timed interval.

The benchmark contains independent dense, RM-transform, table, generated
Paar-circuit, fused, and pruned implementations. It checks every optimized
path against the dense reference before timing.

The selected kernel fuses `A` with the `B` workspace fill, uses the exact
minimum-cost 4+4+4+4 partition for the fixed `A` columns, prunes the RM
transpose to the 29 required correlations, and evaluates the eight
quadratic outputs with a 42-XOR circuit. The GF(2^16) update remains a
four-table bitsliced map.

On Peach's Ryzen 9 7950X, pinned to CPU 15, the selected complete inner has
a **2.332785 ms** median over 21 trials. The AVX2-only fallback has a
2.493797 ms median. See
`../receipts/performance/s16_rm2sub_optimized_peach_7950x.json` and
`OPERATION_LEDGER.md` for the complete comparison and exact count.

## Complete encoder integration

The selected inner is also integrated into the frozen exact-permutation
production evaluator. The production transpose runs epochs in reverse, so
the integration applies the transpose of each GF(2^16) multiplication
matrix. It reuses the existing 24-bit route schedules and unchanged BCH
transpose. Each computed `A` output goes directly from its register to its
route slot while one local copy is retained for `B`; there is no materialized
inner word and no extra 32 MiB reread.

On the same Peach core, the complete outer-permutation-inner encoder has a
**10.821225 ms** median over 21 trials, down from **11.333555 ms** for the
s=64 FieldCheckpoint production path. The optimized reverse inner alone is
2.596211 ms. Correctness is checked both against a dense five-epoch inner
oracle and against a separately staged full-size encoder. See
`../receipts/performance/s16_rm2sub_integrated_peach_7950x.json`.
