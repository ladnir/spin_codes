# ExactPerm-Clos256x32

This evaluator factors each uniform 8192-position region permutation into
local permutations of sizes 32 and 256.  It preserves Riffle ExactPerm
FieldCheckpoint v1 exactly.

The factorization removes the structured-permutation proof problem.  It does
not by itself complete the existing 9% certificate, whose outer-spectrum and
exceptional-word obligations remain open.

The C++ evaluator validates the complete 2,097,152-position encoder output
against the existing exact evaluator.  After sweeping outer handoff tiles of
32, 64, 128, 256, and 512 blocks, the selected 256-block route takes 15.94 ms
end to end on Peach, versus 11.57 ms for the promoted bucket route.  Its phase
medians are 3.09 ms for the checkpoint, 5.49 ms for the region-local
factorization, and 7.03 ms for the region-major coordinate handoff plus outer
transpose.  Tile tuning recovered about 0.4 ms, but did not reach the 14 ms
target.

Fusing checkpoint finalization directly into the first 32-way stage removes
the materialized checkpoint word.  The fused path takes 14.61 ms, leaving a
0.61 ms gap to the target.  Its two main phases are 7.69 ms for the fused
checkpoint and local factorization, followed by 6.95 ms for the selected
256-block handoff and outer transpose.

See `CONSTRUCTION.md`, `receipts/route_validation.json`, and
`receipts/clos256x32_peach_7950x.json`.
The tile sweep is in `receipts/clos256x32_tile_sweep_peach_7950x.json`.
The fusion result is in `receipts/clos256x32_checkpoint_fusion_peach_7950x.json`.
