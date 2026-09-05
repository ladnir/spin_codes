# Full-length designed BA-3 detour

**Status:** inactive and out of scope for the current Riffle outer-replacement
work.

This directory preserves a calculation made after mistakenly treating
designed BA-3 as the complete length-`2^21` code. The calculation studies a
direct sum of constant constituents followed by two permutations and two
accumulators spanning the full codeword.

That is not the active construction. In the active Riffle design, BA is a
candidate replacement for each local `[1024,512]` outer block. Its two
permutations and accumulators are confined to that block. The resulting
length-1024 spectrum must then be composed with the per-block output
permutation, bit transpose, region permutations, and RM2Sub inner.

The detour receipts are retained under `receipts/` only so the calculation is recoverable. They
must not be cited as evidence about the Riffle construction.

The associated exploratory scripts are:

- `scripts/analyze_full_designed_ba_sparse.py`;
- `scripts/analyze_full_designed_ba_bulk.py`;
- `scripts/certify_full_designed_ba_bulk.py`.
