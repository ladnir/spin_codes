# Proof status

## Inherited regular theorem

The parent's regular theorem uses only these properties of `A`:

1. `A` has minimum distance at least 20;
2. every output coordinate is a nonzero linear form; and
3. every two or three distinct coordinate forms are linearly independent.

For the AuxW3 map, complete enumeration of all \(2^{32}\) image words gives
minimum distance 24. The exact coordinate audit finds 128 nonzero, distinct
forms and no dependent triple. Its first three factorial moments match the
required identities exactly. A direct matrix audit proves `BA=0`.

The map therefore meets every inner assumption of the parent's fixed-inner
regular theorem. No numerical recurrence needs to be weakened. At 9%
relative distance, the inherited aggregate regular margin is 80.8139 bits.

## Scope

The result remains conditional on the modeled outer spectrum and nearest
binary64 arithmetic. Mixed occupations containing at least two all-one outer
words and more than 64 total active blocks remain open.

## Evidence

- `../riffle_bchperm_transpose_bitshuffle_splitstate_preaddmul_t128_s32/receipts/lighter_a32/sweep_summary.json`
- `../riffle_bchperm_transpose_bitshuffle_splitstate_preaddmul_t128_s32/receipts/lighter_a32/a32_d2_w3_selection.json`
- `../riffle_bchperm_transpose_bitshuffle_splitstate_preaddmul_t128_s32/receipts/lighter_a32/a32_d2_w3_generators.txt`
- `../riffle_bchperm_transpose_bitshuffle_splitstate_preaddmul_t128_s32/receipts/lighter_a32/a32_d2_w3_histogram.csv`
- `../riffle_bchperm_transpose_bitshuffle_splitstate_preaddmul_t128_s32/receipts/lighter_a32/a32_d2_w3_spectrum_audit.json`
