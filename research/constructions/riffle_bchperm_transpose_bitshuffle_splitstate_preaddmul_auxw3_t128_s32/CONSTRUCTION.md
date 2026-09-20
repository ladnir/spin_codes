# Construction

This construction is identical to
`riffle_bchperm_transpose_bitshuffle_splitstate_preaddmul_t128_s32`, except
for the fixed map

\[
A:\mathbb F_2^{32}\longrightarrow\mathbb F_2^{128}.
\]

The first 32 output coordinates are systematic. Each of the 32 information
coordinates enters three of the next 64 coordinates. Two fixed accumulator
layers, separated by a fixed permutation, mix those 64 coordinates. The
last 32 coordinates enforce `BA=0` for the parent's fixed map `B`.

The selected seed is 2736803764. The exact generator columns are recorded in
the parent's `receipts/lighter_a32/a32_d2_w3_generators.txt` file.
