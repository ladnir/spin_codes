# Riffle BCHPerm-TransposeBitShuffle-StructuredStepConv

This folder records the first coherent structured design after the random
outer and random inner experiments. The family has two inner variants:

- FieldMulStepConv\((t,s)\);
- ToeplitzStepConv\((t,s)\).

Both variants use the same outer code, coordinate permutations, transpose,
and regional bit permutations. They also have the same one-vector inner
transition law. One distance calculation therefore covers both variants.

The current outer spectrum is a conjectured BCH-like
\([256,128,38]\)-shaped model. The saved distance values cover one active
outer block. They are not complete distance certificates.

Files:

- `CONSTRUCTION.md`: exact forward construction and transposed evaluator;
- `RESULTS.md`: parameter sweep and optimized kernel measurements;
- `receipts/one_active_t_s_sweep.json`: numerical distance diagnostics;
- `receipts/algebraic_family_checks.json`: algebraic checks;
- `receipts/bench/`: serial bitsliced benchmark receipts;
- `receipts/direct_kernel_comparison.json`: serial direct-kernel results and
  comparison with the bitsliced baseline;
- `receipts/block_toeplitz_algorithm_comparison.json`: Four Russians,
  diagonal AVX2, and Karatsuba comparisons with block-resident values;
- `receipts/four_russians_cost_profile.json`: movement floor, row-expansion
  tradeoff, and sampled phase costs.
