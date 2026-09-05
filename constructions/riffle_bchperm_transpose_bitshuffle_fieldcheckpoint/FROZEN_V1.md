# Riffle ExactPerm FieldCheckpoint v1

This record freezes the first optimized implementation of the exact-permutation
construction.  Later changes to the permutation distribution are not part of
this version.

The construction has the following stages.

1. Each 128-bit message block enters `ExtendedBch256x128-Eq3`.
2. Each 256-coordinate outer codeword receives an independent permutation.
3. The resulting 8192-by-256 array is transposed.
4. Each region of 8192 elements receives an independent uniform permutation.
5. `FieldCheckpoint` uses 64 state elements and epochs of 256 elements.

The selected transposed evaluator uses one 8 MiB region-major schedule.  It
processes 512 outer blocks per tile, prefetches 48 source positions ahead, and
evaluates two outer transposes in one AVX2 circuit call.

The frozen benchmark is the 21-trial Peach measurement dated 2026-08-27:

- flat staged baseline: 21.759235 ms;
- selected implementation: 19.363925 ms;
- selected-to-baseline ratio: 0.889918;
- checkpoint inner: 3.096018 ms;
- selected permutation and outer: 15.419695 ms.

The benchmark source has SHA-256
`fe65a974b2db968161d8b52b292d969b34dae6b6a074079b9faae6c0c201d1cf`.
The complete receipt is `receipts/optimized_peach_7950x.json`.

The selected implementation is correct relative to the flat staged evaluator.
The benchmark compares the complete output vectors before timing.

