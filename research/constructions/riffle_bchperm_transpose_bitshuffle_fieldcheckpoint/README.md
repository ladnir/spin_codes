# Riffle Transpose-BitShuffle FieldCheckpoint

Frozen implementation name: **Riffle ExactPerm FieldCheckpoint v1**.

This name identifies the exact permutation family and the selected transposed
evaluator.  A construction that changes the permutation distribution requires
a distinct name.  A new evaluator for the same exact permutation may use a
distinct implementation suffix.

This implementation baseline instantiates the first complete transposed
encoder for the spectrum-permuted bit-transpose construction:

1. `ExtendedBch256x128-Eq3`, an explicit binary `[256,128,>=38]` code;
2. an independent permutation of the 256 coordinates of every outer block;
3. the 8192-by-256 bit transpose;
4. an independent permutation inside every transposed region of length 8192;
5. a 64-lane checkpoint accumulator with 256-element epochs and an independent
   nonzero `GF(2^64)` state multiplier at each epoch boundary.

The transposed implementation composes steps 2--4 into one 8 MiB lookup
schedule.  The selected full path applies the checkpoint inner in place, then
evaluates the permutation region-by-region in tiles of 512 outer blocks.  The
tile fits in 2 MiB.  It prefetches 48 region entries ahead and evaluates two
outer transposes at once with an AVX2-packed version of the generated circuit.
A flat staged path is retained as a correctness oracle and baseline.

The explicit outer is the parity extension of primitive narrow-sense
`BCH(255,37)`, giving `[256,131,>=38]`, intersected with the three checks
`x_0+x_1=x_0+x_2=x_0+x_3=0`.  The resulting rank is 128.  Its generated
transposed circuit uses 2961 XORs per outer block.

The initial end-to-end measurement is in
`receipts/end_to_end_baseline_i7_13700h.json`.  It shows that the checkpoint
inner is already below the earlier 7 ms phase target.  The dominant cost is
now the globally scattered structured permutation, not the field checkpoint.

The optimized measurement is in `receipts/optimized_peach_7950x.json`.  On a
fixed Ryzen 9 7950X core, 21 alternating trials measured 21.759 ms for the flat
staged baseline and 19.364 ms for the selected tiled and packed path.  This is
an 11.0% end-to-end reduction.  Every selected output is checked against the
flat staged oracle before timing.

The frozen implementation record is in `FROZEN_V1.md`.  The next implementation
exploration is `ExactPerm-BucketRoute`.  It preserves this construction and
changes only how the transposed permutation is evaluated.

That exploration succeeded.  The fused `ExactPerm-BucketRoute` evaluator uses
an 8 MiB tile and has a 13.767 ms median on Peach over 21 alternating trials.
It clears the 14 ms target while producing exactly the same output as the
frozen evaluator.  Its design and measurement are recorded in
`EXACTPERM_BUCKET_ROUTE.md`.

The evaluator was subsequently promoted and compressed into two immutable
24-bit schedules.  The production implementation uses 12 MiB of setup
schedule, a 40 MiB reusable workspace, and no hot-path allocations.  Its
21-trial Peach median is 11.334 ms with complete-output correctness.  See
`PROMOTED_BUCKET_ROUTE24.md` and
`receipts/promoted_bucket_route24_peach_7950x.json`.
