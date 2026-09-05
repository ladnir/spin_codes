# ExactPerm-BucketRoute evaluator

`ExactPerm-BucketRoute` is a faster evaluator for
`Riffle ExactPerm FieldCheckpoint v1`.  It does not change the code or the
permutation distribution.

The evaluator stores both directions of the exact permutation.  During the
checkpoint accumulation, outputs become final in reverse source order.  The
evaluator immediately appends each final output and its tile-local destination
offset to one of four sequential bucket streams.  It therefore avoids writing
the 32 MiB inner word and reading it again for the permutation.

After routing, the evaluator handles one bucket at a time.  It scatters the
records into an 8 MiB destination tile and evaluates the outer code two blocks
at a time with the AVX2 circuit.  The temporary storage is 32 MiB for values,
8 MiB for 32-bit offsets, and 8 MiB for the destination tile.

The complete output was compared with the frozen staged evaluator before
timing.  The check passed for every tested tile size.  On one fixed Ryzen 9
7950X core, a 21-trial alternating run gave:

- frozen tiled evaluator: 20.841062 ms median in this interleaved run;
- fused bucket evaluator: 13.767435 ms median;
- candidate-to-baseline ratio: 0.660592.

The old evaluator shows a strong cache-order effect in this experiment: its
samples form bands near 19 ms and 21.3 ms.  The bucket samples are mostly near
13.7--13.8 ms.  The directly comparable conclusion is that the exact bucket
evaluator clears the 14 ms target without changing the construction.

The selected tile contains 2048 outer blocks, or 524288 128-bit elements.
Smaller exploratory medians were 15.37 ms for a 1 MiB tile, 15.44 ms for a
2 MiB tile, and 13.79 ms for a 4 MiB tile.  A 16 MiB tile regressed to
18.39 ms.  The 8 MiB tile was the best tested point.

The benchmark source SHA-256 is
`716e47b9963d40908d3bff1f0a2683e1fcdbc261f9e64e9160e8e39818885905`.
The full receipt is `receipts/exactperm_bucket_route_peach_7950x.json`.

