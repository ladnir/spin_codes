# Cache-friendly evaluation of the Riffle permutation

## Objective

`Riffle ExactPerm FieldCheckpoint v1` takes 19.36 ms on Peach.  The checkpoint
inner takes 3.10 ms.  The permutation and outer stage takes 15.42 ms.  A 14 ms
end-to-end target leaves about 10.9 ms for the permutation, outer, and boundary
overhead.

The first objective is to reduce memory latency without changing the
permutation distribution.  This order matters because the existing distance
analysis assumes independent uniform permutations.

## ExactPerm-BucketRoute

Let (P:[N]\to[N]) be the current map from an outer coordinate to its inner
coordinate.  The transposed evaluator computes

\[
    y_q := x_{P(q)}.
\]

Let (Q:=P^{-1}).  Fix a destination tile size (T) that is a multiple of 256.
Write (M:=N/T).  The proposed evaluator uses (M) buckets.

The first pass scans the inner array in increasing source order.  For each
source position (p), it computes

\[
    q:=Q(p),\qquad b:=\lfloor q/T\rfloor,\qquad \ell:=q\bmod T.
\]

The evaluator appends the record ((\ell,x_p)) to bucket (b).  Each append is
sequential within its bucket.  Since (Q) is a bijection, every bucket receives
exactly (T) records.

The second pass processes one bucket at a time.  It reads the bucket records
sequentially and assigns the value in record ((\ell,v)) to tile position
(\ell).  It then applies the packed outer transpose to that tile.

For 512 outer blocks per tile,

\[
    T=512\cdot256=131072
\]

elements.  One destination tile occupies 2 MiB, and (M=16).  The first pass
therefore maintains 16 sequential value streams and 16 sequential offset
streams.  The second pass scatters only within a 2 MiB tile.

This evaluator requires the following temporary storage:

- 32 MiB for bucketed 128-bit values;
- 8 MiB for 32-bit local offsets;
- 2 MiB for one destination tile.

The inverse permutation schedule occupies 8 MiB.  The value buffer can reuse
the 32 MiB workspace of the flat staged evaluator.

## Why the route may be faster

The current tiled evaluator reads random 16-byte values from each 128 KiB
region.  Different destination tiles often consume different values from the
same 64-byte source cache line.  The evaluator reloads that line for each tile.

`ExactPerm-BucketRoute` reads every source cache line once in the first pass.
It converts global random access into sequential input reads and a small number
of sequential output streams.  The second pass confines random writes to one
2 MiB tile.

The route moves more bytes than the current evaluator.  Its advantage depends
on replacing random cache misses with sequential bandwidth.  The benchmark
must therefore include both routing passes and the packed outer circuit.

## Proof status

`ExactPerm-BucketRoute` does not define a new code.  It computes the same map
(P) as `Riffle ExactPerm FieldCheckpoint v1`.  The distance analysis and its
probability space remain unchanged.

The implementation must verify the complete output against the frozen staged
oracle.  This check authenticates evaluation equivalence.  It does not replace
the existing distance argument.

## Construction changes to consider later

If exact bucket routing remains too slow, a later construction can restrict
the permutation family.  Each restriction needs a distinct name and a new
distance argument.

`Riffle TilePerm FieldCheckpoint` would partition each region into contiguous
tiles.  It would permute the tiles and independently permute positions within
each tile.  This family preserves the random tile partition but not the uniform
law on (S_{8192}).  A proof must control active outer blocks that cluster in a
source tile.

`Riffle WindowPerm FieldCheckpoint` would permute only within fixed contiguous
windows.  This family is cheaper, but it exposes the original block positions.
It is the least promising proof candidate because a clustered outer support can
remain clustered after permutation.

A two-pass row-and-column shuffle is an intermediate option.  It can move data
through cache-sized rows and columns.  Its proof would require a contingency
table analysis for the row and column occupancies.

## Next experiment

Implement `ExactPerm-BucketRoute` before changing the construction.  Sweep
destination tiles of 256, 512, and 1024 outer blocks.  Measure the route alone
and the complete encoder on Peach.  Retain the frozen staged evaluator as the
oracle and alternate baseline and candidate trials.

## Result

The exact route works and preserves the complete output.  Fusing the reverse
checkpoint pass with source-major bucket routing was important: each finalized
inner output goes directly to a sequential bucket stream, so the evaluator
does not materialize and reread the 32 MiB inner word.

The best tested tile contains 2048 outer blocks and occupies 8 MiB.  A
21-trial alternating Peach run measured 13.767435 ms for the fused bucket
evaluator.  This clears the 14 ms target.  A 4096-block, 16 MiB tile regressed
sharply, so the useful point lies between too many concurrent bucket streams
and too large a random-scatter working set.

The implementation and exact samples are recorded in the construction folder
under `EXACTPERM_BUCKET_ROUTE.md` and
`receipts/exactperm_bucket_route_peach_7950x.json`.
