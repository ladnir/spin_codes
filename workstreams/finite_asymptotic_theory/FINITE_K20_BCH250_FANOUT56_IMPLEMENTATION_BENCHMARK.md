# BCH250 row-local Fanout-56 implementation benchmark

## Outcome

The packed Fanout-56 kernel adds 2.171294 ms to one post-BCH encoding at the
finite \(k=2^{20}\) parameters. The zero-layer control takes 0.006059 ms, and
the 56-layer path takes 2.177353 ms. These medians come from the same binary,
input, thread, and timing harness.

This result measures the row-local fanout wrapper only. It is not a timing for
the complete proof-model encoder.

| mode | median | p10 | p90 |
|---|---:|---:|---:|
| zero layers | 0.006059 ms | 0.005716 ms | 0.007522 ms |
| 56 layers | 2.177353 ms | 2.067613 ms | 2.280878 ms |
| isolated fanout increment | 2.171294 ms | -- | -- |

The benchmark ran on one pinned logical processor of a 13th Gen Intel
Core i7-13700H. MSVC 19.50 compiled the source with C++20, LTCG, AVX2, and
the optimization flags `/O2 /Oi /Ot /GL`.

## Implemented map

Each encoding contains 8,576 rows of 250 bits. The implementation stores one
row in four 64-bit words. Each row has an independent table of 56 layers.

A layer contains disjoint source and target sets \(S,T\subseteq[250]\), where
\(|S|=31\) and \(|T|=33\). For an input row \(x\), the layer computes

\[
 b:=\bigoplus_{i\in S}x_i
 \qquad\text{and}\qquad
 x_j\mathrel{\gets}x_j\oplus b\quad(j\in T).
\]

The hot path stores \(S\) and \(T\) as eight 64-bit masks. Four masked loads
and one population count compute \(b\). A branch-free all-zero or all-one mask
then updates the four words of \(x\).

The layer is an involution because \(S\cap T=\varnothing\). The inverse of a
56-layer composition therefore applies the same layers in reverse order.

## Correctness checks

The executable performs these checks before timing:

1. It validates every source and target weight and verifies disjointness.
2. It compares the packed kernel with a scalar bit implementation on 32 rows.
3. It applies all layers and their reverse order to all 8,576 rows.
4. It checks that the reverse composition recovers the complete input.

The receipt records distinct deterministic checksums for the input and
Fanout-56 output. The zero-layer checksum equals the input checksum.

## Setup, memory, and throughput

There are \(8576\cdot56=480{,}256\) row-layer applications. The measured
increment corresponds to 221.184 million applications per second.

The schedule stores 64 bytes per layer, or 30,736,384 bytes in total
(29.3125 MiB). The kernel streams this table at an effective 14.156 GB/s.
The physical row buffer occupies 0.261719 MiB. Schedule construction took
257.613 ms and is outside the online timing.

The schedule table dominates the wrapper's memory traffic. The measurement
therefore favors the current explicit-mask representation for online speed.
A compressed representation could reduce memory, but it would add decoding
or regeneration work to the hot path.

## Probability-space boundary

The proof samples each row's 56-layer composition independently from the
uniform layer distribution. The benchmark uses SplitMix64 with a fixed seed
to generate one reproducible schedule table. SplitMix64 is only a benchmark
fixture. It does not instantiate the proof's setup distribution or a
cryptographic setup algorithm.

The measured kernel implements the layer action for any valid schedule table.
A proof-compatible implementation must either sample the table with an
approved uniform sampler or authenticate a table produced by such a sampler.

## Scope of the comparison

Both timed modes begin with the same packed post-BCH rows and copy them into
the working buffer. The zero-layer mode stops after the copy. The 56-layer
mode applies the complete row-local mask table.

Both modes omit:

- the encoder for the fixed \([250,124,\ge38]\) BCH subcode;
- the independent local coordinate permutations;
- the transpose and region permutations; and
- the RM2Sub-S19 inner map.

Consequently, the 359-fold ratio between the two modes is not a useful
full-encoder ratio. The isolated 2.171294-ms increment is the relevant result.

## Remaining implementation obligations

1. Bind the codimension-one BCH basis coordinate and implement the exact
   \([250,124,\ge38]\) encoder.
2. Implement a proof-compatible setup sampler or bind an authenticated
   precomputed schedule table.
3. Implement or fuse the local coordinate permutations.
4. Integrate the transpose, region permutations, and RM2Sub-S19 map.
5. Prove equivalence between the proof-model map and each optimized circuit.
6. Benchmark the complete zero-layer and 56-layer encoders in both directions.

The implementation source is `benchmark_bch250_rowlocal_fanout56.cpp`. The
machine-readable receipt is
`bch250_rowlocal_fanout56_benchmark_windows_i7_13700h.json`.
