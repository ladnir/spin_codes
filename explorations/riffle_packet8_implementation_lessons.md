# Packet-width-eight Riffle: construction and implementation lessons

## Purpose

This note records the lessons from the packet-width-eight implementation pass. Its purpose is to prevent two recurring errors.

First, a proof-motivated change must not enter a frozen construction without review. Second, an implementation optimization must preserve the exact binary matrix unless the construction change is explicit.

The discussion concerns the transposed map used by the PCG application. The input is a vector

\[
c\in (\mathbb F_2^{128})^N,
\qquad
N=2^{21}.
\]

The encoder computes a vector of $K=2^{20}$ field elements. One field element is one 128-bit value. The hot path therefore measures field-element XORs and memory traffic, not packed-bit operations or popcounts.

## Frozen benchmarked map

The benchmarked packet-width-eight map has the following fixed properties.

1. The inner code is the extended BCH code `[128,64,22]`.
2. The outer code is also `[128,64,22]`.
3. The recursive inner state contains 64 field elements.
4. The permutation atom contains eight field elements, or 128 bytes.
5. The packet permutation is one bijection on $N/8$ packets.
6. The outer stage uses 16,386 BCH blocks and 256 punctured coordinates.
7. The graph subcode has codimension 24.

The packet width changes the address schedule. It does not replace either BCH code with a smaller code.

The timed path does not add a packet-internal shuffle. Any such shuffle would require a separate construction and benchmark.

The performance comparison excludes later proof experiments that added full-size permutations, state permutations, or fresh permutations inside packets. Those changes define different matrices and require separate review.

## Efficient transposed encoder

The implementation consumes a mutable input vector `word[0..N-1]` and writes a disjoint output vector `output[0..K-1]`.

The setup phase expands a short seed into the packet bijection and the graph-subcode columns. Setup is outside the timed encoder.

The online encoder has two phases.

### Inner packet chain

Let `order` be the packet bijection. Each logical node owns eight packet indices. The encoder processes logical nodes in descending order.

For each node, the encoder performs these operations:

1. Prefetch the eight packets for the next node with write intent.
2. Obtain eight pointers to the current 128-byte packets.
3. Apply the transposed systematic BCH circuit to those packets and the 64-element state.
4. Write the updated packet elements in place.
5. Replace the state with the 64 transposed accumulator outputs.

The state dependency fixes the node order. The eight physical packets remain pseudorandom addresses.

### Sequential outer transpose

The outer stage interprets the modified input in sequential outer-code order.

1. Transpose the final BCH block to obtain the shared parity adjoint.
2. Transpose the final data block to obtain the 24 graph-adjoint elements.
3. Build the three graph lookup tables from the graph adjoint.
4. Process the $K/64$ message blocks in increasing address order.
5. For each block, combine its parity and graph corrections into one 64-element addend.
6. Run the generated BCH transpose-and-add circuit directly into the output.

The first 256 data blocks omit coordinate zero. The implementation gathers only those 127-element words into a local 128-element buffer. Every later BCH word is read directly from its input address.

This order computes the same linear map as the earlier outer stage. It changes only when the graph correction is applied.

## Cost-model lessons

### Random access dominates simple operation counts

The packet-width-eight inner performs eight pseudorandom packet accesses per node. The BCH state dependency prevents node reordering. A reduction in symbolic XOR count is useful only when it does not increase spills or memory stalls.

The same principle explains why a 1,156-XOR outer circuit beat a 977-XOR alternative. The larger circuit generated better machine code because it used fewer costly spills.

### A pass is not expensive merely because it is a pass

Tungsten performs repeated passes, yet it remains faster. Its passes are regular and stream through memory.

The relevant questions are:

- Are addresses sequential or pseudorandom?
- Does the pass reread cold data?
- Does it require an additional output read-modify-write?
- Can hardware prefetching cover its latency?
- Does the working set evict data needed by the next phase?

Pass count alone did not predict performance.

### The PCG element type rules out packed-bit intuitions

Each logical symbol is a $\mathbb F_2^{128}$ element. An XOR operates on the complete element. Bit packing and popcount do not implement the required transposed map.

An optimization should therefore reduce block XORs, address-generation work, cache misses, or complete memory passes.

### Exact specialization is appropriate in the hot path

The packet-width-eight circuit now receives eight explicit pointers. The hot path does not use a pointer table, dynamic dispatch, or a heap-allocated callback.

This specialization improved performance without changing the general packet-chain interface. The implementation retains generic code outside the packet-width-eight kernel.

## Accepted optimizations

Every accepted optimization passed an element-by-element comparison against the preceding exact implementation.

### Direct outer reads and writes

Unpunctured outer blocks are read at their final addresses. A generated transpose-and-add circuit writes final output blocks directly.

These changes removed a 128-element stack copy and a later parity-add pass for each block.

### Write-intent packet prefetch

Each packet is read and overwritten once. Write-intent prefetch requests cache-line ownership before the hot step.

The isolated inner comparison measured 7.068 ms with write intent and 7.415 ms with ordinary L1 prefetch.

### Graph correction inside the BCH addend

The outer stage computes the graph adjoint before the message blocks. It combines each graph correction with the parity adjoint.

This optimization removes a later read-modify-write pass over the 16 MiB output.

The fused graph path applies when $K$ is divisible by 64. The generic outer path handles a final block that mixes message and graph coordinates.

### Smaller verified BCH schedule

An exact search over 1,024 Paar schedules found seed 1000. The new schedule uses 799 transpose XORs. The preceding schedule used 803.

The interleaved inner comparison measured 7.436 ms for seed 1000 and 7.500 ms for seed 18.

## Rejected directions

### Smaller BCH codes

A `[16,8,4]` experiment was faster, but it changes the code and weakens the intended impulse structure. The deployment path retains both `[128,64,22]` codes.

### Immediate inner and outer fusion

An outer block becomes available after the inner chain has written every packet that intersects it. The parity block completed at node 1,308 for the benchmark seed.

Of the 16,384 complete message blocks, 7,739 completed before parity. A fused implementation therefore needed a later parity pass over 47.2% of the output.

The larger problem was address order. Outer blocks completed in pseudorandom order. Their transposes introduced random reads and output writes into the inner state chain.

The fused path measured 15.566 ms. The separate phases measured 11.793 ms in the same comparison. The implementation retains the separate phases.

### Construction changes hidden as implementation details

Extra full-size permutations, state permutations, lane randomizers, or packet-internal shuffles are not implementation optimizations. Each mechanism changes the binary matrix or its distribution.

Such a mechanism requires a construction delta, a cost estimate, an implementation, and a benchmark before it enters a frozen specification.

## Performance checkpoints

Standalone runs use one pinned Windows thread. MSVC uses C++20, `/O2`, and `/arch:AVX2`. Times are medians over 101 trials.

| Checkpoint | Riffle | Tungsten | Ratio |
| --- | ---: | ---: | ---: |
| Initial matched packet-width-eight path | 15.428 ms | 10.434 ms | 1.479 |
| First exact optimization pass | 11.970 ms | 9.608 ms | 1.246 |
| Final locality pass | 11.670 ms | 9.767 ms | 1.195 |

The final Riffle stages measured 7.768 ms for the inner and 3.636 ms for the outer.

Standalone medians vary with machine state. Interleaved A/B runs provide stronger evidence for local changes.

The final interleaved comparison measured 12.063 ms for the preceding optimized path and 11.543 ms for the final path. The final path used 95.7% of the preceding time.

The rebuilt historical packet-width-64 encoder measured 11.085 ms. The final packet-width-eight path is 0.585 ms slower at the same dimensions.

## Construction and proof discipline

Future work should apply the following sequence.

1. State the exact construction delta.
2. Identify which binary matrix entries or randomness distributions change.
3. Explain why the proof needs the change.
4. Estimate its online cost in the PCG model.
5. Implement the candidate without replacing the frozen path.
6. Check exactness when the candidate claims to preserve the matrix.
7. Benchmark the candidate against the frozen path.
8. Freeze the new construction only after reviewing both proof and performance evidence.

A proof certificate must bind the construction it proves. A proof for a global-lane sampler, added permutation, or different puncture law does not certify the earlier encoder automatically.

## Implementation anchors

The optimized implementation is in the libOTe worktree:

- `libOTe/Tools/RiffleCode/RifflePacketInnerChain.h` contains the packet-width-eight chain and write-intent prefetch.
- `libOTe/Tools/RiffleCode/ExtendedBch128x64SystematicPacket8TransposeCircuit.h` contains the generated seed-1000 circuit.
- `libOTe/Tools/RiffleCode/RiffleCode.h` contains direct outer reads, direct output writes, and graph fusion.
- `libOTe/Tools/RiffleCode/RandomGraphSubcode.h` contains the prepared graph-table interface.
- `libOTe_Tests/RifflePacket8Tungsten_Bench.cpp` contains the matched encoder comparison.

The exact circuit generator is `scripts/probe_systematic_bch_circuit.py` in this repository.

## Next measurement gate

The next benchmark should use the actual PCG/libOTe call site. It should retain the application's buffer ownership, allocation reuse, and surrounding cache state.

That benchmark should answer four questions:

1. Does the caller already provide disjoint $N$-element input and $K$-element output buffers?
2. Can setup tables persist across PCG invocations?
3. Does the caller consume the output immediately enough to benefit from the graph-fused write order?
4. Does a release build preserve the interleaved kernel gains?

No further construction change is justified by the current performance data alone.
