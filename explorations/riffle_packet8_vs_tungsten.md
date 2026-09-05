# Packet-width-8 Riffle compared with Tungsten

The consolidated design record is in `explorations/riffle_packet8_implementation_lessons.md`. This file retains the detailed comparison history and reproduction data.

## Comparison target

This comparison keeps the Riffle inner and outer BCH codes at `[128,64,22]`. Only the permutation atom has width eight.

Both encoders compute a binary linear map on a vector of GF(2^128) elements. The parameters are

```text
K = 2^20
N = 2^21
```

The benchmark uses one pinned Windows thread. MSVC compiles both hot paths with `/O2`, `/arch:AVX2`, and C++20. The two encoders run serially.

## Riffle interface

The Riffle path uses `RifflePacketInnerChain<8>`. Each recursive step retains the 64-element systematic BCH state and gathers eight independently permuted packets of eight elements.

The inner applies the 930-XOR systematic accumulator transpose. It does not use a small BCH code. It also omits the later per-node state permutations.

The outer is the historical sequential `[128,64,22]` BCH transpose, including its graph path. Thus this path measures the intended packet-size experiment rather than the current proof-oriented outer layout.

## Tungsten interface

The Tungsten path uses the production `experimental::TungstenCode` implementation with

```text
ChunkSize = 8
mNumIter = 2
Table = TableTungsten1024x4
```

The coefficient context performs binary addition on GF(2^128) elements. Tungsten reuses a preallocated 16 MiB temporary buffer, matching the Silent OT integration pattern.

Riffle writes a separate `K`-element output, which is also 16 MiB. Allocations and setup are outside both timed loops.

## Measurement

The final run used 101 trials. Times are medians in milliseconds.

| Encoder or stage | Time |
| --- | ---: |
| Riffle packet-width-8 inner | 8.921 |
| Riffle sequential outer | 6.454 |
| Complete packet-width-8 Riffle | **15.428** |
| Complete Tungsten, two iterations | **10.434** |

The complete Riffle path is 1.479 times slower than Tungsten. Tungsten is 4.994 ms faster at this checkpoint.

A separate rebuild of the packet-chain sweep reports 7.672 ms for the packet-width-8 inner. The full comparison reports 8.921 ms because it measures the inner inside the complete encoder harness.

## Exact-map optimization checkpoint

The optimization pass retained the binary matrix, both `[128,64,22]` BCH codes, and the packet-width-eight permutation. It made four implementation changes.

1. The sequential outer circuit now reads each unpunctured BCH word at its input address. The old path copied every 128-element word to a stack buffer.
2. A generated transpose-and-add circuit writes each 64-element output block directly. The old path wrote an intermediate block and read it back to add the shared parity adjoint.
3. The packet chain now prefetches one node ahead. Packet-width-eight work provides enough latency cover without the old two-node distance.
4. The packet-width-eight circuit receives eight explicit packet pointers. This removes repeated indirection through a pointer table inside the generated XOR schedule.

The direct-output outer uses a 1,156-XOR circuit with lower register pressure. The previous 977-XOR transpose-and-add circuit was slower. This result is consistent with the Tungsten comparison: memory traffic and compiler spills matter more than the symbolic XOR count alone.

The optimized circuit was compared with the previous implementation on every output element. The inner and outer comparisons both passed.

The final standalone run used 101 trials:

| Encoder or stage | Time |
| --- | ---: |
| Optimized packet-width-8 inner | 7.765 |
| Optimized sequential outer | 4.017 |
| Complete optimized Riffle | **11.970** |
| Complete Tungsten, two iterations | **9.608** |

The optimized Riffle path is 1.246 times slower than Tungsten. The gap is 2.362 ms.

An interleaved old/new run alternated both Riffle paths on one pinned thread. It used identical inputs and checked equality after every complete map. The old path measured 17.066 ms, and the optimized path measured 13.742 ms. The optimized path therefore used 80.5% of the old time in the cache-stressed A/B harness.

The standalone and interleaved measurements have different working sets. Their absolute times are not directly comparable. The interleaved ratio is the stronger evidence for the optimization gain.

The optimized packet-width-eight result is 0.885 ms above the rebuilt 11.085 ms historical packet-width-64 encoder. Thus reducing the permutation atom from 64 elements to eight no longer causes a large performance loss.

## Locality optimization checkpoint

A second pass tested whether the outer BCH transpose should run as soon as each outer block becomes complete. This fusion preserves the matrix, but it failed as an implementation strategy.

The parity block completed at inner node 1,308. Of the 16,384 complete message blocks, 7,739 completed before parity. The fused path therefore needed a parity pass over 47.2% of the output.

More importantly, completed outer blocks appear in pseudorandom order. Fusing their transposes mixed random outer reads and output writes into the inner state chain. The fused path measured 15.566 ms, versus 11.793 ms for the separate passes. The implementation retains separate inner and outer phases.

Two locality changes did help.

1. The inner now uses write-intent prefetch. Each 128-byte packet is read and overwritten exactly once. In an interleaved kernel comparison, write-intent prefetch measured 7.068 ms, versus 7.415 ms for ordinary L1 prefetch.
2. The outer computes the 24-element graph adjoint before the message blocks. It folds each graph correction into the BCH addend. This removes a later read-modify-write pass over the 16 MiB output.

The graph-fused outer was checked against the previous outer on every output element. An interleaved end-to-end comparison also included write-intent prefetch. The previous optimized path measured 12.063 ms, and the new path measured 11.543 ms. The new implementation used 95.7% of the previous time.

An exact search over 1,024 BCH schedules found a 799-XOR packet circuit at seed 1000. The previous circuit used 803 XORs. An interleaved kernel comparison measured 7.436 ms for seed 1000 and 7.500 ms for seed 18. The production path now uses seed 1000.

The final standalone run used 101 trials:

| Encoder or stage | Time |
| --- | ---: |
| Final packet-width-8 inner | 7.768 ms |
| Final sequential outer | 3.636 ms |
| Complete final Riffle | **11.670 ms** |
| Complete Tungsten, two iterations | **9.767 ms** |

The final Riffle path is 1.195 times slower than Tungsten. It is 0.585 ms above the rebuilt 11.085 ms packet-width-64 checkpoint.

## Abstract memory cost model

The Tungsten comparison suggests a more useful cost model than counting permutations, passes, or element moves. Memory access cost depends on the level that supplies the data, the predictability of the access stream, and the available memory-level parallelism. It is not linear in address distance.

For an access schedule $A$, use the following abstract decomposition:

\[
\operatorname{Cost}(A)
\approx
\operatorname{Stream}(A)
+
\sum_h \frac{m_h(A)\lambda_h}{q_h(A)}
+
\operatorname{TLB}(A)
+
\operatorname{Instructions}(A).
\]

Here, $m_h(A)$ is the number of accesses first served at hierarchy level $h$, $\lambda_h$ is its latency, and $q_h(A)$ is the effective overlap between independent accesses. The stream term accounts for transferred bytes and sustainable bandwidth. This expression is a comparison model, not a cycle-accurate predictor.

The model has threshold behavior. A small change in tile size can move a working set from L1 to L2. A small change in address order can defeat prefetching or increase TLB pressure. Conversely, a permutation can execute cheaply when all candidate elements already reside in L1 or L2.

One PCG element occupies 16 bytes. The relevant locality scales for packet width eight are therefore:

| Scope | Elements | Bytes |
| --- | ---: | ---: |
| One packet | 8 | 128 B |
| 8 packets | 64 | 1 KiB |
| 32 packets | 256 | 4 KiB |
| 128 packets | 1,024 | 16 KiB |
| Complete vector | $2^{21}$ | 32 MiB |

A random gather within a resident 1--16 KiB tile may be inexpensive. The same number of element moves across a 32 MiB vector can incur cache and TLB misses. Sequential access to distant addresses can also beat nearby but unpredictable access because hardware prefetching and bandwidth amortize the distance.

This model changes how proof-motivated neutralizers should be evaluated. The correct question is not merely how many elements a mechanism permutes. The evaluation must identify its locality scope, access order, reuse distance, transferred bytes, and number of outstanding misses.

The main candidate scopes are:

1. A packet-local shuffle randomizes eight resident elements. It can remove fixed slot orientation, but it does not separate elements that already share a packet.
2. A cache-tile shuffle randomizes several adjacent packets. It can remove local co-location and orientation at a cost governed by L1 or L2 traffic.
3. A blocked transpose moves data globally through regular tile transfers. It can change which elements become local without using a packet-at-a-time global gather.
4. A full strided gather randomizes globally, but it exposes the complete vector to irregular cache and TLB behavior.

A promising compromise is a two-level schedule. The encoder shuffles packets and slots inside a cache-resident tile, applies a blocked global remapping, and optionally performs a second local shuffle under the new tiling. The global step supplies long-range movement. The local steps destroy correlations inside each tiling. This construction is only a candidate until its exact distribution, proof effect, and measured cost are established.

The benchmark ladder should compare packet-local shuffles, 1 KiB, 4 KiB, and 16 KiB tile shuffles, blocked transposes, two-level schedules, and the current global packet gather. Each design needs a specialized implementation. A generic permutation routine measures abstraction overhead rather than the attainable construction cost. Each result should report elapsed time and effective bandwidth. Hardware-counter runs should also record cache and TLB behavior when reliable counters are available.

The corresponding proof evaluation should record which relations survive each schedule. At minimum, it should measure outer-block co-location, packet-slot orientation, recursive-state alignment, and the best low-weight attack found under the resulting permutation family.

### Preliminary local-stride calibration

For a tile of $L=8R$ elements, the first calibration uses the deterministic transpose

\[
y_{sR+r}:=x_{8r+s},
\qquad 0\le r<R,\quad 0\le s<8,
\]

inside each aligned tile. This permutation preserves Hamming weight. The standalone kernel materializes every transformed tile in a reusable aligned scratch buffer and copies the result back in place.

Two access orders were optimized and compared. The push kernel streams source packets into eight destination streams. The pull kernel writes one destination stream while gathering resident source elements. A third kernel used 8-by-8 cache-line-oriented microtiles, and a fourth followed precomputed in-place cycles. The latter two were slower at every tested span.

The best 31-trial medians were:

| Tile span | Working set | Best kernel | Separate-pass time |
| ---: | ---: | --- | ---: |
| 64 elements | 1 KiB | pull | 1.796 ms |
| 256 elements | 4 KiB | push | 1.803 ms |
| 1,024 elements | 16 KiB | push | 2.027 ms |
| 4,096 elements | 64 KiB | push | 2.628 ms |

The preferred access order changes between 1 and 4 KiB. This change is direct evidence that one operation count does not describe the locality curve.

The separate pass is not the best implementation for a 1 KiB transform. One Riffle inner node already gathers and prefetches eight packets, totaling 64 field elements. A second implementation transposes this 8-by-8 matrix immediately before the BCH step. It uses 28 explicit swaps, no permutation table, and no additional vector pass.

The corrected harness restores one physical work buffer from a master copy before each timed interval. It alternates the baseline and candidate on one pinned thread. A 31-trial comparison measured:

| Path | Inner | Complete encoder |
| --- | ---: | ---: |
| Frozen packet-width-eight path | 8.537 ms | 12.233 ms |
| Path with fused node transpose | 9.032 ms | 12.701 ms |
| Added cost | **0.495 ms** | **0.468 ms** |

The complete encoder delta is the relevant preliminary cost: about 0.47 ms, rather than the 1.80 ms cost of materializing a full pass.

The fused transform defines a different binary matrix. Its implementation passed an exact 8-by-8 transpose check, but no distance or security claim has been made. The measurement demonstrates the value of design-specific fusion. It does not select the transform for deployment.

### Random local transforms

The next calibration measures random transforms on the 64 elements already gathered by one inner node. The benchmark does not change the two `[128,64,22]` BCH codes or the packet width. Setup generates every finite permutation with rejection sampling, not reduction modulo the group size.

Write the gathered node as an 8-by-8 array $x_{p,s}$. The benchmark implements five randomization families.

1. **Folded packet order.** For each node, setup samples $\pi\gets S_8$. Setup composes $\pi$ into the existing packet-order table. The online kernel is unchanged.
2. **Structured node permutation.** Setup samples independent $\pi,\sigma\gets S_8$. The transform sets $y_{r,j}:=x_{\pi(j),\sigma(r)}$. A 1 KiB scratch tile implements the gather and copy-back.
3. **Full node permutation.** Setup samples one permutation from $S_{64}$. A 1 KiB scratch tile implements the arbitrary local gather.
4. **Independent packet permutations.** Setup samples $\sigma_p\gets S_8$ independently for each packet. The transform sets $y_{p,s}:=x_{p,\sigma_p(s)}$.
5. **Sparse packet mixers.** For each packet, setup samples one or three ordered pairs $(d,s)$ with $d\ne s$. Each pair applies the transvection $x_{p,d}\mathrel{\oplus}=x_{p,s}$. Every transvection is invertible, so every composed packet matrix has full rank.

The mixer benchmark also evaluates a lower-randomness control. The control shares each transvection across all eight packets in one node. The independent variant stores eight times more schedule bytes.

Each invocation processes $2^{21}$ elements of 16 bytes. The benchmark pins one thread to logical core 4. It restores the same 32 MiB work buffer outside every timed interval. It also alternates the baseline and candidate order. The identity control differed by 0.08 ms end to end, which estimates the residual timing noise.

The 31-trial medians were:

| Transform | Online schedule payload | Added inner time | Added complete time |
| --- | ---: | ---: | ---: |
| Folded random packet order | 0 B | -0.209 ms | **0.018 ms** |
| Folded packet order and fixed transpose | 0 B | 0.495 ms | **0.468 ms** |
| Independent packet $S_8$ | 2,097,152 B | 1.151 ms | **1.214 ms** |
| Independent packet $\mathrm{AGL}(3,2)$ | 524,288 B | 2.418 ms | **2.345 ms** |
| Structured $S_8\times S_8$ | 524,288 B | 1.315 ms | **1.297 ms** |
| Full $S_{64}$ | 2,097,152 B | 2.101 ms | **2.209 ms** |
| One shared transvection per node | 32,768 B | 0.350 ms | **0.230 ms** |
| One independent transvection per packet | 262,144 B | 0.381 ms | **0.237 ms** |
| Three shared transvections per node | 98,304 B | 0.705 ms | **0.712 ms** |
| Three independent transvections per packet | 786,432 B | 0.754 ms | **0.709 ms** |

The paired delta is more stable than the absolute time. Baseline complete times ranged from 11.78 to 12.46 ms across these invocations.

Three implementation conclusions follow.

First, packet-order randomness is free online when setup composes it into the existing order table. A separate packet-order pass would measure the wrong design.

Second, independent packet mixers cost almost the same as shared mixers. The independent one-shear mixer adds 0.24 ms. The independent three-shear mixer adds 0.71 ms. Their stronger randomness therefore has negligible extra online cost at these schedule sizes.

Third, regular scratch materialization beats data-dependent in-place swaps. The structured scratch path adds 1.30 ms, while its swap-chain implementation adds 1.64 ms. The full scratch path adds 2.21 ms, while the full swap chain adds 4.26 ms. Scratch traffic is cheaper than unpredictable swap dependencies inside a resident tile.

Fourth, compressed schedules can cost more than their saved bandwidth. The $\mathrm{AGL}(3,2)$ variant stores one 16-bit group index per packet. Its 1,344 permutation rows occupy 10.5 KiB. Eight indexed table reads per node raise the added time to 2.34 ms. Streaming explicit $S_8$ rows is faster despite its 2 MiB schedule.

These results give the following preliminary additive model for this machine:

\[
\Delta T_{\mathrm{complete}}
\approx
0\cdot I_{\mathrm{folded\ packet}}
+0.24\,r_{\mathrm{shear}}
+1.2\,I_{\mathrm{packet\ }S_8}
+0.9\,I_{\mathrm{arbitrary\ cross\mbox{-}packet}},
\]

where times are milliseconds and $r_{\mathrm{shear}}\in\{0,1,3\}$. The final term is the observed gap between independent packet permutations and a full $S_{64}$ permutation. This expression summarizes the measured points. It is not a hardware-independent prediction.

No transform in this section is part of the frozen construction. The proof phase must determine which surviving correlations each family removes. The performance results only define the cost menu for that decision.

## Interpretation

Tungsten retains the stronger performance result. Its chunk width also matches the desired permutation atom.

Riffle performs two expensive classes of work:

1. The recursive inner executes 32,768 full 64-wide BCH steps while gathering eight packets per step.
2. The outer executes 16,384 additional BCH transposes.

Tungsten instead performs two table-driven accumulator passes over the redundant half and one final compression into `K` outputs. Its access pattern and arithmetic are cheaper at this parameter point.

The comparison does not establish equivalent distance or security. Tungsten and Riffle define different binary matrices. The result only compares their production cost for the same input and output dimensions.

The packet-width-8 Riffle path is also not the current certified g=4 construction. It deliberately excludes the state permutations and randomized three-band outer maps whose proof cost is under review.

## Reproduction

The comparison source is `libOTe_Tests/RifflePacket8Tungsten_Bench.cpp` in the libOTe worktree. The command is

```text
out\riffle_packet8_tungsten_bench.exe 101
```

The standalone packet-chain correctness and performance command is

```text
out\riffle_packet_inner_bench_rebuilt_20260813.exe
```

The local-randomization benchmark source is `libOTe_Tests/RifflePacket8RandomLocal_Bench.cpp` in the libOTe worktree. The following PowerShell command pins the benchmark to logical core 4:

```text
$env:RIFFLE_BENCH_CPU=4
out\riffle_packet8_random_local_bench.exe sparse-mixer3 31
```

Replace `sparse-mixer3` with `production-folded-order`, `structured`, `full`, `packet-slots`, `affine-slots`, or `sparse-mixer1` to reproduce another row.
