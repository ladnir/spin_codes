# K20 small outer: 2.962 ms with the unchanged certificate

The selected exact small-outer/s22 implementation now takes **2.962029 ms**
at K=2^20. This is **43.28% less time** than its unchanged K18 kernel and
**10.23% less time** than the larger-outer library baseline in the same campaign.
It retains the full **58.7503767402-bit** setup-failure margin at 10% distance.
No map, routing distribution, randomizer, state dimension, or setup procedure
changed. This report records the research optimization. The subsequent
[library promotion](../../CERTIFIED22.md) retains it as `PacketRsT64S22` and
records a separate matched public-API campaign.

## Matched result

The 2026-10-03 holdout used 128-bit XOR elements, K=1,048,576, and N=2,097,152.
Each process made five warmup calls and 1001 timed in-place transpose calls.
The table reports the median of four process medians, with two processes
per seed. Seeds 113 and 211 were not used to select the optimizations.

| Implementation | Time | Complete margin at 10% distance |
|---|---:|---:|
| Original small outer, cached stores | 5.222449 ms | 58.750377 bits |
| Small outer, streaming stores only | 3.226989 ms | 58.750377 bits |
| **Small outer, packed carry and unrolled outer, streaming stores** | **2.962029 ms** | **58.750377 bits** |
| Larger-outer library baseline | 3.299650 ms | 68.103757 bits |

All timings exclude setup and allocation. All margins concern whole-code
failure under the specified ideal independent setup ensemble. The certificate
and its independent audit are indexed in the
[proof report](../../../research/workstreams/k20_small_outer_screen/README.md).

Both outer choices process the same total number of RS rows and routed records.
At K20, the small outer has 8192 groups of four GF16 RS[16,8] rows;
the larger outer has 4096 groups of eight rows. Both route 32 MiB of data.
The earlier 58% regression was therefore not evidence that smaller groups
intrinsically require more work. The store-policy ablation below isolates
a large implementation cost.

## Selected changes

1. **Streaming route and output stores.** Each route packet writes one aligned
   64-byte cache line. The new route uses non-temporal stores and fences before
   returning to the outer. The outer also uses streaming stores and a final
   fence. Unaligned output falls back to the original cached outer.
2. **Carry state-update outputs between inner steps.** The exact s22 inner
   produces coordinate-interleaved outputs before unpacking them. Retaining
   the two quadratic groups avoids repacking them with cross-vector permutes.
   A byte reversal within each 64-bit lane followed by the original GFNI
   transpose gives exactly the same next-step packed input. Group 0 is unchanged.
   This replaces four cross-vector permutes per step by four local byte shuffles.
3. **Unroll the 16 outer-symbol mixes.** This retains the GF16 randomizer and
   shared-parity circuit. The selected outer uses the same fixed plane finishes
   and streaming output stores. No runtime callback or packet-level policy
   branch was added.

The selected entry points are `reverseRoutePackedNt` in `K20Route.cpp` and
`fieldUnrolledSharedParityNt` in `K20Outer.cpp`. The benchmark mode is
`small-nt-packed-unroll`. The original implementations and the store-only
variants remain available as controls. All policies are selected at compile time
outside the timed loop.

The packed carry does not reinterpret raw GFNI outputs as already packed bits.
Its explicit conversion preserves the original layout, including zero padding
in coordinates 22 and 23. The feedback factor P P^T, independent state updates,
initial state, route order, and final-step convention remain unchanged.

## Ablations and phase diagnosis

The initial seed-41 store screen used 301 calls per process, followed by the
reverse mode order. Median process times were:

| Route stores | Output stores | Time |
|---|---|---:|
| Cached | Cached | 5.219725 ms |
| Streaming | Cached | 3.507418 ms |
| Cached | Streaming | 4.339241 ms |
| Streaming | Streaming | 3.210414 ms |

A separate 501-call schedule screen tested splitting the two outer plane
finishes into noinline functions. That change gained only about 0.4% and was
not selected. Explicit outer unrolling gained about 1.5% in that screen.
The packed-carry route was then compared against its unchanged streaming
control. The selected combination won in both run orders and in the holdout.

The following seed-41 phase measurements use an extra clock between phases.
They are diagnostic, not the headline holdout timings. Phase medians need not
sum exactly to the median of the complete call.

| Implementation | Inner and route | Outer | Complete call |
|---|---:|---:|---:|
| Original cached small outer | 3.357137 ms | 1.864470 ms | 5.222810 ms |
| Streaming stores only | 1.846377 ms | 1.393752 ms | 3.246380 ms |
| Selected packed carry and unrolled outer | 1.613471 ms | 1.343608 ms | 2.958533 ms |
| Larger-outer baseline | 1.692850 ms | 1.585961 ms | 3.279422 ms |

The larger baseline's phase diagnostic calls its existing private stages
directly; its headline measurement still uses `spin::Code`. Compiler inspection
found one ZMM spill/reload in the store-only route and two in the carried route.
Thus the carried version wins despite an extra spill; it is not a spill-free
implementation or an established performance limit.

## Correctness and scope

- The packing conversion passes all 1024 basis bits of its two 512-bit inputs,
  256 deterministic random pairs, and zero-padding checks.
- All eight experiment variants agree with the literal-map scalar route and
  outer, including scratch padding, output guards, four alignments, input
  preservation, and in-place output/tail checks. The scalar forward/transpose
  adjoint check also passes.
- All 12 Release CTests pass. These include K20 seeds 1, 41, and 113 and
  K=2048, 65536, 262144, 264192, 1046528, and 1050624. A separate K20 check
  passes for holdout seed 211.
- All nine store-variant tests also pass with AddressSanitizer and
  UndefinedBehaviorSanitizer. The experiment kernels and scalar oracles were
  instrumented; the unchanged linked production library was not rebuilt with
  sanitizers.
- Every small-outer holdout checksum agrees across implementations for its
  seed. Seed 113 gives `6623330ae9a65586`; seed 211 gives `a4d7de17de341efe`.
- The existing whole-certificate verifier still authenticates all 138 files
  and all 8192 occupancies, and reproduces the exact 58.7503767402-bit union.
  No frozen proof source or receipt was edited.

These tests of neighboring sizes establish implementation correctness, not
new distance certificates. The selected streaming policy was tuned at K20;
do not replace K18's cached policy without a matched small-size comparison.

## Benchmark and reproduction

Hardware and timing conventions match the [initial comparison](README.md).
The host was Peach's Ryzen 9 7950X, CPU 15, 4.5 GHz maximum, boost disabled,
with GCC 15.2, Release optimization, and `znver4` tuning. Input and workspace
both requested large pages; each allocation rounded to 34 MiB. Actual page
backing was not measured. All four shared benchmark locks were held, and
no benchmarks ran concurrently.

The holdout ran these mode orders, each as a fresh process:

- Seed 113: original, large, streaming, selected, selected, streaming, large, original.
- Seed 211: large, selected, original, streaming, streaming, original, selected, large.

Use the build commands in `README.md`, then run the checks before timing:

```sh
ctest --test-dir tmp/k20-small --output-on-failure -j1
taskset -c 15 tmp/k20-small/spin_k20_small_outer small-nt-packed-unroll 1048576 113 1001 1
taskset -c 15 tmp/k20-small/spin_k20_small_outer large 1048576 113 1001 1
```

For phase diagnostics, append `1`. Repeat serially in reversed order with
the benchmark locks held. `small` selects the frozen cached control;
`small-nt` selects only the two streaming-store changes.

The retained base source archive is
`tmp/k20-small-benchmark-source-20261003-v2.tar.gz`, SHA256
`2d6e1ee1e9248bfdd9254d941b906aff79cec52882331fbdb7ffcf676a138b27`.
Overlay `tmp/k20-small-optimization-v4.tar.gz`, SHA256
`de85833dc8a909f224585a50530689f21521eead13754acbfec00d1d9a9a832b`,
to recover the measured experiment. The current report was written afterward.
The measured `K20Route.cpp`, `K20Outer.cpp`, and `bench.cpp` respectively have
SHA256 hashes:

```text
4e83c414611aafaf3725876cd7a0dbf828b7e061c5287bb50ab26c3454935828
33c09b07ad6a48465427947c6247362cc73bcdf0097e748212d82d239b96d8c9
98660a422e66abbfd72fa8b66138c8ffa7aa1b8e166983cffcd8d56ca316e2c5
```

Raw process summaries remain ignored:

| Local path under `tmp/` | SHA256 |
|---|---|
| `k20-small-opt-stores-screen-v1.csv` | `0dedae064b744dabf98d64b29d5cc7551d858303113e3b96f6737790f7a6c40e` |
| `k20-small-opt-stores-phases-v1.csv` | `80d2ec1487d70c5530e2f2937212597a8c4864b6760b5c02d8eafa9795ebd4d5` |
| `k20-small-opt-outer-screen-v2.csv` | `b4fea738c2399988d447544c61eb78dcdb619bcad3503afccd19fe54dbdbbd07` |
| `k20-small-opt-carry-screen-v3.csv` | `ed97188c892dc13b8e39b72f493e4f3efa7781df3c8431f38734c65f35ae77a1` |
| `k20-small-opt-holdout-v4.csv` | `2684b7d59594f07cdd8b31c5a151d6aed8408602e767c3c00959db1df75fbbf5` |

The remote build and logs are in `/tmp/spin-k20-small-opt-6I28Vl`.
The final holdout executable SHA256 is
`d12c0fc14402597b0be5c42950964aaa411f0e80ede0503c515ea367d7237cc5`.

The subsequent [promotion record](../../CERTIFIED22.md) completes the cached
K18 comparison and introduces the exact s22 library profile. K18 selects
packed carry with the cached loop outer; K20 selects packed carry and the
streaming unrolled outer. The larger-outer profile remains available.
