# K65536 kernel index

The selected experimental pipeline is **mode 52** in [bench.cpp](../bench.cpp).
It evaluates the precomputed transpose on 128-bit blocks, using the paired
t64/s15 inner code and the native-field RS16 outer code. Production sources
and the retained reference kernels are unchanged.

## Selected pipeline

| Stage | Entry point | Implementation |
| --- | --- | --- |
| Inner setup | `customizePaired15Shuffle` | [K16Paired15ShuffleSetup.cpp](K16Paired15ShuffleSetup.cpp), [K16Paired15Setup.cpp](K16Paired15Setup.cpp) |
| Outer setup | `customizeNativeField16` | [K16NativeSetup.cpp](K16NativeSetup.cpp) |
| Inner transpose and routing | `reverseRoutePaired15Fold` | [K16Paired15Fold.cpp](K16Paired15Fold.cpp), [T64Paired15ShuffleMap.h](T64Paired15ShuffleMap.h) |
| Outer transpose | `outervariants::fieldLoopSharedParity` | [OuterVariants.cpp](../outer_variants/OuterVariants.cpp) |

Create an `rs::Plan` with `Variant::Rs16Gf16` and
`InnerKernel::RetainedStreaming`. Prepare `PairedTables`,
`PairedOptimizedTables`, and `NativeOuterTables` through the two setup entries.
Evaluate the inner entry first, then the outer entry. Setup is outside timing.

The inner setup rejection-samples invertible 15-by-15 binary matrices. It
embeds each matrix as `diag(GL15,1)` in `Plan.reverseMatrices` for the scalar
oracle. The fast state uses a fixed coordinate permutation; slot 9 stays zero.
The literal map deletes row 10 from the paired s16 map and compacts the other
rows. The kernel reduces raw-input feedback, peels the initial and final
steps, and uses the six-shuffle `totals4` reduction.

Input contains `plan.n` blocks; output contains `plan.k` blocks. Provide
64-byte-aligned scratch with `plan.scratchBlocks()` blocks. Inner input,
scratch, and setup tables must not overlap. Whole-pipeline output may reuse
the input buffer: the inner stage finishes before outer stores begin.

The exact construction and ideal-setup failure statement are in
[PAIRED_S15_CONSTRUCTION.md](../../../../research/workstreams/k16_codesign_100us/proof/PAIRED_S15_CONSTRUCTION.md).
That statement does not certify a particular deterministic setup seed.

## Correctness checks

The independent inner oracles are `reverseRoutePaired15Scalar` and
`forwardInnerPaired15Scalar` in [K16Paired15Scalar.cpp](K16Paired15Scalar.cpp).
Their literal columns are in [T64Paired15Map.h](T64Paired15Map.h).
[ScalarAdapters.cpp](../ScalarAdapters.cpp) exposes the retained outer
scalar transpose and forward map.

The harness compares complete scalar and SIMD outputs, checks the adjoint
identity, and checks alignment offsets, guards, input preservation, and
whole-pipeline in-place evaluation. A check-only invocation is:

```text
spin_k16_codesign 52 65536 1 0
```

Portable symbolic checks, from the repository root:

```text
python -B spin/experiments/k16_codesign_100us/kernel/verify_paired15.py
python -B spin/experiments/k16_codesign_100us/kernel/verify_paired15_shuffle.py
python -B spin/experiments/k16_codesign_100us/kernel/verify_paired_fold.py
python -B spin/experiments/k16_codesign_100us/outer_variants/test_outer_algebra.py
```

These check the literal columns, feedback rows, coordinate permutation,
zero state slot, feedback-fold identity, and outer arithmetic. The fold check
also authenticates the corresponding function bodies in both frozen sources.
[generate_paired.py](generate_paired.py) carries the immutable paired s16
grouping; it does not require an ignored proposal file.

## Retained comparisons and measurement scope

The s16 counterpart is `reverseRoutePairedFold` in
[K16PairedFold.cpp](K16PairedFold.cpp), using `preparePairedShuffle` and the
oracles in [K16Paired.cpp](K16Paired.cpp). Earlier `Basis`, `Wide`,
`Optimized`, `NoAlias`, and `Shuffle` files preserve isolated comparisons.
`Monomial` and `T128` files are experiments, not the selected construction.

The parent reported 319 passing compiled checks before the selection run.
Mode 52 achieved a 98.9095 us median across eight seed/order medians from
four seeds; those medians ranged from 98.454 to 101.419 us. This meets a
100 us median gate, not a guarantee that every run stays below 100 us.
The final report records holdout results and measurement provenance.

All kernel variants are frozen. Run benchmarks serially, never concurrently.
Next, complete holdout confirmation and preserve the selected sources with
the proof and timing receipts; do not promote rejected variants implicitly.
