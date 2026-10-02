# Exact-map native outer alternatives

These isolated research entries preserve the native-field outer map and
consume the existing `NativeOuterTables`. There is no setup change, runtime
type erasure, allocation in a hot path, or production dispatch change.
`OuterVariants.cpp`, `OuterShuffleVariants.cpp`, and `OuterNonTemporal.cpp`
are additional translation units; their corresponding headers expose the
entry points under `spin::research::k16codesign::outervariants`.

| Entry | Difference from the retained native field outer |
| --- | --- |
| `fieldInline` | Inline the two parity/output finishes into each group. |
| `fieldUnrolled` | Also fully unroll the fixed sixteen-symbol mix. |
| `fieldSharedParity` | Unrolled mix; share the x0 GF16 basis multiples. |
| `fieldBatch2` | Interleave two groups' unrolled mix; shared parity. |
| `fieldLoopSharedParity` | Shared parity with the original compact symbol loop. |
| `fieldUnaryPack` | Pack with qword unpacks and unary VBMI. |
| `fieldUnaryUnpack` | Output with unary VBMI and qword unpacks. |
| `fieldUnaryBoth` | Both unary-shuffle changes. |
| `fieldLoopSharedNt` | Shared-loop arithmetic; aligned NT output stores and timed SFENCE. |

The shared-parity circuit replaces the six distinct multiples of x0 with
2*x0, 4*x0, and 8*x0 plus four XOR/ternary-logic operations. With ordinary
common-subexpression elimination, this removes three GFNI affine operations
per RS plane, or twelve per group. The unary permutations were checked as
exact byte-source and byte-destination identities. The NT entry uses regular
stores when output is not 64-byte aligned; other buffer/aliasing requirements
are unchanged.

## Validation and measured selection

All translation units compile locally with GCC13.3, `-O3`, the retained
AVX512/VBMI/GFNI flags, and `-mtune=znver4`. Local hardware lacks AVX512, so
compiled execution was left to the parent task. The portable algebra suite
passes all four tests:

```
C:/Python314/python.exe -B spin/experiments/k16_codesign_100us/outer_variants/test_outer_algebra.py -v
```

The parent ran all compiled scalar/SIMD, adjoint, guards, alignment, and
in-place checks before serialized timing. Its reported selection was:

- Keep `fieldLoopSharedParity`: about 51.5us outer versus about 52.3us native.
- Inlining alone was flat. Full unrolling, shared parity with full unrolling,
  and two-group batching were slower. Local generated group text grows from
  about 3.2KiB to 7KiB with unrolling; stack-memory references also increase.
- Unary packing was slower; unary unpacking was approximately flat.
- Reject NT output: reported whole time rose from about 103us to 126–135us,
  including a higher next-call inner cost. Keep cached output stores.

These are the parent's diagnostic observations, not new benchmark runs by
this task or a claim that the <=100us gate has been met. No benchmark or
remote workload was run by the outer-variant author.

Next, retain the cached shared-loop circuit while consolidating the exact
implementation and proof receipt. Avoid repeating the rejected unrolling,
unary-pack, and NT experiments without a new reason to expect a different
result.
