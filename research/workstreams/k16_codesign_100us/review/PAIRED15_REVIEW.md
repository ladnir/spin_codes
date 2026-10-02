# Paired s15 implementation/proof interface review

The reviewed implementation agrees with the new s15 proof's fixed map and
ideal state-update law. The completed fresh receipt
`../proof/paired-s15-drop10-whole-p256.json` also passes independent assembly
checks. Its whole-code bound is below `2^(-62.04117326862128)`, hence below
the requested `2^-40`. No blocking interface mismatch was found. The bound
applies to the explicit ideal independent setup ensemble, not a particular
deterministic seed or implementation timing.

Both `screen_paired_restriction.prepare((10,))` and the C++ literal map delete
zero-based parent row 10, namely `x0*x2 + x1*x4`, from expansion and transpose
feedback. The remaining rows are compacted without further reordering.
Independent evaluation of all 32768 states reproduces the proof spectrum;
the expansion is injective and `CA=0` still holds. All C++ literal feedback
columns match the freshly generated proof columns.

The initial SIMD representation swaps compact state coordinates 3 and 7. Its
matrix preparation conjugates each literal matrix by exactly that
permutation. The generated expansion ignores padded coordinate 15. The
wide-feedback network sets that coordinate to zero. These claims were
checked by evaluating the actual C++ intrinsic expressions with symbolic
binary supports, for both emitted-feedback and raw-feedback variants, not
merely by comparing comments or a second list of declared coefficients.

`customizePaired15` draws fifteen masked 15-bit rows and rejects the entire
matrix if singular. For independent uniform input words, every binary
15-by-15 matrix has the same proposal probability. Conditioning on full
rank therefore gives exactly uniform GL15. Independent proposal words for
successive epochs give independent accepted matrices. The actual `Words`
implementation is a deterministic reproducibility stream; this argument
does not turn a particular seed, or a distribution on 64-bit seeds, into the
ideal independent ensemble.

The installed 16-by-16 matrix is `diag(M,1)`. There are no cross terms
between the active 15-dimensional subspace and coordinate 15, before or
after the 3-to-7 basis swap. From zero initial state, the dummy coordinate
therefore stays zero after every feedback and update. The scalar forward
map uses the adjoint of the literal reverse matrix. Transposition preserves
uniform GL15, so this forward law is the one required by the proof.

The outer field randomizers are unchanged. Their established fixed-nonzero
uniform-image property still applies to the outer symbols. It is not being
used to substitute GF(2^16) scalar updates for the new GL15 state updates.

The unused `rawFeedback=false` compile-time branch in `K16Paired15.cpp`
contains an older full-s16 emission helper. It is not executed by either
public s15 variant. Variant 1 selects raw feedback inside
`paired15::streamHigh<true>`, whose expansion and feedback were checked.
Changing that separate compile-time flag would require a new review.

## Checks and scope

```
C:/Python314/python.exe -B research/workstreams/k16_codesign_100us/review/test_paired15_interface.py -v
```

All five tests pass: full literal census, generated expansion, generated
wide feedback, basis/dummy invariants, and an exhaustive small-dimensional
analogue of the uniform rejection sampler. These portable tests do not
replace the parent's compiled scalar/SIMD, adjoint, alignment, guard, and
in-place tests. No benchmark or remote workload was run in this review.

The completed-receipt checks also pass:

```
C:/Python314/python.exe -B research/workstreams/k16_codesign_100us/review/test_paired15_receipt.py -v
```

These five checks verify all saved source hashes, exact coverage of q=1..512,
the fresh literal map, and the exact positive dyadic endpoint sum followed by
one upward 256-bit rounding. Independent integer arithmetic reproduces the
saved union endpoint:

```
1758346766487946750625787059791224651167255693421756280934152626859851732555 * 2^-312
```

An exact rational comparison confirms it is below `2^-40`. A fresh q=1
calculation with 32 physical t64 steps per region, rather than 16 proof
macros, agrees with the receipt to relative error below `2^-200`. Its outer
shell counts are independently derived by shortened-support
inclusion-exclusion. The displayed q=1 and q=2 margins are respectively
62.2857139563 and 115.0524190862 bits. This review does not repeat the entire
q=2..512 interval replay; it reviews that source path and checks the final
receipt's assembly.

The receipt SHA256 is
`89bed0c4442349a520d54b1aa6e914eac326beb3530123f8273f78fc764b867d`.

Reviewed implementation SHA256 values:

| File | SHA256 |
| --- | --- |
| `K16Paired15.cpp` | `5186ca9e432cec33122350212c7f3566036be7888caeb49c6a55ed911bf504bc` |
| `K16Paired15Setup.cpp` | `f88430e6fdeb21c0ca03690396b53d556d0f2db27791229386692df7ee2b35e4` |
| `K16PairedBasisSetup.cpp` | `8b5808400598086031273895d1d64ab834b8309e3cd9b5a1e6fdbf0cbaadf70c` |
| `T64Paired15Map.h` | `2233d701e7a1d4710755cf2ecba39a3361b87b35d0182f520a806a22df58765e` |

## Final mode52 audit

The winning `K16Paired15Fold.cpp` kernel also passes independent review.
Its fast-to-literal coordinate order is
`(0,1,2,7,3,4,5,6,8,15,10,12,11,13,14,9)`. The dummy coordinate is therefore
at fast slot 9, not slot 15. `customizePaired15Shuffle` first samples the
same literal GL15 matrix as before. It then overwrites the compact tables
with the correctly conjugated matrix, leaving the literal matrix unchanged.
The permutation is not assumed to be its own inverse.

The new symbolic checks execute the actual generated expansion and folded
feedback expressions. Raw-input and state supports are tested simultaneously
for both feedback variants. Expansion agrees with the certified literal
map; feedback state terms cancel by `CA=0`; dummy slot 9 receives zero.
The folded `totals4` expression is exactly the original four-lane reduction.
Conjugated updates leave dummy slot 9 decoupled and preserve its zero value.

The peeled reverse loop emits the final physical block from zero state,
processes middle epochs in descending order through epoch 1, and finally
emits epoch 0. It preserves the original matrix index and state continuity,
and omits only unused terminal state updates. The one- and two-epoch boundary
cases also agree with the unpeeled recurrence. Mode52 uses this kernel,
the reviewed native-field symbol family, and the exact cached
`fieldLoopSharedParity` outer circuit.

These are representation and scheduling changes, not a new construction.
The reviewed **62.04117326862128-bit certificate applies unchanged**. The
four new checks in `test_paired15_fold.py` pass. The complete review suite now
has 28 passing tests; the parent's compiled checks and timing are separate
evidence, not measurements made by this reviewer.

| Final file | SHA256 |
| --- | --- |
| `K16Paired15Fold.cpp` | `4a1437d996b4e3b20668b0c05b9c7e780befd195d7fd8b56b3ff3c359590042c` |
| `T64Paired15ShuffleMap.h` | `076d4b38792ad0f0c4b9533440d86109a17124c49c06a8121086d1b2becb379f` |
| `K16Paired15ShuffleSetup.cpp` | `d6b902d562fc6971d906394e2cc95e41fd3626db5e570f6be5a301c765959b5d` |
| `OuterVariants.cpp` | `f31c1ca780d5cab7aecb6d83e41186990167ebaf0b7fa78112d0b8e9a74e9006` |

For a fresh replay saved outside the default receipt path, set
`SPIN_PAIRED15_RECEIPT` to that JSON file before running the audit. The s16
audit similarly accepts `SPIN_PAIRED16_RECEIPT`. Defaults remain unchanged.

Next, freeze this exact kernel, setup, and reviewed receipt together during
consolidation. Further changes to the state dimension, map, update law, or
routing distribution require a new proof check.
