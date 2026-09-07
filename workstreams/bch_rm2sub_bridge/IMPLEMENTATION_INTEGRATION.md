# Imported implementation: validation and grid correspondence

The imported code is useful implementation infrastructure, but it is not yet an
implementation of the selected BCH grid cell. This note records the first local
build and the exact constituent comparison on 2026-09-07.

The source snapshot is commit
`3941f2e125fa0b7c9fe8c02694d12f306f82bbed` from worktree 5371.
Both imported directories remain byte-identical to that commit. Build outputs
and the numerical comparison receipt are under `out/implementation_validation/`.
No benchmark was run and no historical performance receipt was changed.

## Local build and correctness

Both standalone projects build with MSVC 19.50.35728, Ninja, C++20, Release,
and the projects' source-local `/arch:AVX2` flags. The host is an Intel i7-13700H.

| Test | Result |
|---|---|
| Golay reused and independent modes | PASS: scalar reference and both historical checksums |
| Golay public header | PASS |
| Structured SPIN public header | PASS |
| Structured SPIN complete correctness executable | FAIL: historical checksum mismatch |

The Golay checksums are `c4108a1282911d32` and `2bd00619556f9a86`.
The Structured SPIN executable reports successful BCH and RM2Sub component
checks. Its staged-versus-fused comparison also completes before the checksum
exception. Its subsequent checked-interface comparison is not reached.

The historical Structured SPIN aggregate is `10793394525839794159`; this build
produces `4573947736370345000`. The setup calls `std::shuffle` and
`std::uniform_int_distribution`. A seed does not specify their cross-library
output sequences. This is a plausible explanation for the GCC/MSVC difference,
not a confirmed diagnosis. The historical checksum remains enforced.

The first Golay build failed to find dependency headers through the long,
unnormalized default include path. Passing the resolved dependency root fixes
the build without changing imported files. From the repository root:

```powershell
$structuredRoot = (Resolve-Path 'constructions/riffle_parityfanout31x33_bchperm_transpose_bitshuffle_splitstate_preaddmul_rm2sub_t128_s19/implementation').Path
cmake -S workstreams/implementation_cleanup/golay_ba3_b240_rm2sub_s19 -B out/implementation_validation/golay -G Ninja -DCMAKE_BUILD_TYPE=Release "-DSTRUCTURED_SPIN_IMPLEMENTATION_ROOT=$structuredRoot"
cmake --build out/implementation_validation/golay --parallel 4
ctest --test-dir out/implementation_validation/golay --output-on-failure
cmake -S "$structuredRoot" -B out/implementation_validation/structured -G Ninja -DCMAKE_BUILD_TYPE=Release
cmake --build out/implementation_validation/structured --parallel 4
ctest --test-dir out/implementation_validation/structured --output-on-failure
```

## The S19 constituent matches spectrally, not as an ordered matrix

The grid fixes a selected map with step width 128 and state dimension 19.
The imported kernel fixes a different map with those dimensions.
`audit_imported_rm2sub.py` authenticates the grid snapshot and compares both maps
using exact integer arithmetic. It enumerates all 524,288 image words and computes
both kernel spectra through the MacWilliams transform.

Of the 128 columns, 120 differ. The first difference is at zero-based index 3:
the imported column is `0x62907`; the selected column is `0x76607`.
Nevertheless, both image spectra are exactly

```text
weight:     0    48      56      64      72    80  128
count:      1  5040  110848  292510  110848  5040    1
```

Both kernel spectra also agree exactly; their minimum nonzero weight is 6.
The imported generator matrix and its transpose satisfy `BA=0`.
These checks do not establish equivalence under coordinate permutations or
state-basis changes; the checker does not attempt that search.

The existing grid engines consume constituent spectra, rather than ordered
columns. These constituent inputs therefore do not distinguish the two maps.
Transferring a complete result still requires matching the outer construction,
inner recurrence, boundary conditions, and setup distribution.

From the bridge directory, replay the comparison with a new output filename:

```powershell
C:/Python314/python.exe -B -m unittest test_audit_imported_rm2sub -v
C:/Python314/python.exe -B audit_imported_rm2sub.py --output ../../out/implementation_validation/constituent_audit_v2.json
```

The original receipt is `constituent_audit_v1.json`. Four checker unit tests pass.
The receipt deliberately leaves `full_margin_bits` null.

## Remaining construction differences

The imported BCH encoder's executed forward map has the following stages:

```text
BCH row encoder -> parity fanout (31 sources, 33 targets)
                -> factored route -> RM2Sub-S19
```

The fanout is visible in both the production transpose and the staged oracle.
It changes outer coordinates by XOR, not merely by permutation. Consequently,
the existing bare-BCH weight bounds cannot be attached to this encoder without
a separate argument. Legacy field-checkpoint routines also exist in its header;
their presence does not mean the current public method executes them.

The imported BCH geometry is fixed at K=2^20. It does not provide the K=2^16
and K=2^18 grid cells, or the other shortlisted inner map t64_s16.
The public setup expands short seeds deterministically; our analysis instead
uses specified independent random permutations and field multipliers.
Correctness for a seeded realization does not prove a margin for that realization.

The Golay package has a different outer code and shortening rule. It remains a
separate comparison implementation, with its documented BA conditioning issue.
Its passing tests do not certify a BCH candidate.

## Recommended next implementation step

Build a separate bare-BCH adapter, retaining these imported baselines unchanged.
Reuse the packed route and unrolled S19 kernel, but omit the parity fanout through
a compile-time choice. Validate the complete transpose against a staged oracle.
Support K=2^16, 2^18, and 2^20 with explicit geometry and buffer checks.

Before attaching numerical results, authenticate the BCH generator and document
the recurrence and setup assumptions. The S19 spectral comparison is already
available; there is no need to rewrite that hot kernel merely to match column bytes.
If an exact selected-map implementation is desired, regenerate its optimized
tables and circuit together rather than replacing only the column table.

Separately, compare explicit routing schedules across GCC and MSVC to resolve the
historical checksum discrepancy. Then benchmark matched implementations serially.
The t64_s16 adapter follows after the S19 path has a reproducible correctness test.
