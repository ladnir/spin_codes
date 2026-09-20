# SPIN with BCH outer and IMT inner: half-rate implementation and performance

This implementation uses the BCH [256,128] outer, randomized bit-transpose
permutation, and **Independent-Map Transvection (IMT)** inner with weight-five
feedback. [IMT](../../../../IMT.md) names only the inner. Historical code and
target identifiers remain unchanged to preserve verification records.

The isolated, optimized transposed encoder takes **10.155 ms** at K=2^20
on Peach. The current supported RM2Sub build takes 11.169 ms in the same
measurement series. The new implementation reduces time by **9.08%**.
Giving RM2Sub the large-page and write-prefetch optimization reduces its
time to 10.658 ms; the new implementation remains **4.72%** faster than
that tuned control.

The selected target is `sparse_pages`: unrolled sparse emission, sparse
transvection updates, large-page advice for owned workspace, and the original
32-entry read-prefetch. It retains packed-24 routing and 2048-row tiles.
The supported library and its defaults are unchanged.

## Proof and implementation binding

This uses the exact IMT feedback map `weight5_seed0` from the [full finite certificate](../README.md),
not a newly screened map. The fixed BCH [256,128] outer has rate 1/2;
the inner preserves length. At K=2^20, N=2^21, t=128, s=19, the certificate
gives relative distance greater than 0.10 with 50.0620882264 margin bits.

`generate.py` authenticates the certificate, reconstructs the exact A and B
columns, and symbolically checks the expansion-feedback circuit and alternative
emission circuits. The transpose computes A^T from raw input, not emitted
output. The exact independent-map recurrence and zero terminal state are
unchanged. The existing paired AVX2 BCH circuit supplies the outer transpose.

`report.py` reruns the exact union verifier, checks the generated map columns,
reconstructs the exact systematic BCH generator rows, and authenticates the
compiled source hashes. The implementation result is in `RESULT.json`;
it supplements the historical numerical certificate without modifying it.

All twelve release tests pass. Nine encoder targets are compared against
the dense oracle at K=2^16,2^18,2^20, covering both routing layouts, in-place
suffix preservation, compaction, alternate tiles, and a second setup with
boundary impulses. Three additional tests check the forward/transpose identity
against a scalar forward model with separate A and B. Those three tests and
the two page-advised sparse encoder variants also pass AddressSanitizer and
UndefinedBehaviorSanitizer. Implementation tests at smaller K do not extend
the distance certificate to those lengths.

## Same-host timings

Each process performs three warmups and 101 timed calls. Calls encode in place,
without copying or resetting the input between calls. All timings are serial
and pinned to CPU 15 on Peach: Ryzen 9 7950X, Linux, GCC 15.2, Release -O3,
znver4/AVX2. The workload is 128 parallel binary instances. Code setup, workspace
construction, buffer initialization, and output hashing are excluded.
Large-page requests occur during workspace construction, outside the interval.

The table reports the median of process medians. Six primary variants have
two interleaved batches of three processes; the last three have one batch.
The runners vary ordering and use the shared benchmark lock plus process checks.

| Variant | Processes | Time (ms) |
|---|---:|---:|
| Supported RM2Sub | 6 | 11.169 |
| RM2Sub + page advice and write-prefetch | 6 | 10.658 |
| IMT, weight five, plain sparse | 6 | 10.488 |
| IMT + page advice, selected | 6 | **10.155** |
| IMT + page advice and write-prefetch | 6 | 10.210 |
| IMT + grouped tables and routing tuning | 6 | 10.207 |
| IMT + masked update and routing tuning | 3 | 10.182 |
| IMT + shared-XOR emission and routing tuning | 3 | 10.269 |
| IMT + 128-entry write-prefetch and page advice | 3 | 10.631 |

The near-10.2 ms variants overlap in run-to-run variation. The data supports
the simple selected implementation, not a robust ordering of every small
kernel variation. The initial repeat batch gave 10.237 versus 11.275 ms;
the confirmation batch gave 10.138 versus 11.055 ms for selected versus
supported. The headline uses both batches, rather than the fastest run.

Retained setup storage drops from 13,828,104 to 12,713,992 bytes. Both use
41,943,040 bytes (40 MiB) of workspace at the selected tile size. Page advice
is best-effort and does not change allocation sizes or system-wide settings.
The benchmark's `setup_ms` field excludes workspace construction, so it does
not include the cost of synchronous page collapse requests.

## Quarter-rate optimizations carried over or tested

- Retained compile-time unrolling, fixed-size SIMD scratch, the pruned AVX2
  zeta transform, and the synthesized A^T finish circuit. No allocations,
  erased callbacks, or dynamic dispatch were added to the encoding loop.
- Retained sparse transvection updates; tested the masked alternative.
- Applied the owned-workspace large-page helper to the half-rate candidate
  and to a matched old-inner control, without changing routing coordinates.
- Tested 32-entry and 128-entry write-prefetch. The selected half-rate kernel
  retains read-prefetch; the quarter-rate write-prefetch winner did not
  establish an extra improvement here.
- Tested optimized nibble grouping and global XOR sharing for B^T emission.
  XOR sharing reduces the symbolic emission count from 512 to 274, but did
  not deliver a repeatable end-to-end win.
- Swept packed-24 and 32-bit layouts with 256, 512, 1024, 2048, and 4096-row
  tiles for four variants: forty cells, each with 31 timed calls. Every
  variant's best cell used packed-24 and 2048-row tiles. Output hashes agreed
  across all equivalent layouts and tiles.
- Reused the existing two-row AVX2 BCH-256 circuit. The quarter-rate BCH-128
  circuit is a different outer and cannot replace it. No permutation law,
  block grouping, or other construction parameter was changed for speed.

## Reproduction and files

`AsymmetricMap.h`, `Weight5Inner.h`, and `Weight5Spin.cpp` are isolated sources.
`BaselineTuned.cpp` is the old-inner routing control. `CMakeLists.txt` selects
all performance policies at compile time. `inner_identity.cpp` supplies the
independent scalar forward check. `correctness.cpp` selects the half-rate
t128_s19 cell of the existing complete-encoder test suite.

On Linux, from the repository root:

```sh
cmake -S workstreams/inner_design/asymmetric/bch256/weight5/implementation -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j3
bash workstreams/inner_design/asymmetric/bch256/weight5/implementation/run.sh "$PWD"
bash workstreams/inner_design/asymmetric/bch256/weight5/implementation/tune.sh "$PWD"
bash workstreams/inner_design/asymmetric/bch256/weight5/implementation/final.sh "$PWD"
bash workstreams/inner_design/asymmetric/bch256/weight5/implementation/confirm.sh "$PWD"
```

Run the scripts sequentially, never concurrently. Build commands assume an
AVX2-capable x86 Linux host; CPU 15 and znver4 are explicit measurement settings.
The optional `screen_final.sh` repeats the short screen against the final
binaries. Timings alone never replace the correctness tests.

For sanitizers, configure another build with
`-DCMAKE_CXX_FLAGS="-fsanitize=address,undefined -fno-omit-frame-pointer"` and
`-DCMAKE_EXE_LINKER_FLAGS="-fsanitize=address,undefined"`. Build and run
`sparse_pages_test`, `sparse_tuned_test`, and `inner_identity_0_test` through
`inner_identity_2_test`, without a concurrent benchmark.

Copy the runner's `measurements/` directory here. Retain the four-test
sanitizer log as `sanitizer.log` and the selected-only log as
`sanitizer-selected.log`. With the local proof inputs present, run:

```text
python workstreams/inner_design/asymmetric/bch256/weight5/implementation/report.py
```

Raw records and `RESULT.json` remain local and uncommitted. The report binds
source, compiler-flag, binary-hash, correctness, and timing records; the
code-and-documentation checkpoint alone is not the complete data bundle.

Next: integrate the certified quarter-rate and half-rate profiles behind an
explicit opt-in, preserving the old baseline. Smaller half-rate certificates
and the asymptotic argument remain separate gates before replacing RM2Sub
across the family. A complete optimized forward encoder still needs its own
implementation and end-to-end measurement.
