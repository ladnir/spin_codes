# Fresh-code target: within 10% of precomputed encoding

The target includes fresh instance parameters, fresh IMT masks, encoding,
and destruction of the per-code plan. The bank and scratch persist across
calls. The reference retains one fully precomputed code and its workspace.
Both flows process K=2^18 output elements, N=2^19 input elements, and
128-bit XOR values using BCH [256,128] and T128S19.

**Latest result:** the screened rotated family meets the median target in
three longer serial comparisons: 9.68%, 9.87%, and 9.56% overhead. See
[TARGET_RESULT.md](TARGET_RESULT.md) for the two flows, final measurements,
reproduction script, memory costs, and heuristic tradeoffs. The original
fixed-code path remains unchanged. No full-code certificate is transferred.
The remainder of this file records the earlier optimization stages.

For the original tested small-bank family, the best matched comparison is 2.211 ms
fresh versus 1.490 ms precomputed: 48.4% overhead. The 10% target is
1.639 ms. Another approximately 0.572 ms must be removed from the fresh flow.

## Changes that helped

`Rolling.h` retains the permutation family from `README.md`. Each region
contains 2,048 inner coordinates. Instead of resolving every address within
the inner loop, the encoder prepares one region's addresses at a time.

1. Stream the row parameters and compute one destination coordinate per row.
2. Apply the region permutation to that temporary coordinate array.
3. Scatter the region's inner outputs using the resulting addresses.

The two temporary arrays occupy 16 KiB. Neither setup nor encoding retains
a full 2 MiB route. Row parameters have a SIMD-friendly representation;
eight-lane AVX2 gathers evaluate the wrappers. For the four-row BCH kernel,
route preparation directly produces its interleaved destination layout.
The emitter therefore avoids another conversion on every output element.

The reusable scratch receives the same owned-page huge-page preparation
used by the production workspace. This is a one-time preparation cost.
The measured experiment also compiles its inline inner with the AVX-512
flags used by the existing fixed-code kernel. Public ISA dispatch is unchanged.

The family, route seeds, mask seeds, inner arithmetic, and BCH arithmetic
are unchanged. Production defaults and libOTe integration are untouched.

## Matched measurements

Ryzen 7950X, CPU 15, GCC 15.2, Linux. Commands below identify the build flags.
Each process initializes a new input for each of 15 calls and discards the
first two timings. The table reports the median of three process medians.
Order reverses on the middle pass; benchmark processes never overlap.
Input initialization, checksumming, bank creation, and reusable scratch
preparation are outside timing. Fresh plans use new route and mask seeds
on every call. Fixed plans retain seeds 17 and 29.

| Flow | Per-code setup | Encode | Total |
|---|---:|---:|---:|
| Fixed, original family | amortized | 1.490082 ms | 1.490082 ms |
| Fixed, bank-derived family | amortized | 1.483902 ms | 1.483902 ms |
| Fresh, SIMD region routing (`vector`) | 0.040446 ms | 2.170453 ms | 2.211419 ms |
| Fresh, region check outside unrolled emitter (`epoch`) | 0.040516 ms | 2.196331 ms | 2.239662 ms |

Independent medians need not add exactly. The fresh winner's process totals
were 2.205688, 2.211419, and 2.237889 ms. Original-family fixed totals were
1.473822, 1.490082, and 1.504429 ms. These results are retained remotely in
`results-bank-optimize-epoch`, outside source control.

The earlier implementation measured 3.592 ms fresh in `FLOWS.md`. The new
result is about 38% lower, but the benchmark protocol and scratch preparation
also changed. Use the matched table above to assess the 10% target.

Other experiments did not improve the winner:

- Separating all inner computation from scatter added a sequential pass;
  the packed-address version measured about 2.376 ms.
- Prefetching destinations 32 elements ahead measured about 2.337 ms.
- Sixteen-lane AVX-512 routing measured 2.246 ms, versus a 1.475 ms fixed
  control in that comparison. It is optional, not the chosen routing path.

The first two comparisons reside in `results-bank-optimize-prefetch`;
the last resides in `results-bank-optimize-wider`. All preserve the family.
They suggest that wider arithmetic alone will not close the remaining gap.

For bank size 16, packed tables occupy 68 KiB plus five gather-padding bytes.
SIMD row parameters add 16 KiB to the previous 59 KiB per-code plan, which
includes masks. Reusable value scratch remains 8 MiB. Size 64 is covered
by routing correctness tests, but these optimized timing results use size 16.

## Validation and scope

Full encoders compare their first K output elements with the existing
precomputed kernel using the same materialized map and two code seeds.
Every timed fresh implementation produced checksum `2d0de1a66683cfce`.
Split scheduling may overwrite the unused N-K suffix; that suffix is not
part of the transposed encoder's output.

`test.cpp` compares scalar, batched, region, and SIMD region addresses
against the independent `Bank` implementation. It enumerates all N addresses
for bank sizes 16 and 64, with both ordinary and four-row-packed destinations.
It checks bijection, inverses, and the documented XOR-family identities.
Focused ASan/UBSan runs passed locally with AVX2 and remotely with AVX-512
routing. These tests do not constitute a sanitizer rebuild of the library.

These are implementation checks, not distance certificates. The unchanged
family retains the empirical evidence described in `README.md`. No new
statistical evidence or proof follows merely from its faster evaluation.

## Reproduce

Portable AVX2 experiment and same-map correctness checks:

```sh
cmake -S spin -B out/bank-opt -DCMAKE_BUILD_TYPE=Release \
  -DSPIN_BUILD_EXPERIMENTS=ON -DSPIN_BUILD_TESTS=ON
cmake --build out/bank-opt --target spin_bank_optimize spin_bank_test -j2
out/bank-opt/spin_bank_test
ctest --test-dir out/bank-opt \
  -R '^spin_bank_(rolling|split|vector|vector-split|epoch)$' --output-on-failure -j1
```

The following build is host-specific: run it only on an AVX-512F/VL machine.
It assumes `build/libspin.a` is an existing Release library with its
AVX-512 backend enabled. `COMPARISON_ROOT` contains that build directory.

```sh
g++ -O3 -DNDEBUG -std=c++20 -mavx2 -mavx512f -mavx512vl \
  -mtune=znver4 -DSPIN_WORKSPACE_ROUTING_OPT=1 -DSPIN_BCH_AVX512=1 \
  -Ispin/include spin/experiments/permutation_bank/optimize.cpp \
  build/libspin.a -o build/spin_bank_optimize512
SPIN_COMPARE_MODES=fixed,fixed-bank,vector,epoch \
bash spin/experiments/permutation_bank/compare_optimize.sh \
  COMPARISON_ROOT 15 spin_bank_optimize512 epoch
```

Add `-DSPIN_BANK_ROUTE512=1` to reproduce the wider-route experiment.
The comparison script acquires the three shared benchmark locks and runs
each process sequentially. No benchmark data is committed.

## Follow-up

The [composed-bank experiment](COMPOSED.md) explores a richer persistent bank
that precomputes more of the composition,
reducing the fresh route's dependent lookups and coordinate arithmetic.
This may trade bank memory for per-code speed and may change the family.
Before adopting a changed family, define its sampled parameters, check
bijection and the one-bit-per-row-per-region property, then rerun fixed-bank
pair, four-point, and bank-adaptive support screens. Do not treat the existing
screens as evidence for a different family. Retain the original fixed-code
baseline and include fresh setup in every candidate's total.
