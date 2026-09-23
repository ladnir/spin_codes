# Fresh-code encoding within 10%: experimental result

The screened heuristic prototype meets the **median** performance target at
K=2^18, N=2^19, with 128-bit elements and BCH [256,128], T128S19.
Three serial measurements of the rotated family gave 9.68%, 9.87%, and 9.56%
overhead over the unchanged original-family precomputed encoder. This is not a
worst-run latency guarantee or a transferred SPIN distance certificate.

## The two flows

The fixed flow retains a fully precomputed code and workspace.
`PrecomputedCode` uses the original permutation family unless explicitly
given a materialized route for a correctness comparison.

The fresh flow retains one bank and an 8 MiB encoding workspace. Each call
generates fresh route parameters and IMT masks, encodes, and destroys those
per-code objects. It prepares only 2,048 route addresses at a time: an 8 KiB
temporary buffer, not a full per-code route. Bank preparation is amortized.

The selected experimental type is
`RollingCode<1,false,false,0,true,ComposedRouting<1,true,true,true,true>>`.
The final `true` enables rotation. Its bank must outlive every code object;
do not move or mutate the bank during use. Each concurrent encoding needs
its own workspace. None of these experiments changes public library or
libOTe defaults.

## What is sampled

The bank contains fixed independent row permutations and one independent
row-label permutation per stored region. See [COMPOSED.md](COMPOSED.md)
for the table construction. A fresh code samples a shared region permutation
and an independent cyclic position shift in each region. In each outer row,
it independently samples the byte permutation

```
F(c) = ((a * rotr8(c,d) + b) mod 256) XOR e,
```

where a is odd, b and e are bytes, and d ranges from 0 through 7. The affine
alternative omits rotation. The seed-driven implementation uses the existing
word generator. The mathematical marginal statements in the quality report
assume the specified independent uniform draws.

All routes remain bijective, with exactly one coordinate from each outer row
in each region. The IMT mask generator preserves the previous generator's
outputs and rejection behavior for the same seed.

## Performance

Ryzen 7950X, Linux, CPU 15, GCC 15.2.0, Release library tuned for znver4.
The native experiment requires AVX2, AVX-512F/VL/DQ/BW and VPOPCNTDQ.
Each process performs 101 encodes and discards 20 warmups. Each fresh encode
uses new route and mask seeds. Processes run serially; comparison order is
reversed on alternating passes. Entries below are medians of process medians.

| Comparison | Fixed total | Fresh setup | Fresh encode | Fresh total | Overhead |
|---|---:|---:|---:|---:|---:|
| Rotated, 9 passes | 1.489170 ms | 0.006552 ms | 1.626687 ms | 1.633289 ms | 9.68% |
| Rotated, independent 15 passes | 1.486345 ms | 0.006542 ms | 1.626557 ms | 1.633009 ms | 9.87% |
| Rotated, clean-build script, 15 passes | 1.490282 ms | 0.006532 ms | 1.626266 ms | 1.632768 ms | 9.56% |
| Affine, 9 passes, 512-bit route arithmetic | 1.484282 ms | 0.006472 ms | 1.606199 ms | 1.612932 ms | 8.67% |
| Affine, independent 15 passes, same build | 1.486405 ms | 0.006482 ms | 1.618131 ms | 1.624663 ms | 9.30% |

The rotated confirmation's median paired ratio was 1.096575; individual
paired ratios ranged from 1.082661 to 1.128676. Some passes exceed 10%.
The affine variant in the selected rotated build measured 8.98% overhead
in its nine-pass comparison and 8.77% in the clean-build script run.
These are machine-specific measurements.

Fresh total includes parameter and mask generation, encoding, and destruction.
Both flows exclude input generation and checksumming. Bank construction,
persistent workspace allocation, and owned-page preparation are outside the
per-call total. Medians of components need not sum to the median total.

Separate setup measurements, using 15 samples and discarding two warmups:

| Component | Median time | Retained storage |
|---|---:|---:|
| Reusable indexed bank | 5.354 ms | 4 MiB |
| Fresh affine route parameters | 3.036 microseconds | 13,320 bytes |
| Fresh IMT masks | 2.935 microseconds | 32 KiB |
| Original fixed plan plus workspace | 10.785 ms | 20,512,780 bytes |

The bank retains a 2 MiB ordinary array as a scalar reference and a 2 MiB
indexed array for the fast path. Setup timings exclude subsequent page advice.
The fresh encoding workspace is 8 MiB, in addition to the bank and per-code
objects. The setup-component measurement used the affine variant; rotated
end-to-end setup is measured directly in the performance table.

## Why the final optimization helps

The routing buffer now stores byte offsets. SIMD routing scales 8 destination
indices together, so each scalar 128-bit store no longer scales its index by 16.
This preserves the sampled map and its seeds. The selected rotated build also
uses four unrolled 256-bit vectors, precomputed row indices, repacked row keys,
table lookahead, and the exact batched IMT mask generator.

Smaller routing chunks, output prefetching, and forced full-region inlining
did not improve the result. The older, short 15-encode runs sometimes crossed
the target in either direction. The final comparisons use the longer protocol
above; the target was not declared achieved from those short-run fluctuations.

## What the screening establishes—and does not

[QUALITY.md](QUALITY.md) records exact finite enumerations and their limitations.
Cyclic row shifts alone admit selected bank-adaptive 32-support alignments
with probability 2^-16. Odd multiply/add/XOR improves the tested cases, but
preserves addition by 128; supports adapted to this invariant aligned with
probability 2^-36 in two enumerated banks. Adding rotation lowers those same
two probabilities to 2^-42. These are selected-support probabilities, not
security margins, worst-case bounds, or probabilities of low-distance SPIN codes.
The selected supports have not been shown to be BCH codewords.

The fresh family is distinguishable from uniform permutations. In particular,
the fixed region tables retain circular row order across instances. No test
here establishes full-code distance or the security properties required by
an SSD application. Using this family in that setting requires an additional
heuristic assumption. Rotation removes one observed invariant; it does not
make the remaining assumptions unnecessary.

Under the specified independent sampling, the original one-active-row routing
marginal is preserved for every fixed bank and fixed row support. That supports
the same one-row first-moment calculation, but does not transfer multi-row
bounds. Retain the original-family precomputed flow wherever its existing
analysis is required. The performance objective is met without promoting the
heuristic flow to a certified or production default.

## Verification and reproduction

The 24 portable bank tests pass. Full-encoder tests compare all K outputs
against a materialized-route reference for two seeds before every timing run.
Routing tests enumerate all N addresses for 32 seeds, compare scalar and SIMD
results, verify bijection and per-region row occupancy, and check both layouts.
They also check byte offsets and 128/512/1,024-address chunks against full-region
preparation. Focused AVX2 and native AVX-512 ASan/UBSan routing tests pass.
Mask tests compare 2,048 seeds, including 15 rejection fallbacks.

From the repository root on a compatible Linux host:

```sh
bash spin/experiments/permutation_bank/reproduce_target.sh \
  "$PWD/out/permutation-bank-target" 15
```

The script was executed successfully from a new build directory. It builds
the library and native experiment, checks routing and masks,
checks both encoder variants, and runs matched comparisons and setup timings.
It acquires the three shared benchmark locks through `compare_optimize.sh`.
Do not run another benchmark concurrently. The optional second argument is
the CPU number; the default is 15. Timing data stay in the selected output
directory and are not source artifacts.

Recorded result directories on the measurement host end in `bytes256`,
`bytes256-confirm`, `bytes`, and `bytes-confirm`, respectively. Setup results
end in `setup-final`. The clean-build run is under
`out/target-reproduce/results-bank-optimize-target` on the same measurement
host. Further work should evaluate the heuristic family for
the intended PCG use before changing a customer-facing default.
