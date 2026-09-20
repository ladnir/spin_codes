# Single-stream forward versus transpose

Follow-up: [FORWARD_SCHEDULE.md](FORWARD_SCHEDULE.md) improves these forward
results by scheduling the same BCH circuit on demand. This note preserves the
earlier matched measurements and routing diagnosis.

The optimized 128-bit forward encoder still trails transpose. Matched timing
shows a smaller gap than the earlier 16-stream workload suggested.
Removing an intermediate copy improves forward by another 5–6% at small K,
but does not close the gap.

## One 128-bit stream per call

Each element is 128 bits. Forward maps K elements to 2K elements; transpose
maps 2K elements to K elements. Both use the same code parameters and setup seeds.
The forward build enables four-row BCH batching. Transpose uses the current
standalone optimized library, not a historical imported transpose implementation.

| K | Forward, default tile | Forward, best tested tile | Transpose, separate output | Transpose, in place |
|---|---:|---:|---:|---:|
| 2^16 | 0.423 ms | 0.406 ms (64 rows) | 0.338 ms | 0.339 ms |
| 2^18 | 1.849 ms | 1.842 ms (256 rows) | 1.446 ms | 1.411 ms |
| 2^20 | 10.267 ms | 9.571 ms (256 rows) | 9.228 ms | 9.102 ms |

Tested tiles: default, 64, 256, and 512 rows. Defaults are 256 rows at K=2^16
and K=2^18, and 2048 rows at K=2^20. No default was changed.
The K=2^16 configuration is (t,s,rounds)=(64,12,2); the others use (128,19,1).

These are warm repeated-call times, not first-call or cold-cache latency.
There is one element stream, no packing, no column assembly, and no division
by a batch size. Setup, allocation, dense-reference checks, and output checks
are outside the timer. Separate-output modes reuse fixed inputs; in-place
transpose evolves its buffer, matching the established transpose harness.

## Where forward spends time

A separately instrumented build gives the following approximate breakdown.
Clock instrumentation can perturb scheduling; these figures diagnose stages
and do not replace the uninstrumented timings above.

| K and tile rows | BCH outer | Tile-to-bucket copy | Inner plus routed reads |
|---|---:|---:|---:|
| 2^16, 64 | 0.183 ms | 0.079 ms | 0.148 ms |
| 2^18, 256 | 0.759 ms | 0.373 ms | 0.739 ms |
| 2^20, 256 | 3.679 ms | 1.891 ms | 3.977 ms |

No single stage explains the whole gap. The intermediate copy costs about
20% of forward time. BCH encoding is also substantial, at about 40% for small K.
The inner timing includes scattered memory reads and does not isolate arithmetic.

## Direct-routing experiment

The direct variant writes four-row BCH outputs into the full workspace, then
lets the inner gather from those outputs. It removes the tile-to-bucket copy.
A setup-time table maps each inner position to its packed BCH output position.
The hot path has no new allocation, virtual dispatch, or runtime callback.
The encoded linear map, permutation distribution, and IMT parameters are unchanged.

| K | Tiled control | Direct forward | Time reduction |
|---|---:|---:|---:|
| 2^16 | 0.405 ms | 0.383 ms | 5.4% |
| 2^18 | 1.801 ms | 1.695 ms | 5.9% |

This is a fresh alternating-order A/B sweep. Its controls differ slightly from
the first sweep, so the percentages use only these paired controls.
Direct routing remains about 13% and 17% slower than separate-output transpose
at the two sizes, using the first sweep's transpose measurements.

The saved copy time exceeds the overall gain. This suggests that removing the
tiled layout loses some locality in the inner's reads; the current measurements
do not isolate that effect from instruction scheduling and index representation.

`SPIN_FORWARD_DIRECT=ON` enables the experiment for eligible AVX-512 forward
instances through K=2^18. Larger instances retain the tiled path. The option
defaults to OFF and requires `SPIN_FORWARD_FOUR=ON`. It adds a 32-bit route entry
per output element: 0.5 MiB at K=2^16 and 2 MiB at K=2^18. Existing workspace
allocation is unchanged. The standalone transpose target and kernels are untouched.
The direct option currently excludes the tiled profiler and prefetch experiment.

## Reproduction and checks

Host: Peach, Ryzen 9 7950X, GCC 15.2, CPU 15, `SPIN_TUNE=znver4`.
Each uninstrumented point uses three fresh processes and 101 measured calls per
process. The reported value is the median of process medians. Direction or
variant order reverses on the second repeat. All benchmarks acquire both shared
locks and reject concurrent benchmark processes. Builds finish before timing.

```sh
bash workstreams/spin_optimized/single_check.sh ROOT
bash workstreams/spin_optimized/single_direct_check.sh ROOT
bash workstreams/spin_optimized/single_direct_validate.sh ROOT
python3 workstreams/spin_optimized/single_summary.py
```

`ROOT/hypercat/native/spin` supplies the frozen upstream snapshot. Raw results
remain ignored under `measurements/single_forward` and `measurements/single_direct`.
Every timed point checks against a dense reference. Direct-routing release tests
also cover both layouts, compaction, small tiles, adjoint identities, and AVX2 fallback.
Both test suites also passed with AVX-512 detection disabled and under
AddressSanitizer/UndefinedBehaviorSanitizer. The object-code audit confirms ISA
isolation for the new direct kernel and the existing forward circuit.

## Next optimization

Keep single-stream 128-bit encoding as the primary target. Next inspect the
forward BCH circuit's scheduling and register pressure, then isolate the inner's
arithmetic from its gather cost. The direct variant is useful evidence, but its
additional route table should not become the default before the remaining gap
and the supported size range are understood.
