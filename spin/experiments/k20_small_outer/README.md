# K20: exact small-outer/s22 comparison

**Update:** the [follow-up optimization](OPTIMIZATION.md) brings the exact
small-outer/s22 code to **2.962029 ms**, versus **3.299650 ms** for the larger
outer in a new matched holdout. Its 58.750377-bit certificate is unchanged.
The [library promotion](../../CERTIFIED22.md) records the subsequent public-API
validation at K18 and K20.
The remainder of this file preserves the initial, unoptimized comparison.

The certified small outer scales to K=2^20 in the proof, but its unchanged
K18 kernel is slower there. In the matched 2026-10-03 campaign, it takes
**5.213749 ms**, versus **3.295944 ms** for the current larger-outer library
profile: **58.19% more time**. Keep the larger outer as the K20 default.

| Configuration | Outer group, input/output bits | State bits | Whole margin at 10% distance | Precomputed transpose |
|---|---|---:|---:|---:|
| Exact small-outer/s22 research kernel | 128 / 256 | 22 | 58.750377 bits | 5.213749 ms |
| Library `PacketRsT64S20` | 256 / 512 | 20 | 68.103757 bits | 3.295944 ms |

Both have rate 1/2 and process 128-bit XOR elements. These margins bound
setup failure under each specified independent ensemble. They are not
measurements of the distance of the benchmark seeds. The
[small-outer certificate](../../../research/workstreams/k20_small_outer_screen/README.md)
and [baseline index](../../../research/CERTIFIED_BASELINES.md) identify the
proofs and their exact maps.

## What was measured

The candidate directly calls the unchanged `k18_certified22` setup and
`reverseRouteFast`, then `fieldLoopSharedParity`. It does not substitute
another map or omit the transformed feedback factor. The baseline calls
the current library through `spin::Code` with an explicit AVX512 backend.
No production source, candidate kernel, or certificate was changed here.

The common templated timing loop uses identical input fill, input allocation,
64-byte buffer offset, five warmups, and 1001 repeated in-place calls.
Setup, allocation, initialization, and checksums are outside the timed interval.
Each process requests large pages for both input and workspace. Each allocation
is 34 MiB, including rounding. The request gives 2 MiB alignment and best-effort
page advice; the driver does not measure actual huge-page backing.

Peach is a Ryzen 9 7950X. All runs used CPU 15, a 4.5 GHz maximum with boost
disabled, GCC 15.2, Release optimization, and `znver4` tuning. All four shared
benchmark locks were held; no benchmarks ran concurrently. Each row below
is a fresh process. Times are medians of 1001 calls, in microseconds.

| Order | Seed | Large s20 | Small s22 |
|---:|---:|---:|---:|
| 1 | 41 | 3311.994 | |
| 2 | 41 | | 5219.820 |
| 3 | 41 | | 5203.810 |
| 4 | 41 | 3256.179 | |
| 5 | 113 | | 5208.389 |
| 6 | 113 | 3303.488 | |
| 7 | 113 | 3288.400 | |
| 8 | 113 | | 5219.109 |

The headline uses the median of each construction's four process medians.
An earlier identical campaign gave 3.318757 and 5.219408 ms. Between campaigns,
only the scratch-allocation reporting and an error-message typo changed.
The timing loop and kernels did not change.

## Correctness and phase diagnosis

All three K20 checks pass, for seeds 1, 41, and 113. The existing checker
compares routed data and padding against the literal-map scalar oracle,
checks final outputs at four alignments, and tests guards, input preservation,
the forward/transpose adjoint relation, and in-place tail preservation.
The map generator also passes its map, basis-inverse, feedback, and conjugation
checks. The common driver's seed-41 checksum after 1001 calls matches the
unchanged original driver: `995a918cb290eca6`.

Separate 301-call diagnostics split the small kernel as follows. These use
an extra clock between phases and are not the headline measurement.

| Seed | Inner and route | Outer | Full call |
|---:|---:|---:|---:|
| 41 | 3.358458 ms | 1.864033 ms | 5.223633 ms |
| 113 | 3.355642 ms | 1.857932 ms | 5.214426 ms |

Thus inner/routing alone costs about as much as the complete larger-outer
encoder. The small kernel uses cached route and output stores; the production
K20 kernel uses non-temporal stores. This comparison does not isolate store
policy as the cause. It also does not establish that the small construction
cannot be optimized at K20.

The next bounded experiment, if further K20 tuning is wanted, is an A/B test
of streaming stores in the exact small-outer kernel, with the required fences
and independent scalar checks. That would preserve the code and certificate.
Do not promote the current slower implementation or change the proof to chase
this implementation cost. K18's existing result remains separate and unchanged.

## Reproduction and retained evidence

From the repository root on the benchmark host:

```sh
cmake -S spin -B tmp/k20-lib -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_CXX_COMPILER=g++-15 -DSPIN_TUNE=znver4 -DSPIN_BUILD_TESTS=OFF
cmake --build tmp/k20-lib -j3
cmake -S spin/experiments/k20_small_outer -B tmp/k20-small \
  -Dspin_DIR="$PWD/tmp/k20-lib" -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_CXX_COMPILER=g++-15
cmake --build tmp/k20-small -j3
ctest --test-dir tmp/k20-small --output-on-failure -j1
```

Acquire the shared benchmark locks before timing. Run only one process at a
time. Use the order and seeds in the table, replacing `small` by `large`
for the baseline:

```sh
taskset -c 15 tmp/k20-small/spin_k20_small_outer small 1048576 41 1001 1
taskset -c 15 tmp/k20-small/spin_k20_exact22_checks 1048576 41 301 1 1
```

The measured source snapshot is retained in ignored
`tmp/k20-small-benchmark-source-20261003-v2.tar.gz`, SHA256
`2d6e1ee1e9248bfdd9254d941b906aff79cec52882331fbdb7ffcf676a138b27`.
This captures the working-tree library sources used for the comparison.
The unchanged candidate `K18Fast.cpp` has SHA256
`a56e2f4be564b36a316ad01f411398e8d25f63590405fde2ee0289d60d696148`.
The final driver has SHA256
`b2bceee8b1e70e1ef1b9f6f66f00fc3545fe62434f7b33f5dfccfb0b7fc79f03`.

Raw process summaries remain ignored:

- `tmp/k20-small-campaign-v2.csv`, SHA256 `b6841133615fa474cab472151b3f48be198cdf87ec8c940af3a349258e230b6d`.
- `tmp/k20-small-campaign-v1.csv`, SHA256 `4993133b0c8f06ec09fd57dff1e08abbfa816fd291f9ef4c5d6d33c65c38715f`.
- `tmp/k20-small-diagnostic-v1.csv`, SHA256 `2cea219ff89016d99279b6f2edcbd310acb0043d5f61caab1121725ddb81dc3e`.

The remote build and logs are in `/tmp/spin-k20-small-certified-iw7EOr`.
No library profile was promoted, and no paper timing was changed.
