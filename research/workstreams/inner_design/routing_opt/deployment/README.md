# Optional workspace routing optimization

The supported encoder now offers an opt-in routing optimization. It uses Linux
large-page advice for owned workspace buffers and write-prefetch for the tile
scatter. It does not change the inner map or coordinate permutation.

Enable it when configuring the supported build:

```sh
cmake -S workstreams/bare_bch_rm2sub -B out/spin-pages -DCMAKE_BUILD_TYPE=Release -DSPIN_ARCH=znver4 -DSPIN_WORKSPACE_ROUTING_OPT=ON
cmake --build out/spin-pages -j3
ctest --test-dir out/spin-pages --output-on-failure
```

The option defaults to OFF. The public API and object layout are unchanged.

## Policy and fallback

Only Linux quarter-rate instances with K >= 2^18 select the optimized path.
K=2^16, half-rate instances, and non-Linux builds use the original kernel
instantiation. Selection occurs once per encoding operation, outside the inner
and scatter loops. A compile-time map tag instantiates the tuned kernel without
changing any map coefficients or inner operations.

Workspace construction requests `MADV_HUGEPAGE` and, where supported by the
headers and kernel, `MADV_COLLAPSE`. These best-effort requests cover only whole
pages inside the two owned vectors. Failure is nonfatal, retains ordinary
allocation semantics, and does not change output. Write-prefetch remains valid
with ordinary pages. The helper preserves `errno` and emits no diagnostics.

No system-wide settings, caller-owned buffers, allocation sizes, or routing
schedules change. There is no hot-path allocation, virtual dispatch, or
per-coordinate policy test. Synchronous page preparation happens during
workspace construction and can add setup latency. The benchmark excludes it;
its `setup_ms` field also excludes workspace construction.

## Measurements

Each cell is the median of three process medians, each containing 101 in-place
encoding calls after three warmups. Input is reused without copying between
calls. On/off runs are interleaved and serialized with the shared benchmark
lock. The host is Peach, Ryzen 9 7950X, CPU 15, Linux, GCC 15.2, znver4.

| Inner implementation | K | OFF (ms) | ON (ms) | Time reduction |
|---|---:|---:|---:|---:|
| Supported | 2^16 | 0.76658 | 0.77371 | -0.93% |
| Supported | 2^18 | 3.77272 | 3.65925 | 3.01% |
| Supported | 2^20 | 17.64652 | 16.54631 | 6.23% |
| Experimental asymmetric | 2^16 | 0.64254 | 0.64886 | -0.98% |
| Experimental asymmetric | 2^18 | 3.33371 | 3.26677 | 2.01% |
| Experimental asymmetric | 2^20 | 16.09406 | 15.25434 | 5.22% |

The optional kernel is not activated at K=2^16; small timing differences between
the builds remain. No speedup is claimed at that size. Larger-size gains are
host-dependent and do not guarantee large-page availability on another machine.
The experimental inner remains separate: the supported build option does not
promote it or select its coefficients.

Retained setup and workspace sizes are identical between ON and OFF for each
cell. Output hashes agree for every corresponding run.

## Validation and provenance

The deployment build has twelve passing test targets:

- Supported ON/OFF quarter-rate and half-rate correctness suites.
- Experimental ON/OFF quarter-rate suites, plus forced unavailable-page and
  simulated portable builds.
- Four page-helper tests: disabled, enabled, failed-advice, and portable modes.

The normal supported CMake build was also configured with the option ON;
its two public correctness targets pass.

The encoding suites compare against independent dense oracles at K=2^16,
2^18, and 2^20. They cover both routing representations, alternate tiles,
compaction, in-place suffix preservation, and independent setups. Page tests
check range rounding, invalid page sizes, overflow rejection, unchanged buffer
contents, preserved `errno`, and no advice for tiny spans. Portable behavior
is simulated on Linux; no native non-Linux validation is claimed.

An initial integration put the size branch inside the inner kernel. Its
experimental small-K build regressed. The final version dispatches between
separate template instantiations before entering the kernel. Initial receipts
are retained under `measurements/deployment/first-integration/` and are not used
for the final table.

The prior balanced distance certificate and generated candidates are unchanged.
Its integrity checker still passes. The checker resolves only the historical
`Spin.cpp` and CMake build inputs to exact digest-named snapshots when the live
supported files differ. Original manifest hashes remain unchanged. Numerical
inputs and live certified candidates do not permit that substitution.
See [source snapshots](../../source_snapshots/README.md).

## Reproduce the deployment comparison

From the repository root:

```sh
python workstreams/inner_design/routing_opt/deployment/generate.py
cmake -S workstreams/inner_design/routing_opt/deployment -B deployment-build -DCMAKE_BUILD_TYPE=Release
cmake --build deployment-build -j3
bash workstreams/inner_design/routing_opt/deployment/run.sh "$PWD"
cmake -S workstreams/bare_bch_rm2sub -B deployment-public -DCMAKE_BUILD_TYPE=Release -DSPIN_WORKSPACE_ROUTING_OPT=ON
cmake --build deployment-public -j3
ctest --test-dir deployment-public -V -j1 > routing-opt/deployment/public-correctness.log
```

The runner writes receipts under `routing-opt/deployment/` in the supplied
build root. Copy them to `workstreams/inner_design/routing_opt/measurements/deployment/`,
then run `deployment/summarize.py` from this directory's parent. The summarizer
checks source hashes, test completion, output hashes, and memory-size equality.
It writes `RESULT.json`; it does not produce a new distance certificate.

The shared implementation is `workstreams/bare_bch_rm2sub/WorkspaceRouting.h`.
`generate.py` applies the same routing hook to the unchanged experimental inner.
No benchmark should run concurrently with another benchmark.

Next: retain this as an optional deployment setting and check a second host
before considering broader defaults. Finish the asymmetric inner's outward
certificate separately from these routing changes.
