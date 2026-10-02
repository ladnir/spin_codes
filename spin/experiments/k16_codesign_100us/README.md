# K16 co-design benchmark

The [frozen research checkpoint](checkpoint/README.md) pins the certified
winner and its proof sources, and records the bounded K18/K20 scaling probe.
It does not promote this construction into the library.

The winning research configuration is **mode 52**: a 128-to-256-bit RS16
outer group, four-bit packet routing, and the certified t64/s15 paired inner.
Precomputed transposed encoding of 128-bit elements at K=65,536 measures
**0.098910 ms**, with a separate holdout result of **0.098950 ms** on Peach.
Its complete 10% distance certificate has **62.04117 bits** of setup-failure
margin. See the [construction, proof, and measurement report](../../../research/workstreams/k16_codesign_100us/README.md).

This directory is an isolated research harness, not a production dispatch.
It links an existing SPIN package and preserves the retained RS sources.
Setup, buffers, and matrix preprocessing are outside the timed call. The
complete inner/routing and outer operations are inside it.

## Build and check

Use a supported x86 host with AVX-512, VBMI, and GFNI:

```sh
cmake -S spin/experiments/k16_codesign_100us -B out/k16-codesign \
  -Dspin_DIR=/absolute/path/to/spin/build/package -DCMAKE_BUILD_TYPE=Release
cmake --build out/k16-codesign -j4
ctest --test-dir out/k16-codesign --output-on-failure -j1
```

There are 319 compiled checks: the retained control plus 53 modes, three
lengths (2,048, 6,144, 65,536), and two seeds. Each process checks literal
scalar equality, forward/transpose adjoints, in-place operation, four
alignments, and guard bytes before benchmarking.

## Run serially

The arguments are mode, K, setup seed, calls, optional stage clocks, and
optional large-page preference. Zero calls performs checks only. Normal
pages and disabled stage clocks are the confirmed measurement settings.

On the shared Peach host, acquire every shared benchmark lock:

```sh
flock /tmp/prindal-addition-encoder-benchmark.lock \
flock /tmp/bare-spin-benchmark.lock \
flock /tmp/hypercat-benchmark.lock \
flock /tmp/hypercat-global-benchmark.lock \
taskset -c 15 out/k16-codesign/spin_k16_codesign 52 65536 1 2001 0 0
```

Never run concurrent benchmarks. Each invocation performs five warmups and
reports median, p10, and p90 call latency. The confirmation protocol compares
eight process medians, using four setup seeds in forward and reverse order.
Do not substitute the minimum sample for this result.

| Mode | Purpose |
|---|---|
| 0 | Retained smaller RS/s16 control |
| 10 | Original s16 map with native-field outer |
| 12 | First completely certified paired s16 map |
| 41 | Paired s16 with shuffle-friendly state basis and shared outer parity |
| 44 | Paired s15 with the corresponding layout |
| 48 / 50 | s16 / s15 with boundary steps peeled out of the loop |
| 51 / **52** | s16 / **s15** with fewer feedback shuffles |

Other modes retain negative controls and intermediate experiments. They are
not interchangeable parameter choices: some have only partial proof screens.
The final s15 scalar oracle uses literal compact coordinates; the fast kernel
uses a fixed state permutation with a zero SIMD slot. Native field16 maps
are used in the outer only. The inner samples independent GL15 matrices.

The final hot path is `reverseRoutePaired15Fold` followed by
`fieldLoopSharedParity`. Neither allocates. The nonaliasing contract concerns
input and routing scratch, which are disjoint even when the full encoder
overwrites its input after routing. Raw receipts and timings are not library
dependencies and remain ignored.
