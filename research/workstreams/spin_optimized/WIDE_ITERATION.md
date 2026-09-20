# Wide forward encoding: on-demand BCH evaluation

The on-demand BCH schedule reduces 512-bit forward-encoding time by 13–19%
relative to the previous wide schedule. The 256-bit improvement is smaller.
Native wide inputs remain useful; packing separate streams can erase the gain.

These changes preserve the binary generator, routing, and IMT parameters.
They do not require a new distance certificate. The transposed kernels and
single-stream 128-bit kernels are unchanged by this iteration.

## Same workload, different element widths

Every entry below encodes **16 independent streams of 128-bit elements** under
one shared setup. Widths 128, 256, and 512 require 16, eight, and four calls.
The times are for the entire batch, not one stream or one encoder call.

Native wide input already interleaves the streams within each element.
Planar input stores streams separately and requires timed packing.
Encoding time includes BCH, routing, and IMT, but excludes packing and column
assembly. Total time includes both, when applicable, and excludes setup.
Column assembly transposes elements, not individual bits.

| K | 128-bit encode | Native 256-bit encode | Native 512-bit encode | 128-bit total | Native 256-bit total | Native 512-bit total |
|---|---:|---:|---:|---:|---:|---:|
| 2^16 | 7.035 ms | 5.897 ms | 5.418 ms | 8.517 ms | 7.391 ms | 6.956 ms |
| 2^18 | 30.620 ms | 28.755 ms | 28.755 ms | 38.066 ms | 35.661 ms | 35.499 ms |
| 2^20 | 144.830 ms | 127.887 ms | 124.488 ms | 174.267 ms | 156.370 ms | 155.093 ms |

All columns use the current optimized 128-bit control where that path applies.
K=2^16 uses (t,s,rounds)=(64,12,2); the other sizes use (128,19,1).
Tiles contain 256 outer rows, except the K=2^20 wide kernels, which use 512.

The 512-bit total is about 18%, 7%, and 11% lower than the 128-bit total.
At K=2^18 and 2^20, however, 256 and 512 bits have similar total times.
The 512-bit workspace is twice as large: at K=2^20 it occupies 136 MiB,
compared with 68 MiB for 256 bits. Input, output, setup, and assembly buffers
are additional. Thus 256 bits remains a reasonable memory-conscious choice.

Packing changes the decision:

| K | 128-bit total | Planar-to-256-bit total | Planar-to-512-bit total |
|---|---:|---:|---:|
| 2^16 | 8.517 ms | 8.109 ms | 8.117 ms |
| 2^18 | 38.066 ms | 39.851 ms | 41.841 ms |
| 2^20 | 174.267 ms | 178.408 ms | 177.318 ms |

At the larger two sizes, retaining separate 128-bit streams is faster than
packing them solely for this encoder. Applications that already produce native
wide elements do not pay that conversion cost.

## What changed

`wide_schedule.py` emits each BCH output after recursively emitting its
dependencies. Shared intermediates are computed once. This changes evaluation
order, not the circuit: the generated code still uses 2,901 XORs.
The generator checks all 256 output linear forms against the BCH matrix
before writing the reordered circuit. Both vector widths use this same schedule.

The previous schedule kept more intermediates live before storing outputs.
Reducing those lifetimes is the motivation; the timings below measure the
whole encoder, not register pressure or spills directly.

| K | Previous native 512-bit encode | On-demand encode | Time reduction |
|---|---:|---:|---:|
| 2^16 | 6.653 ms | 5.418 ms | 18.6% |
| 2^18 | 33.718 ms | 28.755 ms | 14.7% |
| 2^20 | 142.894 ms | 124.488 ms | 12.9% |

For 256-bit input, reductions were 2.6%, 0.3%, and 1.8% respectively.
The middle result is too small to interpret as a meaningful improvement.
The 128-bit implementation was identical in both builds; its timing variation
provides a useful reminder that small differences can be measurement noise.

A separate experiment compiled the 256-bit kernel with AVX-512F/VL enabled,
allowing additional vector registers and EVEX instructions. Runtime dispatch
retained the AVX2 fallback. This was slower on the measured host:

| K | AVX2 native encode | AVX-512VL native encode |
|---|---:|---:|
| 2^16 | 6.009 ms | 6.404 ms |
| 2^18 | 28.799 ms | 29.559 ms |
| 2^20 | 128.518 ms | 136.336 ms |

The extra backend and public backend selector were removed after this test.
The retained implementation uses AVX2 for 256-bit elements and AVX-512 for
512-bit elements. Wider ISA support alone did not improve this 256-bit kernel.

## Reproduction and scope

Host: Peach, Ryzen 9 7950X, GCC 15.2, CPU 15, compiler tuning `znver4`.
Each variant used three fresh processes, with one warmup and seven measured
trials per process. Tables report medians of process medians. Variant order
reversed on the second repeat. All benchmarks ran serially under shared locks.
These are selected operating points, not an exhaustive width or tile search.

Use the same caller-supplied upstream snapshot as the earlier forward work,
located at `ROOT/hypercat/native/spin`. From the repository root:

```sh
bash workstreams/spin_optimized/wide_schedule_check.sh ROOT
python3 workstreams/spin_optimized/wide_iteration_summary.py
bash workstreams/spin_optimized/wide_validate.sh ROOT \
  workstreams/spin_optimized/measurements/wide_dfs_safety \
  -DSPIN_FORWARD_FOUR=ON -DSPIN_FORWARD_DIRECT=ON \
  -DSPIN_FORWARD_SCHEDULE=dfs -DSPIN_WIDE_SCHEDULE=dfs
```

For an application build, enable `SPIN_BUILD_WIDE`, supply
`SPIN_BIDIRECTIONAL_SOURCE`, and set `SPIN_WIDE_SCHEDULE=dfs`.
The comparison harness also enables the current 128-bit forward optimizations.
The schedule option defaults to `global` to preserve existing build behavior.
Forward remains an optional import. The later natural-length integration adds
`MessageLength{K}` and 32-bit wide routing; the timings here precede that change.
The standalone transpose API and defaults are unaffected.

Raw measurements remain ignored under `measurements/wide_schedule` and
`measurements/wide_vl`. The latter records the rejected ISA experiment; its
temporary backend is not retained as an application option.

The release comparison passed the wide, bidirectional, and K16 bidirectional
test suites for both schedules. The retained on-demand configuration also
passed wide tests with AVX-512 detection suppressed and under ASan/UBSan.
Tests cover both selected maps, compaction, small tiles, unaligned buffers,
invalid shapes, and overlap rejection. Sanitized benchmark runs checked column
assembly at K=2^14 and the selected K=2^16 point; their timings are not reported.
The masked column benchmark correctly refused its unavailable AVX-512 path.
The release object-code audit confirmed that EVEX instructions remain confined
to the designated AVX-512 objects, including the 512-bit wide kernel.

The subsequent polish pass rebuilt from `cmake/ForwardRecommended.cmake` and
passed all five selected CTest suites and the ISA audit. The stricter scheduling
generator reproduced the measured header byte-for-byte. New regression tests
reject changed outputs, malformed XOR expressions, cycles, and missing or
duplicate stores, including under `python -O`. Wide tests now also reject
short output spans and partial overlap in either direction.
`wide_polish_check.sh ROOT` repeats these checks against the existing DFS build.

## Next step

The selected build and layout guidance are consolidated in
[FORWARD_USAGE.md](FORWARD_USAGE.md). Keep width and tile selection explicit;
the intermediate-size checks below do not justify a universal 512-bit default.

## Intermediate-size checks

A follow-up used the same host and optimized build at K=2^17 and 2^19, with
T128S19. Each setting used three processes, five measured trials per process,
and one warmup. Tile order reversed on the second repeat. Totals include column
assembly and cover the same 16-stream workload as above.

| K | Tile rows | 128-bit total | Native 256-bit total | Native 512-bit total |
|---|---:|---:|---:|---:|
| 2^17 | 256 | 18.734 ms | 16.596 ms | 17.207 ms |
| 2^17 | 512 | 18.817 ms | 16.291 ms | 18.064 ms |
| 2^19 | 256 | 82.818 ms | 74.063 ms | 74.583 ms |
| 2^19 | 512 | 84.578 ms | 75.003 ms | 76.957 ms |

At K=2^17, native 256-bit input wins. Its 512-row tile improves total time
by about 2% over 256 rows. At K=2^19, both widths prefer 256 rows, and their
totals are close. This supports 256-row tiles as a starting point through
K=2^19, rather than an early switch to 512 rows. It is not an exhaustive sweep.

At K=2^17, packing into 256-bit elements gives a total of 17.785 ms with
256-row tiles, versus 18.734 ms for the 128-bit control. Packing into 512-bit
elements loses. At K=2^19, both packed paths lose: 84.425 and 86.339 ms versus
82.818 ms for 128 bits. Input layout remains part of the performance decision.

Reproduce with `wide_range_check.sh ROOT` after `wide_schedule_check.sh ROOT`.
The script checks that its existing build uses the selected schedules and
direct routing, then runs serially without rebuilding. Results are ignored
under `measurements/wide_range`; `wide_iteration_summary.py` includes them.

The next useful step is integration into an application with native wide input,
measuring its actual layout costs. Additional micro-optimization is less urgent
than confirming that the caller can preserve the favorable input layout.
