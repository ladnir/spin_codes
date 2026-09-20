# Selected K16 inner and K20 wide-tile tuning

Follow-up: [WIDE_ITERATION.md](WIDE_ITERATION.md) records improved BCH scheduling
and compares wide kernels with the subsequently optimized 128-bit path.

The selected K16 inner now supports 256-bit and 512-bit forward encoding.
The 256-bit path gives a 1.31x throughput gain over its 128-bit counterpart,
including column assembly. Separately, reducing the K20 tile from 2048 to
512 outer rows improves the 256-bit path by 1.17x and reduces scratch space.

Both changes preserve the binary generator. K16 retains (t,s,r)=(64,12,2),
including the full 49.3276-bit setup-failure bound at 10% relative distance.
K20 retains (128,19,1). These are forward-encoder results, not replacements
for the earlier transposed-encoder measurements or full PCS timings.

## K16: the selected inner also benefits from batching

Each trial encodes 16 independent 128-bit planes under one shared setup,
then assembles identical column-major output. This is 2048 binary messages.
Widths 128, 256, and 512 require 16, eight, and four encoder calls respectively.
There is no input bit transpose. Row-to-column assembly is included throughout.

The table gives milliseconds for this fixed workload at K=2^16 and a
256-row tile. Native input already has the kernel's interleaved lane layout.
Planar input requires timed packing before each wide call.

| Width and input | Encode | Column assembly | Total | Gain over 128-bit |
|---|---:|---:|---:|---:|
| 128-bit | 8.565 | 1.530 | 10.107 | 1.00x |
| 256-bit, native | 6.194 | 1.495 | 7.692 | 1.31x |
| 256-bit, planar | 6.196 | 1.510 | 8.513 | 1.19x |
| 512-bit, native | 6.731 | 1.708 | 8.420 | 1.20x |
| 512-bit, planar | 6.209 | 1.692 | 8.865 | 1.14x |

The total includes packing where required; phase medians need not sum to it.
The native 256-bit run medians span 7.661–7.713 ms across three processes.
For comparison, the S19 inner takes 7.797 ms at the same width and tile.
The main improvement is batching, not a large additional speedup from the
K16 parameter change. The selected K16 parameters retain the stronger margin.

The generated wide recurrence applies T1, then T2, then one feedback addition.
It emits output before the update, starts with zero state, and does not flush.
These are the same operations as the existing 128-bit forward implementation.
The fixed A/B circuits are lifted from `Map64S12.h`; width lifting performs
no new circuit search.
See [the K16 proof result](../inner_design/k16_parameters_20260919/FULL49_RESULT.md).

## K20: smaller tiles help, but smallest is not best

At 256-bit width, each outer row occupies 256 times 32 bytes in the tile.
The old 2048-row tile therefore occupies 16 MiB. A 512-row tile occupies 4 MiB.
The global routing buffer remains 64 MiB, so workspace decreases from 80 to
68 MiB. Setup, input, output, and column-assembly buffers are additional.

The initial serial screen tested 8, 16, 32, 64, 128, 256, 512, 1024, 2048,
and 4096 rows, in forward and reverse order. Native 256-bit total times ranged
from about 159 ms at 512 rows to 200 ms at eight rows and 214 ms at 4096 rows.
Shrinking the tile alone is not a monotone optimization. The tile also changes
the bucket layout and subsequent gather addresses; this is not an isolated
test of one cache-size threshold.

We confirmed 256, 512, and 2048 rows across all three widths. Each entry is
the pooled median of 21 samples, in milliseconds, for the same 16-plane workload:

| Tile rows | 128-bit | 256-bit native | 256-bit planar | 512-bit native | 512-bit planar |
|---:|---:|---:|---:|---:|---:|
| 2048, old default | 201.110 | 185.945 | 203.265 | 200.069 | 219.244 |
| 256 | 191.603 | 162.435 | 183.856 | 177.922 | 199.993 |
| 512 | 194.319 | 158.685 | 180.478 | 171.179 | 193.313 |

For native 256-bit input, the 512-row choice is 1.17x faster than the old
256-bit default. It is 1.21x faster than the best 128-bit tile in this
confirmation set, and 1.27x faster than the old 128-bit default.
For planar input, packing reduces the gain over the tested 128-bit control
to 1.06x. Prefer native wide input when the application can supply it.

At 512 rows, native 256-bit run medians span 154.15–159.26 ms. At 256 rows,
they span 161.80–163.38 ms. The winner remains 512 rows despite timing variation.
This confirms a useful operating point, not global optimality over every tile
size, seed, CPU, or message length. We did not sweep all tile sizes for 128-bit.

## Use and validation

Build with `SPIN_BUILD_WIDE=ON` as described in [WIDE_FORWARD.md](WIDE_FORWARD.md).
Keep the existing transposed defaults unchanged. For a dedicated K20 forward
setup, select the measured tile explicitly:

```cpp
bare_spin::Spin code(bare_spin::Configuration::T128S19, 20, 1, 2, 512);
code.compact();
bare_spin::WideWorkspace workspace(code, 2);
workspace.forward(input, output); // spans of interleaved 128-bit lanes
```

For K16, use `Configuration::T64S12R2` and exponent 16 with the default tile.
The setup must outlive its workspace and remain unmoved and unmodified.
The old `wideView()` contract stays S19-only; the new `wideForwardView()` exposes
the configuration tag and schedules for the two supported wide forward maps.
Map and ISA selection occur outside the hot loops. No allocations were added
inside an encode, and no runtime round-count loop was introduced.

The release tests compare wide K16 outputs against the dense forward oracle.
They compare different tile schedules against an independent default-tile
128-bit encoder, including K20 at 64 and 512 rows. They also cover compacted
setup, deliberately non-vector-aligned buffers, bad shapes, and overlap rejection.
Every timed workload checks its final column output against independent scalar
assembly. Object-code inspection confines AVX-512 to the designated objects.
All eight release tests passed. The expanded wide tests also passed with
AVX-512 masked off and under AddressSanitizer/UndefinedBehaviorSanitizer.
Sanitizer runs cover both S19 and K16 column-assembly workloads.

The host and workload are unchanged: Ryzen 9 7950X, CPU 15, GCC 15.2,
`-O3 -mtune=znver4`, routing seed 1, coefficient seed 2. Setup, allocation,
page initialization, and output checking are excluded. All benchmarks are serial
and acquire both shared locks. K16 and confirmation cases use three processes,
each with seven measured trials after one warmup; case and tile order vary.

Replay in the Peach snapshot layout:

```sh
bash workstreams/spin_optimized/wide_tune.sh /path/to/snapshot
bash workstreams/spin_optimized/wide_confirm.sh /path/to/snapshot
bash workstreams/spin_optimized/wide_validate.sh /path/to/snapshot \
  workstreams/spin_optimized/measurements/wide_tuning
python3 workstreams/spin_optimized/wide_tuning_summary.py
```

Raw samples, validation logs, source hashes, and binary hashes are under the
ignored `measurements/wide_tuning` directory. The current defaults are unchanged;
the selected K20 tile is explicit until an application-level dispatch policy
is measured. Next test these choices inside the full Brakedown commitment path,
where hashing and data retention may change the benefit.
