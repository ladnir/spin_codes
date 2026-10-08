# Direct blockwise SPIN: forward cost checkpoints

At K=2^20 and rate 1/2, the latest direct-SPIN cost control takes about
6.35 ms for one precomputed forward encoding of 128-bit records. It uses
GRS[256,128] and a width16/state4 inner. It has no full distance certificate
and does not meet the sub-5-ms target. The first width64/state16 control
remains frozen at about9.0ms; its original receipt follows before the newer
experiments below.
The [accumulator checkpoint](ACCUMULATOR_DETOUR_CHECKPOINT.md), certified
high-rate winner, and production library are unchanged.

## Exact construction

The [restart](HD_BLOCKWISE_RESTART.md) now has a native implementation:
GRS[128,64] over embedded AES GF256, one independently sampled uniform row
permutation in each coordinate region, and a direct t64/s16 stateful inner.
Each 128-bit record carries eight GF65536 lanes sharing the same scalar code.
The outer receives four interleaved rows per SIMD tile; this is an input
layout convention, not additional randomization.

The [fixed maps](HD_BLOCKWISE_MAPS.md) interpolate from nodes128..143 to
nodes0..63, with C=A^T. Every packet computes

    Y = X + A Q;
    Q_next = Q + u*(v^T Q) + C X.

Setup samples u uniformly from the nonzero state vectors and samples v
uniformly subject to v^T*u!=1. The benchmark uses a reproducible PRNG for
this sampler and Fisher--Yates routing. Its measured seed is not itself
a distance certificate. State starts at zero, carries across regions,
and is not flushed. There is no input diagonal, pre-accumulator, extra
expansion, puncturing, or uncharged routing pass.

Native sources are
[hd_blockwise_forward.cpp](../../../spin/experiments/char2_field_spin/hd_blockwise_forward.cpp)
and [HdBlockwiseInner.h](../../../spin/experiments/char2_field_spin/HdBlockwiseInner.h).
The outer reuses the tuned sixteen-register leaves in HdFastOuter.h.
The inner uses four explicitly scheduled sixteen-record tiles, with no
heap allocation or full-packet temporary in its hot loop.

## Serial K20 measurements

Host: Peach, AMD Ryzen9 7950X, GCC15.2, one thread pinned to CPU2.
The governor reported userspace, its maximum frequency4500000kHz, and
boost0. Builds use Release, AVX512, GFNI, and `-mtune=znver4`.
Payload allocation requests transparent huge pages under the existing
benchmark policy. All measurements hold the shared benchmark lock.

Setup, allocation, input production, and hashing are excluded. Each timed
full call includes outer evaluation, routed input loads, both fixed inner
maps, rank-one refreshes, and output stores. Input is regenerated before
each sample. Three warmups precede the recorded samples.

| Full forward control | Samples | Median ms | p10--p90 ms |
| --- | ---: | ---: | ---: |
| Initial insert-load implementation, no prefetch | 41 | 10.291 | 10.231--10.384 |
| Initial one-packet lookahead prefetch | 41 | 9.668 | 9.609--9.745 |
| Hardware gather plus prefetch control | 41 | 10.515 | 10.450--10.589 |
| Final source, insert loads, no prefetch | 51 | 10.874 | 10.828--10.973 |
| Final source, insert loads, one-packet prefetch | 51 | 8.958 | 8.915--9.015 |

A preceding final-source run measured8.977ms. The final phase run measured
8.964ms total, with medians3.056ms outer and5.902ms routed inner. Medians
of phases need not sum exactly to the median total. The final source keeps
one-packet prefetch as the default and leaves hardware gather disabled.
Both switches remain compile-time experimental controls. Do not compare
small cross-revision differences as though only the switch changed; the
final no-prefetch/prefetch pair uses the same source and compiler flags.

Component controls explain the current gap:

- Exact pre-routed packet inputs, same mixers:3.047ms for the inner alone.
  Routing is materialized outside that diagnostic timer, so this is not
  an encoding time. Its output checksum matches the full encoder.
- Route-only copying, no inner:4.667ms in the initial control.
- Identity-route raw-inner diagnostic:2.979ms; replaced by the exact
  pre-routed control above to avoid changing the packet inputs.

These timings are not additive. They indicate substantial routing/load
cost and interaction with arithmetic, while fixed maps still cost several
milliseconds. Assembly also contains spilled payload temporaries and
stack-resident constants; no isolated timing has assigned their share.
Lower variable-multiplication counts alone did not produce the target speed.

The retained FLOCK RS baseline is16.622ms forward and has deterministic
distance near50%. This9-ms research control is not yet a certified
replacement for that baseline.

## Correctness and reproduction

Release and ASAN/UBSAN checks pass. They cover every fixed-map basis vector,
the exact identity CA=0, independent dense scalar packet and full encoders,
nonzero incoming states, identity refresh, both uint16/uint32 route types,
raw and full exact aliasing, sixteen-byte offset variants, and buffer guards.
Full scalar comparisons use K=4096,8192,12288, including a non-power-of-two
research size. The pre-routed arithmetic diagnostic is also checked against
the same complete scalar result. The hardware-gather variant passes the
native correctness checks; it is not selected for speed.

The final forward and inner source SHA256 values are respectively

    74456421345f58f49c40b49a94494e8073ca97e7a7972e15e0beb4660458de42
    b1b1665b67c4ba1da20c98afdc1615373f6e5b6f624728dd00d78619f8e1605d

These hashes match the tested remote files. The row-major A matrix has the
independently checked digest recorded in HD_BLOCKWISE_MAPS.md. All190
high-distance Python helper tests pass after the new local/global bounds.

From the experiment directory, a Linux build can be reproduced with

```text
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_COMPILER=g++-15
cmake --build build --target spin_char2_hd_blockwise_forward -j4
build/spin_char2_hd_blockwise_forward --check
flock -n /tmp/spin-char2-benchmark.lock taskset -c 2 build/spin_char2_hd_blockwise_forward --bench 51 full 20
```

The benchmark's optional modes are `phases`, `outer`, `inner`,
`plain-inner`, and `route-only`. `plain-inner` excludes its required
materialization, and `route-only` omits all coding. Only `full` and
`phases` measure the complete encoder. Define `BW_PREFETCH_PACKETS=0`
to disable lookahead; define `HDBW_HW_GATHER=1` for the rejected gather
control. No benchmark was run concurrently with another benchmark.

## Next decision

Combine the [local weighted-state bound](HD_BLOCKWISE_MAPS.md) with the
[capped routing decomposition](HD_BLOCKWISE_GLOBAL.md). Partial-row
shortening and general output supports remain open. Keep the measured
candidate as a control while testing a cheaper fixed-map/state geometry
and better routed-load overlap. Any state reduction must repeat the
geometry and local-bound checks; the current16-symbol state is a search
choice, not a proved requirement. No production promotion, commit, or push
accompanies this checkpoint.

## Small-state cost controls and their rejection

The isolated [small-state driver](../../../spin/experiments/char2_field_spin/hd_blockwise_small_forward.cpp)
adds GRS[256,128] and width64/state4 or8 controls. It does not change the
frozen width64/state16 sources above. All rows remain independently routed
within each coordinate region. The fixed maps use interpolation from
nodes128..128+s-1 to nodes0..63, with C=A^T; the mixer distribution is unchanged.

The following runs used the same host and timing policy, with41 samples.
Padding adds four unused records between scratch regions; it changes only
the memory layout, not the encoded map. Cached and non-temporal outer
stores were measured with the subsequent routed inner included.

| GRS[256,128] control | Measurement | Total ms | Outer ms | Routed inner ms |
| --- | --- | ---: | ---: | ---: |
| s8, unpadded cached stores | phases | 10.736 | 5.938 | 4.786 |
| s8, padded cached stores | phases | 8.376 | 3.323 | 5.052 |
| s8, padded non-temporal stores | phases | 8.548 | 2.289 | 6.265 |
| s4, padded cached stores | phases | 7.827 | 3.304 | 4.518 |
| s4, padded cached stores | full | 7.829 | -- | -- |
| s4, four-packet prefetch | full | 7.809 | -- | -- |

Phase medians need not sum to the median total. The s4 full-call p10--p90
intervals are7.697--7.929ms and7.695--7.906ms, respectively. That difference
does not justify selecting the longer prefetch distance. Padding removes
much of the observed outer slowdown. Non-temporal stores improve the
isolated outer phase but lose overall; the timing does not isolate which
cache or translation effects cause the difference.

These controls are **rejected for40% distance**. The subsequently discovered
[quiet-subcode theorem](HD_BLOCKWISE_QUIET_SUBCODE.md) bounds their actual
relative distance by12.5% for s4 and25% for s8, for every route and mixer.
The larger outer and strong one-active-row bounds cannot fix this many-row
obstruction. Retain the timings to price individual components, not as a
new high-distance performance frontier.

Release checks passed independent scalar encoding at K8192,16384,24576,
all map bases, CA=0, nonzero states, raw/full exact aliasing, both route
index widths, offsets, guards, and exact pre-routed diagnostics. The
small-state controls had not received sanitizer validation at this timing
checkpoint. Their s8 and s4 full-output checksums were respectively
2442595558768910388 and6931791287664064460; layout/store/prefetch variants
preserved the checksum within each geometry.

## Smaller packets and the selected cost control

The next controls shrink both dimensions: GRS[128,64] with t32/s8,
and GRS[256,128] with t16/s4. Their quiet-subcode ceiling is50%, so that
theorem does not reject40%. Both also escape the balanced paired-window
lower bound. These are necessary screens, not certificates. The
[composition bound](HD_BLOCKWISE_COMPOSITION.md) certifies only the complete
one-active-row class, at143 and223 bits respectively.

The recurrence, independent regional routing, and rank-one mixer are
unchanged. For width t, A interpolates the s state values at128..128+s-1
to0..t-1; C=A^T. The identity CA=0 still holds. The smaller kernels use
explicitly unrolled packed transforms:54 GFNI products per t32/s8 packet
and22 per t16/s4 packet, including refresh. Both retain0.5 variable field
products per output record. Counts alone do not predict the full time.

All following full calls include outer evaluation, scratch padding, routing,
inner arithmetic, output stores, and any listed staging copy. One reusable
region buffer is allocated with scratch before timing. Staging copies each
logical region sequentially into that buffer, then executes the fused routed
inner there. State remains continuous across all regions. The padded outer
stores can use non-temporal writes; their fence is inside the full call.

| Complete encoder | Samples | Median ms | p10--p90 ms |
| --- | ---: | ---: | ---: |
| Outer128, t32/s8, cached stores, no staging | 51 | 8.494 | 8.361--8.618 |
| Outer256, t16/s4, cached stores, no staging | 51 | 7.666 | 7.545--7.746 |
| Outer256, t16/s4, cached stores, staged region | 51 | 7.358 | 7.270--7.449 |
| Outer256, t16/s4, non-temporal stores, staged region | 51 | 6.348 | 6.312--6.371 |
| Same selected control, final CMake build | 101 | 6.352 | 6.335--6.386 |
| Staged region, separate gather then raw in-place inner | 51 | 7.407 | 7.384--7.497 |

The selected phase run gives6.343ms total,2.301ms outer, and4.043ms for
staging plus routed inner. Its exact pre-routed inner diagnostic takes
2.631ms, excluding the required materialization. That number is not an
encoder time and is not subtracted from another run to claim an isolated
routing cost. Separate gather followed by a raw in-place inner is slower;
the selected path keeps routed loads fused with packet arithmetic.

The selected driver defaults to outer256, t16/s4, four-record scratch
padding, non-temporal outer stores, staged regions, and four-packet
prefetch. Changing these memory-layout controls preserves the code. Changing
t, s, or outer parameters changes the code and requires its own proof.
The production library and the frozen first native control are untouched.
Follow-up staged runs with no lookahead and one-packet lookahead measured
6.320ms and6.365ms, respectively, with51 samples each. Their checksums and
scalar checks agree. The small differences are not treated as a new speed
result; the longer101-sample receipt above remains the selected checkpoint.

Release and ASAN/UBSAN checks pass on the selected control. They cover the
same scalar, alias, guard, route-width, and nonzero-state checks described
above, now at K8192,16384,24576. All219 high-distance helper tests pass.
The full-output checksum is14884408593612809608. Source hashes, verified
between the local and tested remote copies, are:

| File | SHA256 |
| --- | --- |
| `hd_blockwise_small_forward.cpp` | `a4d1103346d9eb531969bea2b98a18c45d1cb83541cc18693c78d75cb490bedc` |
| `HdBlockwiseSmallInner.h` | `e684496b271219a516943241b967b06dd28afd6bfa9808c1afc5c16931050cd6` |
| `HdBlockwiseOuter256.h` | `22d8163076590ad8d4e2c1e4c2b31ff23ad89ddb6233d3851c660cacf74335eb` |

The final CMake executable has SHA256
`17ee6017b74120e772c20d4d32b082db1840356b1c4ada0c4287cdc8a18753ea`.
From the experiment directory:

```text
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_COMPILER=g++-15
cmake --build build --target spin_char2_hd_blockwise_small_forward -j4
build/spin_char2_hd_blockwise_small_forward --check
flock -n /tmp/spin-char2-benchmark.lock taskset -c 2 build/spin_char2_hd_blockwise_small_forward --bench 101 full 20
```

The main remaining proof task is the many-row class. Do not substitute the
singleton margin for a whole-code margin, or apply width64/state16 routing
receipts to the new geometry. The new [interval-span proof](HD_BLOCKWISE_INTERVAL_SPANS.md)
bounds its dimension-forced complete-row, one-window-per-region trap family
below 2^-528.39, but not partially shortened rows or general output histories.
Keep this6.35ms implementation fixed while closing its proof or identifying
a concrete obstruction.
