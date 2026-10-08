# Native check of the shared beta64 circuit

The first shared-map native experiment does not beat the frozen direct-SPIN
reference. Reducing GFNI products from 82 to 59 increased full encoding time.
A whole-vector 63-product rewrite recovers some of that loss, but remains
slower than the FFT implementation of the same map. Stop tuning this circuit
as the route to 5 ms; the next change should remove work structurally.

These are research controls with no complete distance certificate. None is
promoted into the library.

## Identical map, different schedules

All three shared controls implement the shift64 construction in
[the shared-map note](HD_BLOCKWISE_SHARED_MAPS.md): GRS[256,128], one regional
route, eight FFT8 tiles per 64-record packet, four global-basis state
coefficients, and a genuine GF65536 rank-one state update. The fixed maps
are over GF256. Eight GF65536 lanes share each 128-bit record.

The [cube factorization](HD_BLOCKWISE_SHARED_FACTOR.md) computes local
feedthroughs and shares their Boolean moments with feedback. Its first native
schedule uses 59 packed GFNI arithmetic calls per packet. State emission and
feedback require irregular record packing. The next schedule uses 63 calls
but keeps those two operations in whole ZMM vectors; its emission uses 26 ZMM
XORs instead of 120 individual-record XORs and repacking. Both use a bounded
aligned 1-KiB emission scratch per packet and no per-packet allocation.

The FFT control uses 82 calls. All three include routing, required region
staging, the same outer, and output stores in their full-call timings. Setup
and allocation are excluded. An arithmetic-call count omits shuffles, XORs,
loads, stores, coefficient handling, and register pressure.

## Same-host measurements

Measured serially on Peach's Ryzen 7950X, GCC 15.2, pinned to core 2 at a
4.5-GHz userspace setting with boost disabled. K=2^20 input records, N=2^21
output records, forward direction, precomputed setup. Default layout uses
four padding rows, non-temporal outer stores, charged regional staging, and
four-packet prefetch. A shared benchmark lock excluded concurrent benchmarks.

Full-call measurements use 101 samples:

| Implementation | Median, ms | p10--p90, ms |
| --- | ---: | ---: |
| Frozen identity-feedthrough t16/s4 reference | 6.345797 | 6.319878--6.368519 |
| Shared beta64 FFT, 82 GFNI calls | 6.541121 | 6.500816--6.753146 |
| Shared beta64 cube, 59 GFNI calls | 7.650852 | 7.615306--7.746822 |
| Shared beta64 whole-vector cube, 63 GFNI calls | 7.003654 | 6.976163--7.031106 |

The shared controls all give checksum 5184726652255349648. The old reference
implements a different map and gives checksum 14884408593612809608.

Separate 51-sample phase measurements give:

| Shared schedule | Full phase call | Outer | Staged route plus inner | Pre-routed inner only |
| --- | ---: | ---: | ---: | ---: |
| FFT | 6.514072 | 2.282271 | 4.232021 | 2.331122 |
| Cube | 7.631347 | 2.295776 | 5.330992 | 3.504793 |
| Whole-vector cube | 6.987304 | 2.276169 | 4.708951 | 2.967159 |

All entries are milliseconds. The last column is a diagnostic that excludes
required input materialization; it is not a full encoder result. Medians from
separate measurements need not add. Subtracting them does not isolate an
exact routing cost. The experiments show that lower product counts alone do
not predict latency; they do not supply hardware-counter attribution.

## Validation and reproduction

The research target is `spin_char2_hd_blockwise_shared_forward`, from
`spin/experiments/char2_field_spin/hd_blockwise_shared_forward.cpp`.
It accepts:

    --check
    --bench 101 full 20 fft
    --bench 101 full 20 cube
    --bench 51 phases 20 cube
    --bench 51 plain-inner 20 cube

The default cube selects `HdBlockwiseSharedInner.h`. Compile the same driver
with `-DHDBSH_NEXT=1` to select `HdBlockwiseSharedInnerNext.h` instead. This
compile-time switch changes the schedule, not the map. The FFT choice remains
available in both builds. The next control was built separately with:

    g++-15 -std=c++20 -O3 -DNDEBUG -DHDBSH_NEXT=1 \
      -mavx2 -mavx512f -mavx512bw -mavx512vl -mavx512dq -mgfni \
      -mtune=znver4 -fno-lto hd_blockwise_shared_forward.cpp \
      build/libspin_char2_field.a -o build/hd-shared-next

Both schedules pass all 68 packet basis cases in both extension components,
nonzero-state and dense independent scalar checks, exact alias/offset/guard
checks, 16-bit and 32-bit routes, and full outer integration at K=8192, 16384,
and 24576. Both also pass ASAN/UBSAN checks with the driver built at -O1.
An independent symbolic audit checked the whole-vector emission and feedback,
including arbitrary values in unused moment lanes. GCC required qualification
of the new feedback call to avoid argument-dependent lookup of the old one.

Checkpoint SHA256 values:

- Shared header: `58d2c9001135b50eacc30a26c2f8b49a5cfc268cbac16e6b05e520b28e4acfac`.
- Whole-vector header: `646d8bf139ab35a11bcdeb4f217cdb6b7a78eeb31b78c28e193886fb316d02a5`.
- Driver with compile-time selection: `6806a98d3adee9b4cd97552390b32d12a39fe1d0a683f75946100b794ba9bb68`.

The frozen t16/s4 reference, outer, accumulator source, and accumulator
archive hashes remain unchanged. No production code, commit, or push is
part of this experiment.

## Structural next step

The current feedthrough residual has rank four, but together with its
feedback map requires six independent local quantities. A new state basis
may let the same four quantities drive both operations. That is a construction
change: it must retain strong emission/check distances and be checked against
quiet-subcode and cancellation obstructions before native implementation.

Meanwhile the proof has two new reusable pieces: an
[all-density high-type quiet bound](HD_BLOCKWISE_SHARED_QUIET_COMPOSE.md) and
a [whole-joint-image mixer transfer](HD_BLOCKWISE_SHARED_TRANSFER.md).
The low-type quiet complement and general nonquiet histories remain open.
