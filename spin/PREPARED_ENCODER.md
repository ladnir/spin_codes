# Prepared encoder integration checkpoint

The `PreparedEncoder` API separates full seed replacement from heuristic refresh.
See the [usage guide](README.md#repeated-transposed-encoding-and-seed-updates).
libOTe's integration selects full setup for regular noise and banked
heuristic setup for stationary noise. This API is available in SPIN 0.2.

## Validation

- GCC 13.3/Linux: API tests, original known answers, and complete K=2^18
  route/mask equality with the `row-rotate1` prototype for eight seeds.
- GCC 15.2/Peach: AVX2 and AVX-512 API checks against the independent dense oracle.
  AddressSanitizer and UndefinedBehaviorSanitizer tests also pass.
- MSVC/Windows Release: API and original known-answer tests pass.
- Allocation counters verify allocation-free bank refresh and transposed encoding.
  Tests cover all three inner configurations, minimum lengths, and non-power-of-two lengths.
- libOTe: 34 protocol configurations, three rounds each, both choice layouts,
  regular/stationary noise, malicious-mode consistency, and generic coefficient contexts.

These checks establish implementation consistency, not a distance certificate
for the heuristic permutation family.

## Initial generalized performance checkpoint

Peach, GCC 15.2, Release, `SPIN_TUNE=generic`, automatic backend (AVX-512 BCH),
128-bit records, CPU 15. Three serial passes alternate mode order. Each process
uses 20 warmups and 101 measured rounds; the table gives medians of pass medians.
Inputs are initialized outside the timed region. Bank preparation and workspace
allocation are excluded; banked combined time includes fresh seed refresh each round.
All runs held the three shared encoder/Hypercat benchmark locks.

| K | Fixed full encode (ms) | Bank refresh (ms) | Bank encode (ms) | Bank combined (ms) | Combined overhead |
|---:|---:|---:|---:|---:|---:|
| 65,536 | 0.359 | 0.0281 | 0.484 | 0.513 | 43% |
| 81,920 | 0.440 | 0.0108 | 0.565 | 0.575 | 31% |
| 262,144 | 1.496 | 0.0299 | 1.810 | 1.840 | 23% |

This initial generalized implementation did **not** retain the prototype's 10% target.
The family is unchanged at K=2^18, but the execution implementation differs:
it stores natural destinations instead of the prototype's redundant indexed layout,
uses general region lengths, and batches mask generation through portable word lanes.
Those are candidate causes, not an isolated attribution of the timing difference.
The prototype also used host-specific tuning and additional ISA instructions.

At this checkpoint, K=2^18 used 2,140,280 setup bytes and 8,396,800 workspace bytes.
The region-address workspace is retained across rounds.

Build the checkpoint benchmark with `SPIN_BUILD_BENCHMARKS=ON`, then run
`spin_prepared_bench full 262144` and `spin_prepared_bench bank 262144` serially.
The executable reports setup, workspace creation, refresh, encode, and combined time.
It does not acquire benchmark locks itself; coordinate and serialize runs externally.

## K=2^18 fast-path restoration

The optimized prepared path retains the same bank-family revision, seeds, and
encoded map. It adds three execution optimizations:

- An indexed destination table and packed row keys for the fixed K=2^18 geometry.
  The AVX-512 backend uses a fixed-width epoch loop and a private aligned region buffer.
- Batched mask generation on CPUs with AVX-512 F/VL/BW/DQ and VPOPCNTDQ.
  Rejected words trigger the scalar sampler, preserving the exact stream.
- Linux huge-page advice during bank and workspace preparation, using the existing
  owned-page helper. Unsupported or unsuccessful advice leaves ordinary memory usable.

These changes neither materialize a fresh full route nor allocate during refresh
or encoding. AVX2 and other geometries retain the generalized routing path.
The public API and descriptors are unchanged.

In the first five-pass confirmation, the default `SPIN_TUNE=generic` build measured:

| K | Fixed full encode (ms) | Bank refresh (ms) | Bank encode (ms) | Bank combined (ms) | Combined overhead |
|---:|---:|---:|---:|---:|---:|
| 262,144 | 1.518204 | 0.008466 | 1.656462 | 1.664938 | 9.67% |

Measurements use the same host, CPU, record size, and per-process sampling method
as the initial checkpoint. Five serial passes alternate mode order. Preparation
and workspace allocation remain outside the timed region; every banked round
refreshes its seed. Each column is the median of its five pass medians.
The host is a Ryzen 9 7950X. A second five-pass confirmation measured 1.660840 ms
combined against 1.538271 ms fixed, or 7.97% overhead. Together these runs support
an approximately 8--10% result on this host, not a hard latency bound.

The staged measurements placed the indexed AVX2 implementation near 19% overhead,
then the AVX-512 implementation near 16%. Batched masks reduced refresh from about
31 microseconds to 8.5 microseconds. The private route buffer and huge-page advice
brought combined time near the 10% target. These sequential measurements identify
useful changes, but are not an additive decomposition of runtime.

K=2^18 now retains 4,245,672 setup bytes and 8,396,800 workspace bytes.
The indexed table accounts for most of the setup increase. The specialized kernel
also uses an 8 KiB stack route buffer, which `workspace.bytes()` does not include.
Bank construction and workspace preparation must still occur off the critical path.

This result is specific to the measured K=2^18 configuration and host. The initial
smaller-length measurements above do not establish a 10% target for those lengths.

The restored path passed GCC 15.2 Release, an AVX2-only build, and ASan/UBSan checks.
MSVC Release passed API, known-answer, and mask-stream tests. The API regression
now includes K=2^18 banked encoding against the independent dense oracle.
The mask regression compares 6,144 complete streams; Peach exercised 4,070 SIMD
successes and 2,074 scalar fallbacks. libOTe's `SpinIntegration` and `SpinOptions`
tests also pass against this local library snapshot.

The libOTe integration now measures the stationary path including fresh leaves,
refresh, compression, and hashing: about 84.5 million hashed OTs/s at K=2^18.
The regular-noise path reuses full setup and reaches 89.0 million OTs/s.
Both exclude setup, generation of supplied base correlations, and transport;
see libOTe's `libOTe/Tools/Spin/REGULAR_PERFORMANCE.md` for the timing contract.
