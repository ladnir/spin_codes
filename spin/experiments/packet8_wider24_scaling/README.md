# Byte-packet scaling experiment

**Parked research checkpoint, 2026-10-03.** Keep the certified four-bit code
as the main design. This experiment is retained for reproduction, not promotion.

This isolated experiment measures the frozen wider24 construction at K16,
K18, and K20. It does not change the outer, routing distribution, inner maps,
or state size. The retained flattened nine-product field multiplier stays
packed through the inner and outer. The new alternatives change only stores.
No production kernel is promoted here.

The serial holdout gives about **0.094 ms at K16**, **0.373 ms at K18**, and
**3.389 ms at K20**, for precomputed transposed encoding of 128-bit elements.
K16 uses ordinary pages and cached stores; K18 uses requested huge pages and
cached stores; K20 uses requested huge pages and streaming route/output stores.
Setup and allocation are excluded. The matched older four-bit K20 control is
**3.268 ms**, so the byte construction is 3.7% slower there.

See the [performance record](../../../research/workstreams/packet8_codesign/iteration7/PERFORMANCE.md)
for every holdout run median, policy screens, and source identities.
Only K16 has a complete distance certificate. The
[larger-size proof screen](../../../research/workstreams/packet8_codesign/iteration7/README.md)
does not establish a whole-code failure probability at K18 or K20.

## Build and check

```sh
cmake -S spin/experiments/packet8_wider24_scaling -B build/packet8-scaling \
  -Dspin_DIR=/absolute/path/to/spin-build/spin -DCMAKE_BUILD_TYPE=Release
cmake --build build/packet8-scaling -j4
ctest --test-dir build/packet8-scaling --output-on-failure
```

All ten native tests pass: K=2048, 6144, 65536, 262144, and 1048576, each
with seeds 1 and 17. Each tests all six modes against the literal scalar
encoder, 128 adjoint lanes, four pointer alignments, in-place operation,
unchanged input/suffix, guards, and route padding. Non-temporal output has
a cached fallback for pointers that are not 64-byte aligned. The route's
aligned scratch is internal. Each streaming stage fences before returning.

The executable exits 77 on unsupported hardware. CLI:

```text
[K=65536] [seed=1] [calls=0] [mode=1] [normal|huge] [phases=0]
```

| Mode | Route stores | Outer/output |
|---|---|---|
| 0 | Cached | Frozen original randomizer |
| 1 | Cached | Retained flattened randomizer, cached output |
| 2 | Streaming | Retained flattened randomizer, cached output |
| 3 | Cached | Copied flattened circuit, streaming output |
| 4 | Streaming | Copied flattened circuit, streaming output |
| 5 | Cached | Copied flattened circuit, cached output |
| 99 | All modes | Correctness only; calls must be zero |

Modes 1 and 5 have the same arithmetic; the latter controls for copying the
outer into a translation unit with store-policy specialization. Their noisy
small timing differences are not a new construction improvement.

## Serial measurements

Never run two benchmarks concurrently. The runners acquire all four shared
host benchmark locks, pin one CPU, record binary hashes, and require fresh
output paths. Both include five untimed warmups. Stage clocks must remain
off for headline timings; when enabled, route/outer phases are measured in
the same full call, not in separate cache-warmed loops.

```sh
bash spin/experiments/packet8_wider24_scaling/run_serial.sh \
  build/packet8-scaling/spin_packet8_wider24_scaling /tmp/byte-k20-screen.txt \
  1048576 41 201 15 huge 0 1 2 3 4 5
bash spin/experiments/packet8_wider24_scaling/compare_k20.sh \
  build/packet8-scaling/spin_packet8_wider24_scaling \
  /absolute/path/to/spin_rs16x8_border /tmp/byte-k20-holdout.txt 301 15
```

For K16/K18, run modes `1 5 5 1` for seeds 271 and 557, and `5 1 1 5`
for seeds 419 and 863, each with 1001 calls, separately for normal/huge
page policies. Do not select from individual fastest calls.

`anon_huge_kib` is process-wide Linux `smaps_rollup` residency after warmup,
not per-buffer evidence. The memory field records the requested policy.
No claim is made that every intended allocation received huge pages.
Logs and generated receipts remain ignored.
