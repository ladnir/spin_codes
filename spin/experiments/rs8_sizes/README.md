# One RS8 outer at K16, K18, and K20

The K20 RS8 outer now runs with the retained K16 paired-s15 inner and K18
exact-s22 inner. One unchanged outer implementation serves all three sizes.
K16 performance is effectively unchanged; K18 improves by about 1% in the
matched holdout. All three selected configurations have complete 10%-distance
certificates. This is a research port, not a public-library promotion.

## Selected configurations

Each group maps 128 message bits to 256 code bits. Eight parallel GF(16)
RS[8,4] rows form eight 32-bit columns. Independent nonzero GF(2^32)
multipliers randomize those columns. The forward outer uses their binary
adjoints; the transposed kernel uses ordinary field multiplication.
The 64 four-bit packets per group pass through the original 64-region route.
Setup is sampled once and precomputed.

| K | Inner | Precomputed transpose | Complete setup-failure margin at 10% |
|---|---|---:|---:|
| 2^16 | t64 / paired-s15, GL15 updates | 0.098599 ms | 61.8840972766 bits |
| 2^18 | t64 / exact-s22, GL22 updates | 0.461422 ms | 58.2456253319 bits |
| 2^20 | t64 / exact-s22, GL22 updates | approximately 2.87 ms | 51.8224664153 bits |

The outer and route geometry are uniform; the inner maps and state dimension
are selected by size. All inner states start at zero, persist across region
boundaries, and end without a flush. The margins concern the specified ideal
independent setup ensembles, not a measured distance for a particular seed.

The [K16 certificate](../../../research/workstreams/k16_outer_routing/RS8X8_CONSTRUCTION.md)
already covers the exact RS8 outer with paired-s15. Source and algebra checks
confirm the new flattened field circuit preserves its basis, generator,
adjoint orientation, coordinate order, and all-randomized setup distribution.
The [new K18 certificate](../../../research/workstreams/rs8_sizes/K18_CERTIFICATE.md)
refolds every saved sparse support endpoint and retains authenticated bounds
for every larger occupancy. The
[K20 certificate](../../../research/workstreams/k20_two_ms/RS8_TRANSFER.md)
is unchanged.

## Matched measurements

Times below are microseconds for single-thread, in-place transposed encoding
of 128-bit XOR records on Peach's Ryzen 7950X, core 15, GCC 15.2 Release.
The experiment uses the unchanged fixed-width, unrolled RS8 kernel from
`../k20_two_ms/rs8`; it does not replace the hot path with generic loops.
Setup, allocation, initial fill, and five warmups are outside timing.
Each call includes both inner/routing and outer work. The input is not reset
inside the timing loop. Every benchmark ran serially under the four shared locks.

| K | Outer | Median of eight process medians | Process-median range |
|---|---|---:|---:|
| 2^16 | Retained RS16 | 98.5740 | 98.383--106.138 |
| 2^16 | RS8 | 98.5985 | 97.722--99.756 |
| 2^18 | Retained RS16 | 466.3405 | 464.317--492.629 |
| 2^18 | RS8 | 461.4215 | 459.558--462.573 |

Four setup seeds (113, 211, 2027, 7919) run in both comparison orders.
Each process times 2001 calls at K16 and 1001 at K18. K16 uses normal pages;
K18 prefers explicit huge pages. Both candidates use the same policy at each
size. No outliers are removed. K16 is a tie, not a claimed speedup. The K18
gain is modest, and the earlier short screen also measured a near tie.

A separate 501-call stage screen gives RS8 inner/routing and outer medians
of 48.019 / 50.675 us at K16 and 256.569 / 203.200 us at K18.
Independently sorted stage medians need not sum to the whole-call median.

A 301-call K20 regression check measures 2.865380 ms in this harness versus
2.859969 ms in the original harness, with identical result checksums. That
single-seed check confirms consistency; the approximately 2.87 ms table entry
retains the earlier multi-campaign result, not a new public-API measurement.

## Implementation and validation

`rs8-15` and `rs8-22` compile the same existing RS8 outer sources. No outer
source was copied or changed. `Route15.cpp` exposes the exact production
paired15 route body using unchanged production helpers. Its baseline outer
also uses the production helper. The s22 path links production stages directly.
All geometry adaptation and setup are outside encoding; the hot path allocates
nothing and uses compile-time call composition.

At K16, the actual production paired15 plan is retained unchanged. Extra
research lengths adapt geometry and regenerate the same GL15 distribution
only during setup. Such lengths are implementation tests, not new certificates.
In particular, `rs8-15` at K18 and `rs8-22` at K16 are not selected certified
configurations in this report.

Validation passed:

- 32 Release CTests: four modes, K2176/K6144/K16/K18, and seeds 17/41.
  K2176 places region boundaries inside physical inner steps.
- The same 32 tests with AddressSanitizer and UndefinedBehaviorSanitizer.
  New experiment translation units and the reused RS8 sources are instrumented;
  unchanged linked production objects are not rebuilt with instrumentation.
- An additional complete K20 RS8 check in the new harness.
- Independent scalar routes and complete forward/transpose adjoint identities,
  four legal alignments, guards, scratch-padding independence, and in-place
  prefix output with the original suffix preserved.
- Eight K18 certificate tests, independent mathematical review, and a parent
  replay reproducing the exact certificate. Four retained K16 certificate
  tests also pass, including the direct physical-step q1 replay.

The new K18 receipt is `tmp/k18-rs8-gf32-s22-v1.json`, SHA256
`d9ca35896b4963c5ac4eaceb1ee5ba54ad7b7e3a9310c92518755906c9eed6f0`.
Raw evidence stays ignored under `tmp/`: `rs8-sizes-tests-v3.txt`,
`rs8-sizes-sanitized-v3.txt`, `rs8-sizes-holdout-v2.txt`, and
`rs8-sizes-stage-regression-v2.txt`.
The measured source overlay is `rs8-sizes-v2.tar.gz`, SHA256
`b809c1c416a0129b3c3373d35c43723dd78f21558e83f07b49284e8cc63fd39c`.
The later v3 overlay adds this report and K2176 test registrations only;
it does not change the measured C++ sources.

## Reproduce

Build a current SPIN package, then configure this isolated experiment:

```sh
cmake -S spin/experiments/rs8_sizes -B build-rs8-sizes \
  -Dspin_DIR=/absolute/path/to/spin-build \
  -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_COMPILER=g++-15
cmake --build build-rs8-sizes -j2
ctest --test-dir build-rs8-sizes --output-on-failure -j1
```

CLI: `mode K seed calls phase`, where mode is `baseline15`, `rs8-15`,
`baseline22`, or `rs8-22`. Zero calls performs correctness checks. Phase 0
times the full encoder; phase 1 adds stage clocks; phases 2 and 3 are
diagnostics for route and outer separately, not complete encoder timings.
The research driver accepts natural multiples of 128 from K2048 within the
retained route-offset limit. It does not extend any public sizing contract.

For measurements, acquire the shared locks and run only one process at a time:

```sh
flock /tmp/prindal-addition-encoder-benchmark.lock \
flock /tmp/bare-spin-benchmark.lock \
flock /tmp/hypercat-benchmark.lock \
flock /tmp/hypercat-global-benchmark.lock \
taskset -c 15 build-rs8-sizes/spin_rs8_sizes rs8-15 65536 113 2001 0
```

This common outer is now promoted as
[`Parameters::PacketRs8T64`](../../RS8.md), with paired-s15 through K16 and
exact s22 above. The public API preserves these research maps and matches their
timings in a fresh serial comparison. Earlier identifiers, research controls,
and certificates remain available; the default paper-SPIN profile is unchanged.
