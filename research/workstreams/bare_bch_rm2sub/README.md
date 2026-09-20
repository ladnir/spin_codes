# Bare BCH-256/RM2Sub encoder sampling

Start with [PERFORMANCE.md](PERFORMANCE.md) for measurements and recommendations.
The selected (128,19) map now has full distance/setup certificates through the
five tested proof sizes K=2^16, 2^18, 2^20, 2^22, and 2^24. See the
[paper handoff](../bch_rm2sub_bridge/PAPER_HANDOFF.md). This does not extend
the implementation's accepted sizes below or add new runtime measurements.
This implementation has no parity fanout. The imported implementations remain unchanged.

Four exact selected inner maps share one setup, routing, and BCH implementation:
`t64_s20`, `t128_s19`, `t64_s16`, and `t256_s14`.
Message exponents 16 through 20 are accepted; exponents 16, 18, and 20 are tested
and measured. The public operation maps `2K` input blocks to `K` output blocks.

The optional `Outer::Bch128x32` adds a quarter-rate `4K` to `K` operation with
the same selected t128_s19 inner. See the [quarter-rate implementation and report](../rate_quarter_bch/implementation/README.md).
The default constructor and the historical rate-half measurements remain unchanged.

## Build and run

The standalone build uses the imported cryptoTools header subset. No external
numerical packages or regeneration step is needed to compile the checked-in code.

```sh
cmake -S workstreams/bare_bch_rm2sub -B out/bare-spin -G Ninja -DCMAKE_BUILD_TYPE=Release -DSPIN_ARCH=znver4
cmake --build out/bare-spin -j 3
ctest --test-dir out/bare-spin --output-on-failure
```

MSVC uses `/arch:AVX2` and ignores `SPIN_ARCH`. On other GCC hosts, supply an
appropriate architecture such as `-DSPIN_ARCH=native`. There is no runtime ISA dispatch.

### Optional workspace routing optimization

On Linux, add `-DSPIN_WORKSPACE_ROUTING_OPT=ON` to opt in to large-page advice
for owned workspace buffers and write-prefetch in the quarter-rate tile scatter.
The option is off by default. It applies only to `Outer::Bch128x32` with
K >= 2^18; smaller messages and half-rate instances retain the original path.
Non-Linux builds retain the original path even when the option is enabled.

Workspace construction requests `MADV_HUGEPAGE` and, when available,
`MADV_COLLAPSE`. Failed requests are tolerated; the existing allocation remains
usable with ordinary pages. This is a best-effort hint, not a guarantee of
large-page backing. It changes neither system settings nor caller-owned memory.
Page preparation can add workspace-construction latency and is excluded from
encoding timings. Copying or replacing workspace vectors does not repeat the hints.

The option does not change the linear map, public API, retained setup size, or
workspace size. It adds no allocation or policy branch inside the coordinate
loop. See the [deployment validation](../inner_design/routing_opt/deployment/README.md)
for correctness, forced-fallback tests, and measurements. The experimental
asymmetric inner remains a separate implementation; enabling this option does
not select it.

Run benchmarks only when no other benchmark is active. Linux runs pin CPU 15
and refuse to start if another benchmark executable is found. Do not run the
benchmark on Windows to reproduce the published numbers; its affinity guard
currently targets Linux.

```sh
# Odd trial count; tile=0 selects the measured default; layout=0 is packed24.
./out/bare-spin/spin_benchmark 101 0 0
# Explicit tile and 32-bit routing, K=2^20 only.
./out/bare-spin/spin_benchmark 31 1024 1 20
```

Each command runs its cells serially. Reported online time excludes setup and
allocations. Output is JSON Lines. `report.py RECEIPT_DIRECTORY` validates the
three final runs against correctness output and source hashes, then rebuilds
PERFORMANCE.json and PERFORMANCE.md.

## API

```cpp
#include "Spin.h"
using namespace bare_spin;
Spin code(Configuration::T128S19, 20);
code.validateSetup();                 // optional setup diagnostic
code.compact();                       // retain packed24, discard oracle storage
Spin::Workspace workspace(code);      // allocate once per active call
std::vector<block> input(code.codeBlocks()), output(code.messageBlocks());
// Populate input before calling.
code.encode(input.data(), input.size(), output.data(), output.size(), workspace);
```

`encode` checks buffers, geometry, and overlap. `encodeInplace` checks one `N`-block
buffer and replaces its first `K` blocks without an input copy; the suffix is unchanged.
`encodeUnchecked` omits checks; input/output may be disjoint or exactly equal.
The caller must provide buffers disjoint from workspace storage, a matching workspace, and
a retained routing layout. `compact(Indices32)` instead keeps the 32-bit layout;
subsequent encode calls must explicitly request that layout.

Setup can be shared by concurrent callers after construction and compaction.
Each active call requires its own workspace and input/output buffers. Do not
compact while another call is using the setup. `reference` and `validateSetup`
require un-compacted diagnostic schedules. Reinitialization means constructing
a new object. The default tile policy uses 256 rows up to exponent 18 and 2048
rows above it for the rate-half outer; an explicit power-of-two tile count of
at least two overrides this choice. Quarter-rate defaults are documented separately.

## Generated code and provenance

`generate.py` authenticates local selected-map snapshots, reconstructs their
columns, checks `BA=0`, chooses irreducible field polynomials, and generates
specialized quadratic circuits. It symbolically verifies the pruned zeta
schedule and both BCH circuit directions. The output manifest records the
exact selected maps, field polynomials, generator hash, and source dependencies.

```sh
python -B workstreams/bare_bch_rm2sub/generate.py
```

Regeneration requires the existing bridge input snapshots and the repository's
BCH and circuit-generation scripts. It is deterministic, but intentionally
rewrites generated files; inspect the diff before accepting regenerated data.
The generated source is approximately hundreds of KB, not experimental datasets.

The runtime oracle uses ordinary dense multiplication and independent polynomial
reduction. The optimized kernel retains fixed-width SIMD, compile-time dispatch,
explicit table construction, and paired outer rows. No virtual interface,
`std::function`, or per-call allocation is introduced in the hot path.
