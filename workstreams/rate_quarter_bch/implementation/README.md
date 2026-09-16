# Quarter-rate BCH/RM2Sub implementation

Start with [PERFORMANCE.md](PERFORMANCE.md) for measured results and certificate scope.
The [follow-up optimization notes](OPTIMIZATION_NOTES.md) record rejected variants
and the proposed next step: a jointly optimized, equivalent state representation.
This implements the fixed [128,32,32] outer with the selected RM2Sub `(t,s)=(128,19)` inner.
The complete transposed operation maps `4K` input blocks to `K` output blocks.
Each 128-bit block carries 128 independent binary instances.

The implementation extends the [shared bare encoder](../../bare_bch_rm2sub/README.md).
It reuses the optimized inner unchanged. The additional outer circuit is generated here;
the routing, workspace, checked interfaces, and benchmark executable are shared.
The old constructor still selects the rate-half outer by default.

## Build and verify

Run from the repository root. Generated C++ and the required vendor headers are checked in.
No Python package is needed for the C++ build.

```sh
cmake -S workstreams/bare_bch_rm2sub -B out/quarter-spin -DCMAKE_BUILD_TYPE=Release -DSPIN_ARCH=znver4
cmake --build out/quarter-spin -j2
ctest --test-dir out/quarter-spin --output-on-failure -j1
```

Use an appropriate `SPIN_ARCH` on other GCC hosts. MSVC uses `/arch:AVX2`.
The implementation has no runtime ISA dispatch. The published measurements use Linux, not MSVC.

## In-place API

```cpp
#include "Spin.h"
using namespace bare_spin;
Spin code(Configuration::T128S19, 20, 1, 2, 0, Outer::Bch128x32);
code.validateSetup();                  // optional diagnostic
code.compact();                        // retain only packed24 routing
Spin::Workspace workspace(code);       // allocated once
std::vector<block> buffer(code.codeBlocks());
// Populate all 4K input blocks.
code.encodeInplace(buffer.data(), buffer.size(), workspace);
// The first K blocks now contain the result; the suffix is unchanged.
```

The two seeds select the route and nonzero field coefficients. Setup uses the
same row/region permutation structure as the certificate. State continues across regions.
All input is consumed into workspace before output begins, so the in-place call needs no input copy.
The caller's buffer must be disjoint from the workspace buffers.

For separate buffers, use the existing checked `encode` interface; it rejects overlap.
`encodeUnchecked` accepts separate buffers or exact input/output pointer equality,
but omits geometry and retained-layout checks. Use the checked interface unless those invariants are already established.
After setup and compaction, a code object can be shared; concurrent calls need separate workspaces and buffers.

Message exponents 16 through 20 are accepted. Exponents 16, 18, and 20 were tested and measured.
Quarter-rate default tile counts are 256 for exponent 16, 2,048 for exponents 17–18,
and 4,096 for exponents 19–20. Only the three measured sizes informed this policy.
An explicit power-of-two tile count of at least two overrides the default, capped at the number of outer rows.
Quarter-rate inner choices other than t128_s19 are rejected.

## Reproduce timing

Do not run two benchmarks at once. The Linux executable takes the shared benchmark
lock, checks for other visible benchmark processes, and pins itself to CPU 15.

```sh
# trials, tile (0=default), layout (0=packed24), exponent (0=all), quarter, config (2=t128s19), inplace
for run in 1 2 3; do
  ./out/quarter-spin/spin_benchmark 101 0 0 0 1 2 1
  ./out/quarter-spin/spin_benchmark 101 0 0 20 0 2 1
done
```

All cells run sequentially. The second command supplies the same-session rate-half reference.
Setup, allocation, and input initialization are excluded; in-place timing performs no reset or copy between calls.
The timer covers `encodeUnchecked` with equal input/output pointers. The public `encodeInplace` interface adds argument validation.
For a tile/layout sweep, supply the explicit tile count and layout `0` or `1`.

Sanitizer check:

```sh
cmake -S workstreams/bare_bch_rm2sub -B out/quarter-spin-asan -DCMAKE_BUILD_TYPE=Debug \
  -DCMAKE_CXX_FLAGS="-O1 -g1 -fsanitize=address,undefined -fno-omit-frame-pointer"
cmake --build out/quarter-spin-asan -j2
ASAN_OPTIONS=detect_leaks=1 UBSAN_OPTIONS=halt_on_error=1 \
  ctest --test-dir out/quarter-spin-asan --output-on-failure -j1
```

## Regeneration and proof identity

```sh
python workstreams/rate_quarter_bch/implementation/generate.py
python -m unittest discover -s workstreams/rate_quarter_bch/implementation -p 'test_*.py' -v
python -m unittest discover -s workstreams/rate_quarter_bch -p 'test_*.py'
python workstreams/rate_quarter_bch/implementation/report.py
```

The generator reconstructs the certified outer, row-reduces its generator without
changing its span, synthesizes the transpose, and checks every output symbolically.
The [manifest](generated/MANIFEST.json) records the basis, selected inner, source hashes, and generated circuit hashes.
It does not regenerate the inner or modify any certificate-bound proof source.
The proof tests use the numerical dependencies listed in the parent workstream.

The [outward certificates](../SMALLER_OUTWARD_CERTIFICATE.md) apply at `K=2^20`:
16.5% relative distance with setup failure below `2^-40`, or 19% with failure below `2^-30`.
These are two interpretations of the same configuration, not different timed encoders.
