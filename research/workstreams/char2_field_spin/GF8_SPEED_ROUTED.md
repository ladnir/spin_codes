# Exact-map optimization of the 83.33% GF256 code

The [subsequent phase audit](GF8_SPEED_PROFILE.md) measured the separate
inner, route, and outer costs and tested sixteen further schedules. It found
no material improvement, so this remains the selected checkpoint.

The selected K20 implementation now takes about **3.2 ms**, versus about
4.0 ms for the previous cached-eight implementation in matched tests.
Restoring the original input before every timed call gives 3.33--3.36 ms,
versus 4.13--4.19 ms. These are precomputed transposed-encoding times for
128-bit records, each containing sixteen parallel GF256 lanes.

The [construction and certificate](GF8_SPEED_DESIGN.md) are unchanged:
GRS[64,40], four input and state coordinates, three retained outputs per
step, 10% symbol distance, and 57-bit whole-code setup-failure margin at
K20. Nominal rate is 83.3333%; requested rate is 83.2520% after padding.
The requested dimension is 1,048,576 and retained length is 1,259,520.
Production profiles and the [previous checkpoint](GF8_SPEED_NATIVE.md)
remain unchanged.

## Selected schedule

Four changes preserve the exact sampled linear map:

1. Apply the input diagonal during inverse routing. Setup reorders its
   coefficients so that `routedDiagonal[route[i]] = inputDiagonal[i]`.
   This removes diagonal work from the recurrence loop.
2. Process sixteen packets with a staggered schedule. Only two four-packet
   tiles are live at once. Completed tiles write directly to the regional
   temporary; their registers then prepare earlier tiles. The recurrence
   still runs in the original reverse order across all region boundaries.
3. Fuse the outer's right leaves with their ancestor butterflies and retained
   stores. Hot scratch shrinks from 64 to 48 SIMD vectors. The arithmetic
   uses GFNI multiplication, not the slower affine alternative measured here.
4. Use a scratch-column stride of `groups+4` records. Since `groups` is a
   multiple of 64, this stride is an odd number of cache lines. It avoids
   giving every column the same page offset. Final output stores remain
   non-temporal when alignment permits them.

The research entry point is `rs4TransposeRouted16` in
`spin/experiments/char2_field_spin/RateSpeed4RouteNext.h`.
Its `RateSpeed4RoutedPlan` adds one reordered coefficient byte per raw
record. It does not allocate the gauge or pair-lookahead experiment tables.
The hot path allocates nothing. Plan construction and scratch allocation
remain outside the timed encoding operation.

Input and output may alias exactly. Scratch must be disjoint from both.
The entire inner/routing phase completes before output stores begin.
The implementation retains natural-size shortening, unaligned output,
and the 32-bit route fallback. These implementation checks do not extend
the certificate to new dimensions.

## Final entry-point measurements

The final driver compares only the previous cached-eight entry point and
the selected routed-sixteen entry point. Both use the same original payload,
input address, scratch addresses, and sampled code. It allocates only the
selected plan, without unused gauge or pair tables. Candidate order rotates
between batches of five. Each launch measures 151 calls per candidate.

| Protocol / launch | Previous median | Selected median | Selected p90 |
|---|---:|---:|---:|
| Restore before each five-call batch, 1 | 3.997 ms | 3.184 ms | 3.354 ms |
| Restore before each five-call batch, 2 | 4.003 ms | 3.201 ms | 3.408 ms |
| Restore before every call, 1 | 4.134 ms | 3.332 ms | 3.392 ms |
| Restore before every call, 2 | 4.194 ms | 3.361 ms | 3.443 ms |

The reduction is about 20% under both protocols. Payload restoration,
setup, allocation, and checksums are excluded; all encoding and output
stores are included. Restoring the payload changes cache pressure. It does
not establish cold-memory latency or application throughput.

A separate 29-candidate phased screen measured 2.428 ms for selected
inner/routing and 0.820 ms for the outer. Its full median was 3.246 ms;
the matched old schedule took 4.036 ms. These phase values describe that
screen, not a decomposition of the later two-candidate timings.

Measurements use Peach's Ryzen 9 7950X, GCC 15.2, AVX512/GFNI, CPU 2,
and disabled boost. Benchmark allocations of at least 16 MiB request
huge pages. Every benchmark ran serially under the same host lock.

## What did not win

The screen isolated layout, arithmetic, and their combinations. Representative
matched values from the final broad screen are below; all values are medians.

| Schedule | Full time | Restore every call |
|---|---:|---:|
| Previous cached-eight | 4.036 ms | 4.116 ms |
| Hardware gather | 4.516 ms | 4.537 ms |
| Relocated diagonal only | 3.832 ms | 3.906 ms |
| Fused outer only | 3.852 ms | 3.962 ms |
| Relocated diagonal + fused outer + padded stride | 3.389 ms | 3.495 ms |
| Same, with sixteen-packet pipeline | **3.246 ms** | **3.375 ms** |
| Same, with gauge-normalized recurrence | 3.458 ms | 3.550 ms |
| Same, with precomputed pair lookahead | 3.409 ms | 3.498 ms |
| Sixteen-packet pipeline + batched route coefficients | 3.323 ms | 3.479 ms |

Gauge normalization moves the running state scales into setup. The carried
state then uses an XOR recurrence, but normalization adds independent field
multiplications and coefficient loads. Pair lookahead halves the chain of
dependent state multiplications at a smaller arithmetic cost. Both are
correct and improve some isolated schedules; neither improves the selected
complete schedule. Their implementations remain comparison controls.

Routing prefetch was neutral or worse. Batching four routing diagonals into
one load also lost in the selected schedule. Hardware gather was consistently
slower than four scalar 128-bit loads. Affine outer multiplication and outer
prefetch did not beat the multiply-fused outer.

Scratch-stride effects varied across launches. Therefore the recommendation
uses the measured complete schedule and matched controls, not an assumed
sum of independent percentage improvements or a claimed hardware-counter
explanation.

## Validation and reproduction

Release checks pass for all 29 schedules. Small cases compare against the
independent dense-Cauchy/direct-polynomial oracle. Large cases compare the
sampled linear map with the previous implementation at K16, K18, and K20.
Checks cover both outer dimensions 39 and 40, all 128-bit alignment offsets,
shortened tails of one through three records, the coefficient-32 boundary,
input preservation, exact aliases, unchanged alias suffixes, and guards.
Every candidate also runs with the compact route removed, forcing the
32-bit fallback.

An additional 256-trial local test supplies arbitrary nonzero terminal dual
states. Gauge, pair lookahead, and the sixteen-packet pipeline match the
previous kernel's outputs and final state. Gauge state is converted back
using the setup routine's returned initial scale. Independent review checked
the reordered diagonal, normalization algebra, route bounds, and unchanged
sampling assumptions.

The final ASAN/UBSAN replay passes both the independent small checks and
the full K16/K18/K20 checks, including all 29 alternatives and the selected
entry point. Leak detection is enabled.

```text
cmake -S spin/experiments/char2_field_spin -B build/field -DCMAKE_BUILD_TYPE=Release
cmake --build build/field --target spin_char2_rate_speed4_next -j2
build/field/spin_char2_rate_speed4_next --check
build/field/spin_char2_rate_speed4_next --large-check
taskset -c 2 build/field/spin_char2_rate_speed4_next 1048576 151 selected
taskset -c 2 build/field/spin_char2_rate_speed4_next 1048576 151 selected-fresh
```

The `full`, `fresh`, and `phases` modes retain the complete 29-candidate screen.
Never run benchmark commands concurrently.

Next, profile inner arithmetic separately from routing, using this complete
schedule as the control. Inner/routing still accounts for roughly three
quarters of full time. If reducing that work requires a construction change,
retain the 80--84% rate target and close a new 10% distance certificate.
