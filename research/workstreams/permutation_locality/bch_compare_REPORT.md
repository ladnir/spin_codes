# Fresh BCH Implementation Comparison

2026-09-30. GFNI wins this comparison on 128-bit XOR elements. At `K=2^20`,
its isolated BCH phase takes 22.1% less time than the nonaliasing XOR circuit;
the controlled full encoder takes 15.7% less time. A bounded polynomial-based
prototype loses to both. These results compare concrete implementations,
not the best possible algebraic BCH algorithm.

Production sources and defaults are unchanged. The retained research encoder
remeasures at approximately 5.83--5.84 ms. The new comparison harness
has a roughly 0.39 ms inner/routing overhead discussed below; its full
timings must not replace that retained optimized timing.

## Three Code-Equivalent Implementations

All three compute the same binary BCH[256,128] transpose on four rows of
128-bit elements. Input is coordinate-major across those four rows; output
is row-major. One call consumes 16 KiB and produces 8 KiB.

- `bchTranspose4`: production's generated AVX-512 XOR circuit.
- `bchTranspose4Restrict`: the existing generated nonaliasing XOR variant,
  with the same output forms and a different evaluation schedule.
- `bchTranspose4GfniBlend`: the existing GFNI tile implementation, including
  input bit transposes, prepared-input stores, all GFNI arithmetic, output
  transposes, and systematic/parity additions.

All objects were rebuilt with the same GCC 15.2 optimization and ISA flags.
The headers for the inner, selected maps, setup sampling, and SIMD block
match the current checkout byte-for-byte. The matrix header matches too.
The comparison therefore does not reuse an older optimized XOR object
against a freshly compiled GFNI object.

GFNI does not change the alphabet to GF(2^128). It rearranges bits locally
so one hardware instruction evaluates an 8-by-8 binary matrix in each byte.
The 128 bitplanes of each element remain independent. The XOR circuit and
GFNI tiles implement the same dense binary coefficients in different ways.

## Matched Results

Each entry below is the median of six process medians: seeds 1 and 17,
three rotated execution orders per seed, and 101 timed calls per process.
Each process has three warmups. Runs are serial and hold all three shared
benchmark locks. No compiler work runs concurrently with the measurements.

### Repeated Cache-Hot BCH Calls

Each timed batch evaluates 2,048 four-row tiles, reusing a pool of 1, 4, or
16 tiles. Values are nanoseconds per four-row tile. “Cache-hot” does not
mean that the complete source, destination, and workspace fit in L1.

| Reused tile pool | Production XOR | Nonaliasing XOR | GFNI |
|---|---:|---:|---:|
| 1 tile | 1185.376 ns | 1115.290 ns | 803.064 ns |
| 4 tiles | 1210.349 ns | 1156.191 ns | 803.823 ns |
| 16 tiles | 1214.437 ns | 1165.823 ns | 804.772 ns |

For one tile, the corresponding median invariant-TSC counts are 5324.183,
5009.381, and 3606.965. These are TSC ticks, not measured core cycles.
The TSC/clock ratio is approximately 4.49 ticks/ns on this host.

### Complete Bulk BCH Phase

The input tiles use the same 1,028-element stride as the current shared
packet route. Output tiles are contiguous. There is no route or inner
work in this table. Times include each BCH kernel's complete preparation.

| Message length | Production XOR | Nonaliasing XOR | GFNI |
|---|---:|---:|---:|
| `K=2^16` | 0.160223 ms | 0.151862 ms | 0.106087 ms |
| `K=2^18` | 0.639613 ms | 0.610239 ms | 0.424920 ms |
| `K=2^20` | 3.711607 ms | 3.462278 ms | 2.697653 ms |
| `K=2^22` | 17.012992 ms | 16.512159 ms | 11.608602 ms |

At `K=2^20`, process medians range over 3.676741--3.732656 ms,
3.438960--3.478824 ms, and 2.687604--2.704475 ms, respectively.

The growth between `2^18` and `2^20` is not an arithmetic change. The
source and destination footprint grows from about 12 MiB to 48 MiB.
This is consistent with cache-capacity effects, but no hardware-counter
experiment attributes the slowdown to a particular cache or memory event.

### Complete Shared-Packet Encoder

The construction, setup seeds, GF16 packet maps, and IMT parameters remain
fixed. Only the BCH implementation changes. All sizes use IMT(128,19) with
two updates for this controlled comparison; these are not separately tuned
per-size configurations or new distance certificates.

| Message length | Production XOR | Nonaliasing XOR | GFNI |
|---|---:|---:|---:|
| `K=2^16` | 0.471612 ms | 0.459323 ms | 0.405718 ms |
| `K=2^18` | 1.856488 ms | 1.819408 ms | 1.537758 ms |
| `K=2^20` | 7.632447 ms | 7.432078 ms | 6.267352 ms |
| `K=2^22` | 31.737313 ms | 30.942628 ms | 26.074956 ms |

At `K=2^20`, process medians range over 7.604254--7.913793 ms,
7.427023--7.889918 ms, and 6.236259--7.152821 ms. Outlying processes remain
in the retained records and reported ranges; they were not discarded.

Separate instrumented runs use 31 calls plus warmups. At `K=2^20`, the
medians of their per-process phase means are:

| Phase | Production XOR | Nonaliasing XOR | GFNI |
|---|---:|---:|---:|
| Inner and route | 3.572626 ms | 3.580987 ms | 3.561212 ms |
| BCH | 3.992511 ms | 3.861375 ms | 2.661422 ms |

These diagnostic means are not an exact decomposition of the independent
101-call headline medians. In particular, one GFNI process had a larger
inner/routing mean; the median phase statistic avoids silently replacing
the headline measurement with that instrumented outlier.

## Why the Fresh Full Timing Is Not 5.86 ms

A direct control calls the original `joint.cpp` encode template inside the
fresh executable, rather than the comparison harness's equivalent wrapper.
It gives about 6.22 ms, matching the wrapper. The retained original binary
still gives 5.842809 and 5.830406 ms for seeds 1 and 17 on the same host.

The difference lies in inner/routing: approximately 3.15 ms in the retained
driver versus 3.54--3.56 ms in the fresh driver. BCH remains approximately
2.68 ms in both. The discrepancy is therefore not evidence that GFNI or
the BCH map became slower. This bounded investigation does not identify
the cause. No synthetic “corrected” full
timing is presented, and the optimized retained timing remains unchanged.

## Bounded Algebraic Screen

The separate `bch_algebraic.py` generator constructs a fixed-coefficient
polynomial multiplication circuit with Karatsuba leaves of size 32, then
transposes its XOR graph. Its exact variant includes the basis conversion
and rank-five correction needed to reproduce the production generator.
Its raw variant uses a different message basis spanning the same code.

This is one generated prototype, not a tuned polynomial-encoder campaign.
The following table uses a separate, smaller screen: one 31-call process
per seed, with seeds 1 and 17. Each displayed value is the midpoint of
the two process medians. The three established kernels were rerun in the
same screen; they are not imported from the 101-call table above.

| Implementation | Cache-hot tile | Bulk BCH, `K=2^20` | Full encoder, `K=2^20` |
|---|---:|---:|---:|
| Production XOR | 1189.409 ns | 3.694565 ms | 7.545135 ms |
| Nonaliasing XOR | 1128.475 ns | 3.451106 ms | 7.431458 ms |
| GFNI | 804.693 ns | 2.694226 ms | 6.205814 ms |
| Algebraic, exact production map | 1760.148 ns | 4.724814 ms | 8.903406 ms |
| Algebraic, alternate message basis | 1306.835 ns | 4.124672 ms | Not measured |

Exact algebraic bulk results are 4.697533 and 4.752095 ms; raw-basis bulk
results are 4.138442 and 4.110901 ms. Exact full results are 8.949116 and
8.857695 ms. The raw variant is not a drop-in replacement for callers
expecting the production message basis, so its full encoder was not timed.

The generator formally verifies every output form. It also verifies that
the raw and production row spaces agree and both have rank 128. Thus the
slower exact result is not explained by comparing unrelated codes.

## Compiled Footprint

GCC's stack-usage reports and object symbol sizes provide a concrete reason
not to equate fewer source-level XORs with faster execution.

| Selected kernel path | Text bytes | Compiler-reported stack frame |
|---|---:|---:|
| Production XOR helper plus wrapper | 65,811 | 18,888 B helper; 8 B wrapper |
| Nonaliasing XOR | 60,945 | 13,000 B |
| GFNI entry plus eight selected tile helpers | 18,098 | 8,320 B entry; 8 B each helper |
| Algebraic, exact | 104,244 | 31,880 B |
| Algebraic, alternate basis | 75,112 | 25,288 B |

GFNI also uses a 2 KiB coefficient table. Its entry frame includes the
deliberate 8 KiB prepared-input array. The non-inlined tile helpers retain
their accumulators in registers. Static disassembly contains 2,594
stack-addressed vector-move instructions in the production helper and
1,246 in the nonaliasing XOR helper. These are static instruction counts,
not measured dynamic spill traffic or a hardware roofline model.

The algebraic generator emits 5,867 vector XOR operations for the exact
map and 4,327 for the alternate basis. Those are source-graph counts, not
CPU instruction counts. The larger code and frames are consistent with
the measured losses, without proving that no better algebraic schedule
or algorithm exists.

## Validation and Reproduction

All 131,072 physical input basis vectors of the four-row 128-bit-element
kernel were checked. The nonaliasing and GFNI variants match production.
The exact algebraic variant also matches production; the raw variant
matches its independent dense polynomial-basis reference. Three additional
dense inputs check that reference. All tested operations are fixed binary
linear maps, so the complete basis comparisons establish map equality.

Complete shared-R2 encoder and in-place suffix checks pass for all three
primary backends: two seeds and three input patterns at `K=2^14`, plus
three patterns at `K=2^20`. The exact algebraic full encoder passes the
three-pattern `K=2^20` check too. Every timed process checks its result
before timing. No new sanitizer campaign was run for this research comparison
harness; it reuses the previously checked routing and encoder kernels.

Machine: Peach Ryzen 7950X, pinned core 15, GCC 15.2. Common flags:

```text
-O3 -DNDEBUG -std=gnu++20 -mavx2 -mavx512f -mavx512vl -mavx512bw
-mgfni -mtune=znver4 -fno-strict-aliasing
```

The three BCH objects all receive those flags. The driver does not itself
require `-mgfni` because the GFNI intrinsics are in the separate object.
Both XOR kernels therefore use the same permitted ISA as GFNI; the fresh
production-source object is not claimed to be bit-identical to a shipped
library object.

Setup, allocations, input initialization, reference checks, and checksums
are outside timing. The encoder is in-place and is not reset between
calls. Isolated BCH inputs are read-only; each call writes the complete
output. Both source and output are disjoint for the nonaliasing kernel.

Files: `bch_compare.cpp`, `bch_compare_run.sh`, `bch_algebraic_perf.cpp`,
and `bch_algebraic.py`. The existing `schedule_bch.py` and `gfni_bch.py`
generate the two established alternatives. On Peach:

```text
root=/tmp/spin-bch-compare-2e2UPS
reference=/tmp/spin-joint-9n57TT
bash "$root/bch_compare_run.sh" "$root" build "$reference"
bash "$root/bch_compare_run.sh" "$root" check "$reference"
bash "$root/bch_compare_run.sh" "$root" full "$reference"
bash "$root/bch_compare_run.sh" "$root" hot "$reference"
bash "$root/bch_compare_run.sh" "$root" bulk "$reference"
bash "$root/bch_compare_run.sh" "$root" profile "$reference"
bash "$root/bch_compare_run.sh" "$root" control "$reference"
bash "$root/bch_compare_run.sh" "$root" legacy-control "$reference"
bash "$root/bch_compare_run.sh" "$root" build-algebraic "$reference"
bash "$root/bch_compare_run.sh" "$root" algebraic "$reference"
```

The build directory contains copies of `joint.cpp`, its two packet headers,
the three generators, and current `BchAvx512.cpp` / `BchCircuit.h` sources.
The third argument supplies the unchanged SPIN headers. All measured
objects are freshly generated and compiled in the comparison directory.
Raw logs, `.su` files, disassembly, and symbol sizes are retained remotely
there and locally under ignored `tmp/bch-compare-results/`.

Object SHA256 values:

```text
production: bdcaca623b62a9b23215e471887f951078459fe4cc52dc1ea3e3c3a100d49d4c
nonaliasing: 3f45cc50d56283d4ce66b37dbecccdea5da5f7e333ab83e54826d5b93b8f8fba
GFNI: 6caa3115bd985d056cc91892f0387b0a77582accb6c2f1bdddd0c85753d1ac5d
algebraic: 9eaeb1b52943449eddd728397972d5ad35759155277562633ba5ca13dec54696
```

The initial full-run executable hash is
`5b0279376ae610c345625b3a4cfccbe2275f2bf455b1a082915154794b8547da`.
Adding the original-driver control produced
`c2ec0c6c5fa39cabd696d16bc1e4585c7d55a2873a8d8c0c0d29fd3c5c27e255`,
used for hot, bulk, phase, and control runs with the same BCH objects.
The separate algebraic executable hash is
`e1081983df19865ba4707a832e7fe0ceadd1d4ffcf9a02689a5c39b20445ac48`.

Recommendation: retain GFNI for the packet research path. Investigate a
new algebraic algorithm only with a specific reason to beat these measured
baselines; source-level sparsity alone is insufficient.
