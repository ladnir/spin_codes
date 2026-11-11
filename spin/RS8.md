# RS8 packet profile

`Parameters::PacketRs8T64` selects a precomputed half-rate binary code, with
stable descriptor identifier 8. It is optimized for transposed encoding of
128-bit records. The public interface is `Code`, `Workspace`, and the full
mode of `PreparedEncoder`.

```cpp
spin::Code code({.message_size=1u<<18,
                 .parameters=spin::Parameters::PacketRs8T64,
                 .route_seed=17, .inner_seed=41});
auto work=code.make_workspace();
auto buffer=code.make_buffer();
// Populate buffer with 2K records before encoding.
code.transpose_inplace_bytes(buffer.bytes(),work);
```

Setup and allocation occur before encoding. The in-place call overwrites the
first K records and preserves the remaining K. Separate-buffer transpose is
also supported. Forward encoding supports 128/256/512-bit records through a
scalar fallback; this promotion does not provide a tuned forward kernel.
Generic XOR records, bit-packed forward encoding, and banked heuristic setup
are not supported by this profile.

## Construction and sizing

Each group maps 128 input bits to 256 output bits using eight RS[8,4] rows
over GF(16). Eight independently randomized 32-bit columns connect this outer
to the four-bit packet route. The reverse column maps are multiplication by
nonzero GF(2^32) elements; forward encoding applies their binary adjoints.

The outer, routing structure and inner are shared across sizes: the inner is
the exact s22 map at every K (until 2026-11-11, K <= 65,536 used paired-s15,
which was only numerically certified at K16; seeded maps at those sizes changed).
The identifier therefore determines the inner without another caller option. `route_seed` determines the route and outer maps in separate
sampling domains; `inner_seed` determines the inner updates.

Accepted message sizes are positive multiples of 128 up to 132,152,832.
Call `valid_message_size` before allocating. `tile_rows` accepts 0 or 8.
The AVX-512/VBMI/GFNI kernel is selected at runtime, with a baseline SSE2
fallback. Owned buffers use the existing automatic large-page policy.

## Certificate scope

The completed bounds below give 10% relative distance under the corresponding
random-setup model. Margin lambda bounds the probability of sampling a bad
code by 2^(-lambda). The implementation preserves the checked research maps
and setup domains. Accepting another natural size does not certify it.

| Message size | Inner | Distance lower bound | Margin (bits) |
|---|---|---:|---:|
| 2^14 | Exact s22 | 3,277 | 64.25152 |
| 2^16 | Exact s22 | 13,108 | 62.34042 |
| 2^18 | Exact s22 | 52,429 | 60.36266 |
| 2^20 | Exact s22 | 209,716 | 58.36785 |
| 2^22 | Exact s22 | 838,861 | 56.3691 |
| 2^24 | Exact s22 | 3,355,444 | 54.36978 |
| 2^26 | Exact s22 | 13,421,773 | 52.32161 |

Every margin is a Lean theorem (`Spin.Turbo.G{e}Final.g{e}_endpoint` in
`research/lean`; regeneration inputs in `research/lean/ladder_plans/`), using
only Lean's standard axioms.

Earlier profiles remain available with their original identifiers and seeded
maps. In particular, RS8 does not inherit the larger margin of the RS16 s22
profile: its separate certificate accounts for the changed outer.

## Validation

The 2026-10-03 promotion comparison measured these precomputed, single-threaded
128-bit in-place transpose times through the public API:

| K | Public API | Frozen research | Public/research difference |
|---|---:|---:|---:|
| 2^16 | 0.099526 ms | 0.099050 ms | +0.48% |
| 2^18 | 0.461416 ms | 0.460890 ms | +0.11% |
| 2^20 | 2.861814 ms | 2.869859 ms | -0.28% |

Each entry is the median of eight process medians: four seeds, each tested in
both comparison orders. K16/K18/K20 used 2001/1001/301 timed calls per process,
after five warmups. Setup, allocation, and input filling were excluded; inputs
were not reset between calls. The machine was a Ryzen 7950X, pinned to core 15,
with GCC 15.2 and `SPIN_TUNE=znver4`. All benchmarks ran serially.

K16 used normal memory. K18 and K20 explicitly requested `PreferHugePages`
for data and scratch, matching the research campaign. These are memory-policy
settings, not a guarantee that the operating system supplies huge pages.
The default `Automatic` policy requests them only for allocations of at least
16 MiB; callers reproducing the K18 timing should request them explicitly.
All matched checksums agree. No process was excluded; one K18 public process
had a 0.486107 ms median, versus 0.460119--0.461882 ms for the other seven.
The small aggregate differences show no material promotion regression.

`spin_rs8_test` exercises public API behavior, portable/optimized agreement,
natural-size boundaries, forward adjoints at all widths, in-place suffix
preservation, workspace reuse, and full prepared encoding.
`SPIN_BUILD_RESEARCH_ORACLE_TESTS=ON` adds byte-for-byte comparisons with the
frozen RS8 research implementation. `SPIN_BUILD_BENCHMARKS=ON` builds the
public API timing driver. Neither option is needed by consumers.

Promotion validation passed all 12 Linux Release CTests and all nine Windows
Release CTests. Installed-package consumers passed on Linux and Windows.
The RS8 API suite also passed with AVX-512 disabled. A full-library Linux
AddressSanitizer/UndefinedBehaviorSanitizer rebuild passed both RS8 test targets,
including the optional frozen oracle. Fixed forward/transpose known answers
are checked in ordinary builds. All eight K18 certificate tests and four
retained K16 certificate tests passed. Independent source review found no blocker.

The optional comparison driver is built with both switches enabled:

```sh
cmake -S spin -B build-rs8 -DCMAKE_BUILD_TYPE=Release -DSPIN_TUNE=znver4 \
  -DSPIN_BUILD_TESTS=ON -DSPIN_BUILD_RESEARCH_ORACLE_TESTS=ON \
  -DSPIN_BUILD_BENCHMARKS=ON
cmake --build build-rs8 -j2
ctest --test-dir build-rs8 --output-on-failure -j1
# Run these sequentially under the same affinity and with no competing benchmark.
taskset -c 15 build-rs8/spin_rs8_bench public 262144 113 1001
taskset -c 15 build-rs8/spin_rs8_bench frozen 262144 113 1001
```

Local raw evidence remains ignored under `tmp/`. The measured source archive
`rs8-promote-v2.tar.gz` has SHA256
`315f266fd11fbd108e9b843ecd9a174112b55d9594295333cde30281dca4bb33`.
The timing log `rs8-promotion-holdout.txt` has SHA256
`43ddb7606521e0b6d6214fbda1e1a5f407cbb6e5e203c81a9d1416e1826100a9`.
The later test-only overlay `rs8-promote-tests-v3.tar.gz` adds fixed known
answers without changing kernels; its SHA256 is
`b228200e022233dc16390752f1f6898335b9263d4f4d63de868e4bd7fba943ac`.
