# Wider24 packed-kernel experiments

Selected research kernel: mode4, the exact-map flattened field multiplier.
It measures 93.895 us at K=2^16 versus97.8175 us for the frozen control in the
same campaign. It retains the frozen 10%-distance /68.893749-bit certificate.
See the [research report](../../../research/workstreams/packet8_codesign/iteration6/README.md)
for proof scope and the complete performance record.

This directory does not modify or replace `../packet8_wider24`, its frozen
certificate, or the production backend. The scalar oracle and retained
implementation are compiled alongside the candidate kernels.

## Build and test

Use a Release SPIN build/package with the required x86 kernels enabled:

```sh
cmake -S spin/experiments/packet8_wider24_opt -B build/packet8-opt -Dspin_DIR=/path/to/spin/package -DCMAKE_BUILD_TYPE=Release
cmake --build build/packet8-opt -j4
ctest --test-dir build/packet8-opt --output-on-failure
```

The standalone executable is research-only. It exits77 on unsupported
hardware. CTest runs all modes at K=2048,6144,65536, with two seeds each.
The CLI is `[K] [seed] [calls] [mode] [stage-clocks]`. Mode99 with zero calls
checks every variant; use explicit mode4 for the selected implementation.

## Modes

| Modes | Purpose |
|---|---|
| 0 | Frozen control |
| 1,2,3 | Half-payload, parity-first, combined smaller-workspace layouts |
| 4,5,6 | Flat9 byte multiply, flat9 affine tables, direct16 affine blocks |
| 7 | Flat9, first symbol randomizer omitted |
| 8,9 | Parity-first flat9, without/with one omitted randomizer |
| 10–15 | Route write-prefetch distances1,2,4,8,16,32 |
| 20,21 | Structured MDS sandwich, GFNI or bitwise fixed middle |
| 99 | Test all modes, no timing |

Mode4 changes no mathematical map. Modes7/9 change the construction but
preserve its aggregate first-moment proof. Modes20/21 instantiate the
seven-diagonal MDS family and use its separate transported certificate.
Their reference Plan contains literal rows for the changed family; tests
do not compare them against the old fast outer.

The route and outer stay byte-packed between stages. `RandomizerPlan` and
`StructuredPlan` are setup-only tables. All dispatch is outside the hot
group loops. There are no hot heap allocations or runtime callback wrappers.

## Serial measurement

Never run benchmarks concurrently. The shell runner takes the shared
benchmark locks, pins one CPU, records the binary digest, and requires a
fresh output path. For example:

```sh
bash spin/experiments/packet8_wider24_opt/run_serial.sh build/packet8-opt/spin_packet8_wider24_opt /tmp/packet8-holdout271.txt 2001 15 271 0 4 7 20 20 7 4 0
```

Repeat serially for seeds419,557,863, reversing mode order on419 and863.
Use normal pages and no stage clocks for comparison with the recorded
campaign. Setup, allocation, algebra checks, and five warmups are outside
the timed region. Generated logs and binary archives remain ignored.
