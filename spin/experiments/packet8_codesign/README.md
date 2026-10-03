# Byte-native packet prototype

This is an isolated research encoder, not a library promotion. Its K16
precomputed transpose takes about 91 us, but it has no complete distance
certificate. The [research report](../../../research/workstreams/packet8_codesign/README.md)
records the proof screens and preserves the certified four-bit controls.

## Construction and representation

Each outer group encodes 128 bits into 256 bits using four GF16 RS[16,8]
rows and independent nonzero-transitive 16-bit symbol maps. The implementation
reuses the frozen native field16 outer setup and its literal scalar oracle.

Each group independently permutes its 32 consecutive byte packets into
32 regions. Each region independently permutes the groups. This is a new
byte route, not adjacent entries of an old four-bit route.

For eight input bytes x and state (a,b), the literal forward inner is

```
y_i = x_i + a + alpha_i*b,  alpha_i = i in AES GF256, i=0,...,7
next_state = M*(a,b) + (sum_i x_i, sum_i alpha_i*x_i).
```

Each physical 64-coordinate step has an independent sampled GL(2,GF256)
matrix M. Setup is fixed across all messages. The initial state is zero,
continues across region boundaries, and is discarded at the end.
The ideal sampling model and numerical bounds are specified in the research
report; a reproducible seed does not itself certify the realized code.

The public prototype input and output remain arrays of 128-bit XOR elements.
Each eight-coordinate packet is internally packed into two adjacent ZMM
vectors. The inner and route retain this representation. The outer consumes
it directly, avoiding a second input transpose, and unpacks the final output.
Binary adjoints are required for every GF256 scalar in the transposed inner.

`Plan` owns setup tables and addresses. Hot entry points allocate nothing.
Scratch is disjoint and 64-byte aligned; input/output need only 16-byte
alignment. Output may overwrite the input prefix after scratch materialization.
The prototype retains unused setup portions of its reference `rs::Plan`;
it is not an optimized setup-time implementation.

## Build and check

Use a Release SPIN package with the required x86 kernels enabled:

```sh
cmake -S spin/experiments/packet8_codesign -B /tmp/packet8-build \
  -DCMAKE_BUILD_TYPE=Release -Dspin_DIR=/path/to/spin/package
cmake --build /tmp/packet8-build -j8
ctest --test-dir /tmp/packet8-build --output-on-failure -j1
```

The measured compiler was GCC 15.2 on a Ryzen 7950X. The target enables
AVX512, VBMI, GFNI, and write-prefetch, tunes for Zen 4, and disables LTO.
The CMake file carries the full flag list. The new translation units also
compile with MSVC, but native execution in this study was on Peach.

Six native configurations pass: K=2048, 6144, and 65536, each with seeds
1 and 17. Every configuration checks all enabled variants, scalar/SIMD
agreement, full forward/transpose adjoints, four pointer alignments, padding
and output guards, input preservation, and in-place encoding. Independent
portable tests live in the linked research directory.

## Benchmark interface

```
spin_packet8 [K=65536] [seed=1] [calls=0] [mode=0] [phases=0]
```

Zero calls runs correctness checks only. Every timed invocation checks first,
then allocates its timing buffers and performs five untimed warmups. Setup,
allocation, and checks are outside the measured interval. Mode 3 with phases
zero is the selected headline path. Stage clocks are diagnostic only.

| Mode | Work |
|---|---|
| 0 | Whole encoder, initial outer circuit |
| 1 / 2 | Route-only paired stores / initial outer only |
| 3 / 4 | Whole encoder / outer only, frozen shared-parity circuit |
| 5 / 6 / 7 | Whole shared-outer encoder, write-prefetch 8 / 16 / 32 packets ahead |
| 8 / 9 / 10 | Route only with those prefetch distances |
| 11 / 12 | Whole encoder / outer only, half-payload shared circuit |

`run_ab.sh` holds the shared host locks throughout a serial ABBA campaign.
Never run another benchmark concurrently. For the selected path:

```sh
PACKET8_MODE=3 bash spin/experiments/packet8_codesign/run_ab.sh \
  /path/to/frozen/spin_k16_codesign /tmp/packet8-build/spin_packet8 \
  /tmp/fresh-packet8-holdout.txt 2001 15 41 113 257 997
```

## Results, 2026-10-02

The holdout used CPU 15 fixed at 4.5 GHz, normal pages, no stage clocks,
four fresh seeds, and 2,001 calls per run. Each seed ran control/candidate/
candidate/control. No sample was removed.

| Variant | Median of eight run medians | Range of run medians |
|---|---:|---:|
| Frozen certified four-bit mode 52 | 99.3405 us | 98.264--106.760 us |
| Byte-native mode 3 | 91.1200 us | 90.258--97.091 us |

The time reduction is 8.3%. The candidate is uncertified, so this is a
comparison of implementations of different constructions, not a replacement
claim. Earlier two-seed screens gave about 92.4 us for mode 0 and 91.3 us
for mode 3, against roughly 98.6--98.7 us for the control.

Stage diagnostics show about 44.5 us for inner plus route and 46 us for
the selected outer. The inner has no stack spills in the inspected GCC
assembly. The outer uses a 4 KiB intermediate and thirteen spilled ZMM
results. A same-circuit, half-payload variant cuts temporary storage to
2 KiB, but doubles output-store instructions. It loses: 101.9--102.4 us
whole-call versus 90.4--91.0 us for mode 3 in that matched screen.

Write-prefetch improves route-only time from roughly 43.1 to 40.1 us.
All three distances leave full inner/route time near 45 us and do not
improve whole encoding. Neither prefetch nor half-payload mode is selected.

A separate serialized 100,001-call `perf stat` run of mode 3 gives
90.699 us median and about 2.02 instructions per cycle. Multiplexed counters
report many cache misses but only about 1.3 dTLB load misses per call.
These whole-process counters include startup and are diagnostic, not a
per-stage attribution or proof that one particular cache level dominates.

The measured final executable SHA-256 is
`3422a6b80abb97384d38b4246d8b32a2c588d880882057aa148523b759e5ba28`.
The frozen control executable is
`5678ca1235eade2b86843ff422d450638f39d2b889ff529d465635fa89024d33`.
The local new-source archive is `tmp/packet8-first-study-source-20261002.tar.gz`,
SHA-256 `41b3ef53f76d6e49c67035e5cb40468b78d5efc922194b50d8988037c55faafa`.
Raw logs and this archive remain ignored. The archive contains the new
implementation files before this README; retained dependencies are the
unchanged frozen K16 sources, not a copy of the dirty production kernels.

Next: change and screen the construction before more low-level tuning.
The present middle-occupancy bound cannot support a complete certificate.
