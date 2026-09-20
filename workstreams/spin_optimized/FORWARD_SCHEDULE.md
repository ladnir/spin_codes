# Forward BCH scheduling

Follow-up: [INNER_DIAGNOSIS.md](INNER_DIAGNOSIS.md) separates routed-read costs
from contiguous-inner work and records a rejected expansion-circuit experiment.
The settings in this note remain the best measured choices.

Computing BCH outputs on demand improves single-stream forward encoding by
4–9% over the preceding direct-routing implementation. The circuit still uses
2,901 XORs and computes the same linear map. Output-local recomputation loses.

## Matched schedule comparison

| K | Previous schedule | On-demand schedule (`dfs`) | Time reduction |
|---|---:|---:|---:|
| 2^16 | 0.385 ms | 0.349 ms | 9.2% |
| 2^18 | 1.756 ms | 1.628 ms | 7.3% |
| 2^20 | 9.509 ms | 9.097 ms | 4.3% |

These are complete forward-encoding times for one 128-bit element stream.
They include BCH, routing, and IMT, but exclude setup and validation.
The first two sizes use direct routing; K=2^20 uses tiled routing with 256 rows.
The configurations are (t,s,rounds)=(64,12,2) at K=2^16 and (128,19,1) otherwise.
The timings are warm repeated-call latency, not cold-cache latency.

A fresh alternating-direction confirmation measured:

| K | Forward (`dfs`) | Transpose, separate output | Transpose, in place |
|---|---:|---:|---:|
| 2^16 | 0.348 ms | 0.340 ms | 0.340 ms |
| 2^18 | 1.601 ms | 1.460 ms | 1.419 ms |
| 2^20 | 9.031 ms | 9.289 ms | 9.135 ms |

Forward is now near parity at K=2^16 and K=2^20. It remains about 10% slower
than separate-output transpose at K=2^18. The default transpose tiles are
unchanged; forward uses direct routing below K=2^20 and a 256-row tile at K=2^20.

All variants ran on Peach's Ryzen 9 7950X, CPU 15, with GCC 15.2 and
`SPIN_TUNE=znver4`. Each point uses three processes, each with 101 measured calls.
The table reports the median of process medians. Variant order reverses on
the second repeat. Builds and correctness tests finish before serial timing.
Every process holds both benchmark locks and checks for another benchmark.

## What changed

The old source creates many intermediates before storing outputs. The new
generator visits each output's dependencies recursively, emits each shared
intermediate once, and stores that output before visiting the next one.
It preserves the existing XOR expressions and fixed-width AVX-512 operations.
It adds no runtime scheduler, allocation, indirect call, or per-element dispatch.

The compiled function's stack reservation falls from 37,640 to 28,296 bytes.
Static assembly lines referencing the stack fall from 4,094 to 3,808; total
instruction lines remain similar (6,669 versus 6,682). Static counts are not
retired-instruction measurements. The lower stack traffic supports, but does
not independently establish, improved intermediate-value locality as the cause.

Two alternatives isolate groups of 16 or 32 outputs in non-inlined helpers.
They shorten intermediate lifetimes but repeat shared work across helpers:

| K | 16-output groups | 32-output groups |
|---|---:|---:|
| 2^16 | 0.497 ms | 0.444 ms |
| 2^18 | 2.201 ms | 1.994 ms |
| 2^20 | 11.352 ms | 10.362 ms |

Neither grouped variant is selected. The 16-output variant uses 6,729 XORs
instead of 2,901. Smaller live sets alone do not imply better performance.

## Correctness and reproduction

Before emitting code, `forward_schedule.py` evaluates every output as a 128-bit
linear form and compares all 256 forms with the committed BCH generator matrix.
The schedule changes neither the permutation nor the inner. No new distance
certificate is needed. Release tests compare the complete encoder with dense
references and check both layouts, compaction, small tiles, and adjoint identities.
The winning schedule also passes both test suites with AVX-512 detection disabled
and under AddressSanitizer/UndefinedBehaviorSanitizer. The object-code ISA audit
passes. The standalone transpose archive is byte-identical to the earlier
single-stream build (SHA-256
`76e7a8c5b2c0a089ebd66cda44478fa7b3e41bff7bcf7017cc32a6d44031fb69`).

Use the optional forward build with:

```sh
-DSPIN_BCH_AVX512=ON -DSPIN_FORWARD_FOUR=ON \
-DSPIN_FORWARD_DIRECT=ON -DSPIN_FORWARD_SCHEDULE=dfs
```

`SPIN_BIDIRECTIONAL_SOURCE` still identifies the frozen upstream source.
`SPIN_FORWARD_SCHEDULE` defaults to `global`; neither the public default nor
the standalone transpose kernel changes in this experiment.

Reproduce from a repository root with `ROOT/hypercat/native/spin` available:

```sh
bash workstreams/spin_optimized/forward_schedule_check.sh ROOT
bash workstreams/spin_optimized/single_direct_validate.sh ROOT dfs
bash workstreams/spin_optimized/forward_schedule_confirm.sh ROOT
python3 workstreams/spin_optimized/forward_schedule_summary.py
```

Raw measurements, disassembly, and build logs stay ignored under
`measurements/forward_schedule`. Generator receipts stay in each build directory.
The preceding single-stream results and stage profiles are in
[SINGLE_FORWARD.md](SINGLE_FORWARD.md).

The remaining priority is the K=2^18 gap, followed by consolidation of the
winning forward settings. Avoid another broad scheduling sweep until the inner's
arithmetic cost has been separated from its gather cost.
