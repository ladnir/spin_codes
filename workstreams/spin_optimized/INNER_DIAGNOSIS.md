# Remaining forward-inner cost

Replacing the S19 expansion tables with a generated XOR circuit did not improve
single-stream encoding. Keep the table implementation. The earlier BCH scheduling
and direct-routing improvements remain the best measured settings.

## Isolating routed reads

The diagnostic evaluates the same S19 inner on identical values, using either
scattered reads from packed BCH positions or an already-contiguous input buffer.
A third mode performs only the gather and writes a contiguous buffer.

At K=2^18:

| Isolated operation | Median time |
|---|---:|
| Inner with contiguous input | 0.632 ms |
| Gather only | 0.445 ms |
| Inner with routed reads | 1.007 ms |

The contiguous case still includes linear input reads and output writes; it is
not a pure arithmetic measurement. Nor are these timings additive components
of the full encoder: instructions, memory access, and cache residency interact.
They show substantial cost in both the contiguous inner and its routed reads.

The driver reconstructs the exact direct-route positions from the public tiled
view and then supplies the same values to both inner entry points. It initializes
the raw input independently of BCH encoding. The inner's control flow depends on
the fixed setup masks, not on the input values. Both outputs are compared before
warmup and after the final timed call.

Each mode runs in a separate process, with three warmup calls and 101 measured
calls. This avoids a cache artifact in the first diagnostic, which alternated
large buffers and scanned expected outputs between timed calls. The superseded
`measurements/inner_probe` results must not be used. The summary script reads only
`inner_probe_isolated` and the subsequent circuit experiment.

The K=2^20 probe uses direct routing as a diagnostic, although the selected
encoder uses tiled routing there. Its large scattered-read cost is not the
production encoder's inner time. Full-encoder measurements remain authoritative.

## Exact expansion-circuit experiment

The candidate replaces the state-to-output expansion for `Map128S19` only.
It preserves the existing feedback circuit, transvections, raw-input contribution,
routing, and all code parameters. A deterministic synthesis produces 400 XORs
for the 19-to-128 expansion, plus 128 XORs adding the raw inputs.
The generator verifies all 128 output linear forms against `Map128S19.columns`.

Full single-stream forward encoding, measured in an alternating-order A/B sweep:

| K | Existing tables | Expansion circuit |
|---|---:|---:|
| 2^18 | 1.614 ms | 1.619 ms |
| 2^20 | 9.005 ms | 9.169 ms |

The K=2^18 difference is within process variation; K=2^20 regresses by about 1.8%.
The isolated contiguous inner is also effectively unchanged: 0.632 versus
0.636 ms at K=2^18. This candidate does not support replacing the existing tables.
It does not establish that every circuit-based implementation must lose.

The complete release tests pass against dense forward and transpose references,
including both layouts, compaction, small tiles, and adjoint identities.
The rejected circuit was not taken through a new sanitizer campaign or promoted
to a supported default. `SPIN_FORWARD_EXPANSION` remains OFF, and wide lifting
is explicitly excluded for that research option. Transpose kernels are unchanged.

## Reproduction and next step

Host: Peach, Ryzen 9 7950X, GCC 15.2, CPU 15, `SPIN_TUNE=znver4`.
Results are medians of three process medians, each based on 101 measured calls.
Benchmark modes or variants reverse order on the second repeat. The shared locks
and process checks prevent concurrent benchmarks. Build work finishes first.

```sh
bash workstreams/spin_optimized/inner_probe_check.sh ROOT
bash workstreams/spin_optimized/forward_expansion_check.sh ROOT
python3 workstreams/spin_optimized/inner_summary.py
```

The full-encoder control for the second script comes from the earlier
`forward_schedule_check.sh ROOT` run. Source generation imports the frozen
`ROOT/hypercat/native/spin` snapshot without editing it. Raw results and generated
receipts remain ignored; they do not belong in the source commit.

The next useful step is to consolidate the winning forward settings into a
small, tested configuration policy: on-demand four-row BCH, direct routing for
cache-resident sizes, and tiled routing for larger sizes. Keep the remaining
K=2^18 gap documented rather than add another unproven inner setting.
