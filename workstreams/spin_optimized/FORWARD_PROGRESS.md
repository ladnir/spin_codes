# Forward encoding: batching and gather experiment

For current wide-element timings after BCH rescheduling, see
[WIDE_ITERATION.md](WIDE_ITERATION.md). The wide tables below are historical.

For matched **single-stream** timing and the subsequent direct-routing experiment,
see [SINGLE_FORWARD.md](SINGLE_FORWARD.md). The measurements below use 16 streams.

Four-row BCH batching improves the 128-bit forward path by 4.5–6.8% in this sweep.
The code, routing permutation, and IMT parameters are unchanged. This is an
implementation change, so it requires no new distance certificate.
The standalone transposed library is unchanged.

## Measured improvement

The workload encodes 16 independent 128-bit element streams. The table divides
encoding time by 16; it excludes setup, packing, and row-to-column assembly.
These are amortized times within that workload, not isolated-call latency.

| K | IMT (t,s,rounds) | Original two-row forward | Four-row forward | Time reduction |
|---|---|---:|---:|---:|
| 2^16 | (64,12,2) | 0.533 ms | 0.501 ms | 5.9% |
| 2^18 | (128,19,1) | 2.383 ms | 2.220 ms | 6.8% |
| 2^20 | (128,19,1) | 10.369 ms | 9.902 ms | 4.5% |

The generator lifts the existing XOR circuit to four independent rows per
AVX-512 vector. Routing consumes the corresponding interleaved tile layout.
The option retains the required BCH offset table after setup compaction.
The AVX2 fallback keeps its original two-row circuit.

## Wide elements

The second sweep compares four-row 128-bit encoding against the existing
256-bit kernel. Both encode the same 16-plane workload. Native wide inputs
already contain pairs of 128-bit planes; planar inputs require timed packing.

| K | 128-bit encode / 16 | Native 256-bit encode / 16 | 128-bit total | Native 256-bit total | Planar-to-256-bit total |
|---|---:|---:|---:|---:|---:|
| 2^16 | 0.503 ms | 0.377 ms | 9.555 ms | 7.511 ms | 8.333 ms |
| 2^18 | 2.246 ms | 1.747 ms | 43.435 ms | 34.886 ms | 39.955 ms |
| 2^20 | 9.777 ms | 8.215 ms | 186.049 ms | 160.051 ms | 181.841 ms |

Totals include all 16 planes and row-to-column assembly, plus packing where
needed. Assembly is an element transpose, not a bit transpose. Native 256-bit
inputs remain useful; converting planar inputs is only marginally beneficial
at K=2^20. The 128-bit control varies slightly between the two sweeps.

An explicit prefetch 64 elements ahead regressed the 128-bit encoding time by
about 13–25%, and native 256-bit time by about 2–9%. Do not enable
`SPIN_FORWARD_PREFETCH=64` for these configurations. The research switch defaults
to zero; no prefetch instructions are added to the normal build.

## Scope and reproduction

Host: Peach, Ryzen 9 7950X, GCC 15.2, CPU 15, `SPIN_TUNE=znver4`.
Each sweep uses three fresh processes per variant, seven measured trials per
process, and one warmup. The reported statistic is the median of process medians.
Variant order reverses on the second repeat. Benchmarks acquire both shared
locks and never run concurrently. Tile sizes are 256 rows at K=2^16 and 512
rows at K=2^18 and 2^20. The upstream snapshot is imported, never edited.

From the repository root, with the pinned upstream source at
`ROOT/hypercat/native/spin`:

```sh
bash workstreams/spin_optimized/forward_check.sh ROOT
bash workstreams/spin_optimized/forward_prefetch_check.sh ROOT
bash workstreams/spin_optimized/forward_validate.sh ROOT
python3 workstreams/spin_optimized/forward_summary.py
```

The prefetch script uses the control and four-row builds from `forward_check.sh`.
Raw samples and build logs remain ignored under `measurements/forward_*`.
No measurement data belongs in the source commit.

The optional build requires `SPIN_BIDIRECTIONAL_SOURCE`, `SPIN_BCH_AVX512=ON`,
and `SPIN_FORWARD_FOUR=ON`. Enable `SPIN_BUILD_WIDE` for wide workspaces.
The default supported target remains `spin_half_transpose`; forward still uses
the imported exponent-based interface, not the generalized transpose interface.

Release checks passed for dense forward references, transpose references,
adjoint identities, both index layouts, compaction, small tiles, and wide inputs.
The same three suites passed with AVX-512 detection disabled and under
AddressSanitizer/UndefinedBehaviorSanitizer. The object-code audit confirmed
that EVEX instructions remain confined to the designated AVX-512 objects:

```sh
python3 workstreams/spin_optimized/isa_check.py build-forward-four --optional-only
```

## Next step

Profile the remaining outer circuit and inner gather separately. Four-row
batching narrows the forward/transpose gap but does not close it, particularly
at small K. Prioritize reducing tile traffic or improving the forward BCH
schedule; the measured prefetch variant does not justify further tuning yet.
Keep the transposed API and kernels isolated from this work.
