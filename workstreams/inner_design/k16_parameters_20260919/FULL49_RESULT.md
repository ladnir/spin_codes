# K16: 49.33 bits at essentially the original runtime

Integration follow-up: the reusable libraries now expose `T64S12R2`, including
block forward and bit-packed forward encoding. See
[the integration documentation](../../spin_optimized/README.md#two-round-k16-option-4933-bits).
The prototype build and timings below record the earlier isolated experiment;
use `workstreams/spin_optimized` for current builds.

The selected `(t,s,r)=(64,12,2)` IMT inner closes a full **49.3275773868-bit**
setup-failure bound at `K=65536`, rate 1/2, and relative distance 10%.
Its measured transposed encoding time is **0.342–0.343 ms** on Peach.
The old `(128,19,1)` configuration takes 0.344–0.345 ms in the same batch.
Treat these runtimes as essentially equal, not a demonstrated speedup.

## Construction and bound

The outer is the existing deterministic BCH [256,128] code, repeated over
512 rows. The independent row and region permutations are unchanged.
The fixed expansion and feedback maps are exactly those in
`../../spin_optimized/k16_map.py` and the one-round subspace result.

For input group x and state a, the inner emits x+A a and then updates
a to T2 T1 a+B x. T1 and T2 are independently sampled transvections for
that update. Setup is sampled once and reused for all messages.
The state starts at zero, persists across regions, and has no final flush.
There is only one feedback addition and no output between the transvections.

The probability is over this setup distribution. The bound covers the event
that some nonzero message produces an output of weight at most 13107,
among 131072 output positions. It uses deterministic BCH shell bounds,
not a guessed BCH spectrum.

| Component | Active outer rows | Margin, bits |
|---|---:|---:|
| Sharp one-row transfer | 1 | 49.3295959765 |
| Sparse transfer | 2–63 | 72.9999999997 |
| Dense cover | 64–512 | 58.8098666781 |
| Exact sum | 1–512 | **49.3275773868** |

The sparse calculation reevaluates all 62 occupancies with 11 retained
witnesses. The dense calculation reevaluates the old 1112 rectangles and
adds one split, giving 1113 leaves. Retuning and splitting the largest
contribution explain why the final dense margin exceeds the earlier
ten-leaf spot-check minimum.

All components passed 512-bit replay. The one-row replay uses linear
region products. The final assembler authenticates source and receipt
hashes, checks the dense partition, reconstructs each contribution, and
sums exact rational upper bounds. It checks that the sum is below 2^-49.
The assembled receipt is
`measurements/target49_sub12_r2_full.json`.

The transfer justification is in [TARGET49.md](TARGET49.md). On nonzero
states, two rounds replace P by P²=P/2+J/2. The new fixed-weight transfer
combines positive majorants using the same source-measure representation.
The Bernoulli alternative uses a conservative whole-row comparison.
No numerical one-round failure bound is assumed to remain valid.

## Matched performance

These are medians of five process medians, each with 101 timed encodings
and three warmups. Runs were serial, pinned to CPU 15 of Peach's Ryzen
7950X, using GCC 15.2, AVX-512 BCH, and the direct K16 routing path.
Each encoding operates in place on 128-bit blocks. Setup, allocation,
and input copying/resetting are excluded.

| Configuration | Full margin, bits | Seed 1, ms | Seed 17, ms | Retained setup |
|---|---:|---:|---:|---:|
| Old `(128,19,1)` | 41.818 | 0.344764 | 0.344523 | 520 KiB |
| Selected `(64,12,1)` | 45.382 | 0.332751 | 0.334294 | 528 KiB |
| Selected `(64,12,2)` | **49.328** | **0.343311** | **0.342019** | 544 KiB |

All three use 3 MiB of workspace. The second round adds 16 KiB of
retained setup to the one-round selected map. Its time penalty is 2.3–3.2%
against that faster option, while recovering about 3.95 bits of margin.
Against the old configuration, it gains about 7.51 bits at similar cost.

This is an isolated standalone-transpose prototype under `r2bench/`.
The second transvection is explicitly expanded in the hot loop; there
are no new allocations or runtime abstractions inside the inner.
The transpose applies T2^T and then T1^T before adding A^T times its input.

Five release tests and the same five ASan/UBSan tests passed. They include
full dense-oracle comparisons, AVX2/AVX-512 paths, layouts, tile sizes,
compaction, in-place behavior, and a separate forward/transpose adjoint
test for 20 sampled setups. Exact A/B constants match the assembled
certificate. The prototype does not yet expose a complete bidirectional
public API; the separate forward recurrence is a correctness oracle.

Raw timings, source and binary hashes, and test logs are ignored under
`../../spin_optimized/measurements/k16-r2/`.

## Reproduction

Run from this directory with Python, NumPy, SciPy, and python-flint.
Existing seed receipts are local research inputs, not part of the
code-only submission artifact. Outputs must not already exist.

```powershell
python target49_rounds_cover.py --mode sparse --maps measurements/subspaces_v1.json --seed measurements/sub12_sparse_final.json --output measurements/target49_sub12_r2_sparse.json
python target49_rounds_cover.py --mode sparse --maps measurements/subspaces_v1.json --seed measurements/sub12_sparse_final.json --output measurements/target49_sub12_r2_sparse.json --verify
python target49_rounds_cover.py --mode dense --maps measurements/subspaces_v1.json --seed measurements/sub12_dense_cover.json --output measurements/target49_sub12_r2_dense.json
python target49_rounds_cover.py --mode dense --maps measurements/subspaces_v1.json --seed measurements/sub12_dense_cover.json --output measurements/target49_sub12_r2_dense.json --verify
python target49_assemble.py --q1 measurements/target49_sub12_r2_q1.json --sparse measurements/target49_sub12_r2_sparse.json --dense measurements/target49_sub12_r2_dense.json --output measurements/target49_sub12_r2_full.json
```

The earlier Q1 producer/replay is `target49_rounds_q1.py`. Build the
isolated encoder with CMake using `r2bench/` as the source directory,
`SPIN_BCH_AVX512=ON`, `SPIN_TUNE=znver4`, and `SPIN_BUILD_BENCHMARK=ON`.
Run CTest before invoking `r2bench/measure.py`. The measurement script
expects the local experiment's named control and candidate build directories.
It alternates case order and never runs two benchmarks concurrently.

## Recommendation

Prefer this candidate over another small one-round speed improvement when
the application benefits from a 49-bit margin. The proposed reusable
bidirectional integration has now been implemented. Next, select this
explicit option in the downstream application and measure its effect.
Do not extend the bound to other message lengths.
Defaults, submitted paper/artifact, and K20 configurations remain unchanged.
