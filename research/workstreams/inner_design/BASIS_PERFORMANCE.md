# Exact-encoder basis search

No repeatable material gain was found. The baseline remains the supported
implementation. This closes a practical optimization tier, not all possible
equivalent implementations.

## Equivalence

Let `r` be the original reverse state and put `z=V r`, where `V` is invertible.
The original emission, feedback and transition maps are `A`, `A^T` and `M^T`.
The replacement maps are respectively

```text
A V^-1,       V A^T,       V M^T V^-1.
```

Feedback is **not** reset to the transpose of the new emission map. That would
generally change the encoder. The setup conjugates every transition; the hot
kernel uses generated emission/feedback circuits. The original dense oracle
still uses the original maps and field coefficients.

## Search and verification

Two searches used 30 restarts each, with 30,000 and 100,000 moves per restart:
3.9 million moves altogether. Starting bases included grouped, feedback-RREF
and emission-RREF. The objective traded emission-table lookups against feedback
circuit cost; only full-encoder timings determined adoption.

Nine representative bases were compiled. All nine passed the full original
dense-oracle suite at message sizes 2^16, 2^18 and 2^20. Python tests independently
check emission/feedback identities and transition conjugacy on a bilinear basis.
The search receipts retain their historical producer hashes; tests validate
their matrices rather than claiming a replay of an earlier script version.

## Full-encoder timings

Peach Ryzen 9 7950X, CPU 15, GCC 15.2, release `-O3 -DNDEBUG -march=znver4`.
Quarter-rate [128,32,32], t=128, s=19, K=2^20, 128-way bitsliced transposed
encoding, in-place, packed24 layout, tile size 4096. Setup is outside the timer.
Each row below contains three serial runs of 101 trials. Each run block measured
the baseline followed by the four finalists. No benchmarks ran concurrently.

| Basis | Three run medians (ms) | Median of medians (ms) |
|---|---|---:|
| Unchanged baseline | 17.596, 17.475, 17.766 | 17.596 |
| Feedback-RREF | 17.734, 17.693, 17.719 | 17.719 |
| Joint, weight 0.5 / restart 11 | 18.097, 18.006, 18.151 | 18.097 |
| Joint, weight 1 / restart 13 | 17.839, 17.909, 17.921 | 17.909 |
| Joint, weight 2 / restart 10 | 17.733, 17.587, 17.550 | 17.587 |

The apparent 0.05% gain in the last row is far smaller than run-to-run variation.
Every repeated run has the same output hash. Setup conjugation also increased
setup time; no online improvement justifies that extra work here.

The historical reported baseline is 17.996 ms. These fresh control runs do not
represent a code speedup. Prior failed low-level experiments are recorded in
[OPTIMIZATION_NOTES.md](../rate_quarter_bch/implementation/OPTIMIZATION_NOTES.md).
Their stage breakdown includes routing with the inner; it is not a pure
inner-arithmetic profile.

## Reproduction

```text
python workstreams/inner_design/search_basis.py
python workstreams/inner_design/search_basis.py --steps 100000 --restarts 6 --seed-offset 10 --mixed-starts --output BASIS_SEARCH_SECOND.json
python workstreams/inner_design/generate_candidates.py
python workstreams/inner_design/generate_candidates.py --search BASIS_SEARCH_SECOND.json joint_w0.5_r11 joint_w1_r13 joint_w2_r10
cmake -S workstreams/inner_design -B /tmp/inner-basis-build -DCMAKE_BUILD_TYPE=Release
cmake --build /tmp/inner-basis-build -j 4
ctest --test-dir /tmp/inner-basis-build --output-on-failure
```

Example **single** timing, only after confirming no other benchmark is running:

```text
taskset -c 15 /tmp/inner-basis-build/basis_joint_w2_r10_bench 101 0 0 20 1 2 1
```

The executable acquires `/tmp/bare-spin-benchmark.lock` itself and checks for
other benchmark processes. Do not wrap it in a second lock on the same file.
Do not launch that command concurrently with any benchmark. Baseline and
candidate must use the same arguments and machine. The compact retained
[measurements](measurements/) can be summarized without running a benchmark:

```text
python workstreams/inner_design/summarize_measurements.py
```
