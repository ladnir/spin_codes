# A cheaper inner with full finite-length coverage

The balanced-image, one-round mixer retains both quarter-rate operating points
at K=2^20. Its full in-place transposed encoder takes about 17.14 ms on Peach,
versus 17.85 ms for the unchanged baseline in the same serial comparison.
That is about 4% less online time, with 2.125 MiB less retained setup.

The 256-bit outward calculation passes all occupancies. A 512-bit replay
encloses every recomputed term within the saved upper bounds. The isolated
implementation is not yet the production default.

## Construction and certified statement

Use the fixed [128,32,32] outer from
[`smaller_outer.py`](../rate_quarter_bch/smaller_outer.py). There are L=32768
outer rows, K=2^20 message bits, and N=2^22 output bits. Independently permute
each outer row's 128 coordinates into the regions, then independently permute
the L positions within each region. These are the retained SPIN permutations.

The inner has t=128 and s=19. Its B columns are exactly those in
[`NO_CONSTANT_MAP.json`](NO_CONSTANT_MAP.json), and A=B^T. Each epoch takes
input X and state q, emits Y=X+Aq, and updates q'=Mq+BX. The state starts at
zero and continues across regions. The output excludes the final state.

For each epoch, sample u uniformly from the nonzero 19-bit vectors, then sample
v uniformly from u's orthogonal complement, including zero. Set M=I+uv^T.
Samples are independent across epochs and independent of the permutations.
One sampled setup defines the linear encoder for all messages.

For this ensemble, let E_d be the event that some nonzero message has output
weight at most floor(dN). The all-occupancy first-moment calculation gives:

| Relative-distance target d | Bad-weight cutoff | Full margin, approximately | Certified target |
|---|---:|---:|---:|
| 16.5% | 692,060 | 41.048169 bits | Pr[E_d] < 2^-40 |
| 19% | 796,917 | 30.033492 bits | Pr[E_d] < 2^-30 |

Both rows describe the same encoder, not different parameter choices. Exact
upward dyadics, rather than the displayed decimal margins, establish the
inequalities. The statement is finite-length at K=2^20; the runtime correctness
tests at smaller K do not extend this certificate automatically.

The [full certificate](NO_CONSTANT_MARGIN_CERTIFICATE.json) and
[higher-precision replay receipt](NO_CONSTANT_MARGIN_CERTIFICATE_REPLAY.json)
retain the numerical evidence.

## Why the proof closes

One mixer round sends any fixed nonzero state to half its original point mass
and half the uniform nonzero distribution. This exact marginal law is proved
in [PROOF_REQUIREMENTS.md](PROOF_REQUIREMENTS.md). It is not a claim that
different messages acquire independent state trajectories.

The unchanged map contains an all-one image. Its interaction with the lazy
half of the mixer made the dense bounds difficult. The new map has all nonzero
image weights in [48,80], including every nonzero difference of two image words.
[SPECTRUM_REDESIGN.md](SPECTRUM_REDESIGN.md) gives the exact construction.
The kernel distance decreases from six to five; the proof recomputes the kernel
spectrum and cancellation bounds rather than ignoring that cost.

Fix a nonzero message and weight each trajectory by z to its emitted weight,
where 0<z<1. The transfer bounds this weighted state measure, including state
zero, arbitrary nonzero mass, and pointwise caps within image-weight classes.
It retains weight-class information on empty epochs. On occupied epochs it
explicitly bounds lazy cancellation using syndrome fibers. Thus the proof does
not apply a uniform-state cancellation estimate after the emitted weights have
biased the state law. [FIBER_BOUNDS.md](FIBER_BOUNDS.md) derives the fiber caps
and the alternative transfer for the dense Bernoulli surrogate.

The union over nonzero messages is partitioned by the number Q of active outer
rows. The certificate covers Q=1 separately, Q=2 by band compositions, and
Q=3..64 or Q=3..128 by adaptive sparse bounds. The dense covers then include
every remaining Q through L, using 135 and 308 integer count boxes respectively.
Exact geometry checks exclude gaps and overlaps. The proof never resets the
inner state at a region boundary or charges un-emitted terminal weight.

Q=1 dominates both final unions. Its separate margins are 41.048245 and
30.034069 bits, while the dense contributions have margins above 107 and
99 bits. Further dense-bound tightening would barely change the final margin.
The 19% target still has only about 0.033 bits beyond its requested 30-bit
margin; this result does not claim large slack at that operating point.

`dense_cover.py` searches for witnesses in binary64. `certify_no_constant.py`
imports only fixed witnesses and recomputes transfers, counting costs, region
powers, box bounds, and the complete union in Arb. The dense bound may use
either of two complete valid moment bounds; it does not mix incompatible
state representations. The exact partition and all source hashes are checked.

The main failure modes are therefore addressed explicitly: zero-syndrome
inputs, lazy cancellation, weight-biased states, the all-one image, dense
occupancies, and unflushed termination. The old-map one-round candidate still
has no full certificate. That failed bound is not a counterexample to its code.

## Performance of the exact certified map

Peach, Ryzen 9 7950X, CPU 15, GCC 15.2, `-O3 -DNDEBUG -march=znver4` with
AVX2 and carry-less-multiply flags. The workload is the complete 128-way
bitsliced transposed encoder, in place, packed24 routing, and 4096-row tiles.
Setup, allocation, input copying, and initialization are outside the timer.

| Implementation | Three 101-trial run medians (ms) | Median of medians |
|---|---|---:|
| Unchanged baseline | 17.847, 17.633, 17.928 | 17.847 ms |
| Balanced map, sparse update | 17.140, 17.135, 17.085 | 17.135 ms |
| Balanced map, fixed masked update | 17.145, 17.181, 17.041 | 17.145 ms |

Every candidate run improves on its paired baseline. Run order changes across
repetitions, and each executable holds the shared benchmark lock. Sparse and
masked updates differ by only 0.010 ms in the aggregate; these measurements
do not establish a meaningful winner. Fixed masked updates offer predictable,
unrolled control flow at effectively the same measured cost.

Retained setup falls from 27,656,200 to 25,427,976 bytes. Workspace stays at
72 MiB. [BALANCED_PERFORMANCE.json](BALANCED_PERFORMANCE.json) and the
[raw receipts](measurements/balanced/) retain the timings and source bindings.
The earlier original-map mixer timings are a separate experiment, not evidence
substituted for this final map.

`generate_balanced.py` symbolically checks the complete feedback circuit against
the map columns and verifies the grouped lookup coordinates. `MixerInner.h`
reuses the emission table for the rank-one dot product. The implementation has
no per-epoch allocation or dynamic dispatch. Two isolated C++ variants pass the
dense oracle at K=2^16, 2^18, and 2^20, with independent setups, boundaries,
linearity, layouts, tile sizes, compaction, and in-place suffix preservation.
The oracle uses raw columns and raw u/v samples, not the optimized circuit.

## Reproduction and disposition

From the repository root, with Python, NumPy, SciPy, and python-flint installed.
The tested versions are recorded in [requirements.txt](requirements.txt)
(CPython 3.14.4, NumPy 2.4.4, SciPy 1.18.0, python-flint 0.9.0).

```text
python -m unittest discover -s workstreams/inner_design -p "test_*.py" -v
python workstreams/inner_design/certify_no_constant.py --verify
python workstreams/inner_design/generate_balanced.py
python workstreams/inner_design/summarize_balanced.py
python workstreams/inner_design/verify_balanced_artifact.py
```

The replay uses 512-bit arithmetic against the saved 256-bit receipt. Omit
`--verify` only when deliberately producing a replacement certificate. Witness
discovery is unnecessary for replay; the saved covers are part of the artifact.
The final command binds the certificate and replay to the exact generated map,
outer code, tested sources, and timing receipts. It also checks the dyadic union
with exact integer arithmetic on a common dyadic scale and verifies every
weight-band partition.

On the benchmark host:

```sh
cmake -S workstreams/inner_design -B out/inner-design -DCMAKE_BUILD_TYPE=Release
cmake --build out/inner-design --target basis_balanced_sparse_test basis_balanced_masked_test basis_balanced_sparse_bench basis_balanced_masked_bench -j2
ctest --test-dir out/inner-design -R balanced --output-on-failure -j1
```

`run_balanced_bench.sh` accepts an isolated root containing `baseline/spin_benchmark`
and the candidate `build/` directory. It records all nine runs serially. Do not
hold an external copy of the lock when starting these executables.

The exact-equivalence search is complete at its practical stopping criterion;
it found no repeatable improvement. This second track has a full finite-length
certificate, a successful higher-precision replay, and a repeatable gain for a
changed inner. The two-part investigation is complete. Next, integrate the
selected variant and proof statement deliberately. Larger t, different state sizes, and wider K ranges are
future investigations, not prerequisites for this operating point.
