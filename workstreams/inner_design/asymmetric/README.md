# Independent expansion and feedback maps

The selected independent-map inner now has an outward certificate for the
quarter-rate BCH [128,32,32] instance at K=2^20, t=128, s=19. The 256-bit
calculation and 512-bit replay certify both operating points: 40 margin bits
at 16.5% relative distance and 30 bits at 19%. See `CERTIFICATE_RESULT.md` for
the precise claim, verification commands, and remaining migration gates.
The supported default remains unchanged.

The collision-aware weight-three feedback map `greedy3_2`, with its sparse
update kernel, reduced measured in-place transpose time from 17.037 ms to
16.100 ms (5.50%) in the original comparison. These historical timings exclude
the later routing optimization; see `../routing_opt/deployment/README.md`.

See `RESULT.json` for the bound and performance summary, `measurements/` for
raw timings and correctness logs, and `TRANSFER_ARGUMENT.md` for the changed
mathematics. `RESULT.json` retains the historical binary64 search results;
the outward result and its verification receipt are separate artifacts.

## Goal and initial search budget

Keep the quarter-rate [128,32,32] outer and initially fix the certified balanced
expansion A at t=128, s=19. Search independent feedback maps B with cheap sparse
columns. Retain the one-round transvection mixer as the first control.

The first stage compares the symmetric baseline, weight-two columns with a
rank-restoring odd column, random weight-three columns, collision-aware
weight-three columns, and mixed weight-two/three columns. Use two random seeds
for the weight-three families. This bounded screen is not an optimality claim.

For each candidate, audit rank, zero/duplicate columns, low-weight kernel words,
the full spectrum of B^T, and exact cancellation for one/two input positions.
Screen the single-active-row bound before more expensive occupancy calculations.
Advance a small shortlist to isolated C++ kernels and serial measurements only
after the mathematical screen. All measured claims require complete-encoder
timings; arithmetic counts alone are discovery metrics.

## Exact interface

Over F_2, A maps s state coordinates to t output coordinates and B maps t
input coordinates to s state coordinates. At epoch i,

    Y_i = X_i + A q_i,
    q_(i+1) = M_i q_i + B X_i,       q_0 = 0.

Use the same independent per-epoch mixer setup and SPIN permutations as the
balanced baseline. Do not impose B=A^T or BA=0 in this exploration.
The input/output map is still causal with identity diagonal blocks, hence
invertible as a length-preserving inner map; this does not require BA=0.

For transpose input U_i and reverse state r_(i+1), initialized to zero at the
last epoch, the exact adjoint recurrence is

    V_i = U_i + B^T r_(i+1),
    r_i = A^T U_i + M_i^T r_(i+1).

Thus sparse B columns make transpose emission cheap. A^T acts on the raw U_i.
If an implementation instead applies A^T to the emitted V_i, it must subtract
(XOR) A^T B^T r_(i+1). The old kernel omitted this correction because BA=0.
Reusing that shortcut for independent maps would implement a different code.

## Proof requirements

The weighted-state envelope retains classes defined by wt(Aq). Its emitted
moment depends on the A spectrum. Its zero-syndrome probabilities and fiber
caps instead depend on B. Compute the kernel spectrum by MacWilliams inversion
of the separately enumerated B^T spectrum, not the A spectrum.

The generic fixed-input-weight transfer remains valid using these separate
inputs. Its exact low-input cancellation term is wt(X+A BX); it does not
require symmetry or BA=0. The existing all-one-aware dense Fourier shortcut
does require replacement: wt(B^T a) and intersections with Aq no longer follow
from the A spectrum alone. The binomial mixture of valid fixed-weight transfers
remains an available dense bound.

Failure screens include never-activated states, repeated syndromes, short
kernel words, lazy cancellation, and full-input/end-of-stream cases. A rank
deficiency is flagged as outside the initial full-rank search, not misreported
as a counterexample to the full SPIN code. A failed upper bound likewise does
not prove that the candidate lacks distance.

## Search results

Seven exact map audits are retained in `FEEDBACK_SCREEN.json`. All selected
maps have rank 19 and distinct nonzero columns. Some representative results:

| Feedback map | Direct transpose state-sum XORs | Kernel words of weight 3 / 4 | Minimum B^T image weight |
|---|---:|---:|---:|
| Symmetric control | 1016 | 0 / 0 | 48 |
| Weight two plus one odd column | 127 | 396 / 3578 | 1 |
| Half weight two, half weight three | 192 | 93 / 779 | 11 |
| Random weight three, seed 1 | 256 | 0 / 591 | 16 |
| Collision-aware weight three, seed 2 | 256 | 0 / 291 | 19 |

The XOR column counts a direct sparse implementation, excluding the common
input XOR. The symmetric baseline uses optimized lookup tables, so 1016 is
not its actual executed XOR count. These are not speedup measurements.

The weight-two-plus-one map illustrates why full rank is insufficient: one
dual direction observes only one input coordinate. The mixed and collision-aware
weight-three maps are the first shortlist; the cheapest map needs separate
scrutiny before treating its low operation count as useful.

The fixed-tilt Q=1 diagnostic is essentially unchanged for all seven maps.
`OCCUPANCIES_Q8.json` retains the small-occupancy control comparison.
`OCCUPANCIES_Q128.json` checks the complete sparse range for both shortlisted
maps: Q=1..64 at 16.5% distance and Q=1..128 at 19% distance.

The original dense witnesses did not pass unchanged. `DENSE_SCREEN.json`
records that diagnostic failure, not a distance counterexample. Retuning and
subdividing the cover closed every remaining occupancy, as recorded in the
four `DENSE_<candidate>_<distance>.json` files:

| Feedback | Distance | Dense cover leaves | Dense-only margin | Full-union margin |
|---|---:|---:|---:|---:|
| Mixed weight two/three | 16.5% | 141 | 115.70 bits | 41.048168 bits |
| Mixed weight two/three | 19% | 320 | 85.46 bits | 30.033491 bits |
| Collision-aware weight three, seed 2 | 16.5% | 139 | 117.78 bits | 41.048168 bits |
| Collision-aware weight three, seed 2 | 19% | 318 | 107.89 bits | 30.033491 bits |

The entries in this search table are binary64 evaluations of upper-bound formulas. Each
dense cover is an exactly checked partition, with zero unresolved leaves.
Combined with the sparse calculation, it covers Q=1..32768. This is materially
stronger evidence than sampling occupancies. The subsequent outward replay
certifies the selected `greedy3_2` map, not the mixed map.

Tests check the exact transpose identity on every bilinear basis pair of a
small three-epoch example with BA nonzero. They also verify that the old
shortcut fails on that example. Separate tests compare low-kernel counts to
enumeration and the asymmetric weighted transfer to exact rational state laws
for every small-example input weight. Bernoulli-transfer tests include small,
central, and large input densities and arbitrary source atoms/uniform shells.

## Implementation and performance

Four isolated kernels combine the two shortlisted maps with sparse or masked
transvection updates. All passed the independent dense-oracle correctness suite
at K=2^16,2^18,2^20, including layouts, tiling, in-place output, suffix
preservation, compaction, and linearity. They remain outside the supported
encoder. Generated source and map hashes are in `IMPLEMENTATION.json`.

All benchmarks ran serially on Peach, Ryzen 9 7950X CPU 15, GCC 15.2.0,
Release -O3 with znver4/AVX2 flags. Each row below is the median of three
101-trial run medians, with order varied between repeats. The workload uses
K=2^20, the quarter-rate outer, 128-way bitslicing, packed24 routing, and
4096-row tiles. Setup is excluded; encoding is in place without a timed copy.

| Inner map/update | Complete encoder | Time reduction vs balanced |
|---|---:|---:|
| Certified balanced / masked | 17.037 ms | control |
| Mixed weight two/three / sparse | 16.334 ms | 4.13% |
| Mixed weight two/three / masked | 16.775 ms | 1.53% |
| Collision-aware weight three / sparse | 16.100 ms | 5.50% |
| Collision-aware weight three / masked | 16.525 ms | 3.00% |

The selected kernel's three run medians were 16.098, 16.177, and 16.100 ms.
All candidates retained the same 25,427,976 bytes of setup and 75,497,472
bytes of workspace as the balanced control. These are one-host measurements,
not a claim about every processor or parameter choice.

Standalone inner diagnostics exclude routing and the outer but retain real
contiguous output stores. For the selected map, the smaller 32,768-block test
drops from 0.035706 to 0.027812 ms, about 22%. At 4,194,304 blocks (64 MiB
input), it drops from 5.210277 to 5.049708 ms, about 3%. The streaming runs
show more run-to-run variability. These standalone times are not additive
components of the fused encoder, so they must not be subtracted from total
time to infer exact routing or outer costs.

The redesign therefore offers a useful but modest end-to-end improvement.
Cheaper emission alone does not remove feedback computation or memory traffic.
The mixed map uses fewer direct XORs but loses to the triple map in the complete
encoder; operation counts are not enough to rank candidates.

## Completion and next step

The bounded search, exact map audits, independent-map transfer derivation,
all-occupancy numerical screen, isolated implementation, correctness tests,
serial full-encoder timings, and kernel diagnostics are complete. The CWC
writing workflow keeps the exact identities, numerical bounds, and remaining
certification obligations separate in `TRANSFER_ARGUMENT.md`.

The selected `greedy3_2_sparse` candidate has now passed outward production,
higher-precision replay, exact union and coverage checks, and binding to its
tested map and outer implementation. This closes the first finite-instance
gate only. Next extend certificates to the paper's BCH-256 cells, adapt the
asymptotic argument, and measure a complete forward encoder, as specified in
`MIGRATION_PLAN.md`.

## Reproduce

For the new outward certificate, follow the commands in
`CERTIFICATE_RESULT.md`. The commands below reproduce the earlier search and
performance exploration; rerunning witness discovery is not required for
certificate replay.

From the repository root, using the dependencies in `../requirements.txt`:

```text
python workstreams/inner_design/asymmetric/search_feedback.py
python -m unittest discover -s workstreams/inner_design/asymmetric -p "test_*.py" -v
python workstreams/inner_design/asymmetric/screen_occupancies.py
python workstreams/inner_design/asymmetric/screen_occupancies.py --maximum 128 --names mixed2_3 greedy3_2
python workstreams/inner_design/asymmetric/screen_dense.py
python workstreams/inner_design/asymmetric/refine_dense.py --name greedy3_2 --distance 33/200 --nodes 48
python workstreams/inner_design/asymmetric/refine_dense.py --name greedy3_2 --distance 19/100 --nodes 48
python workstreams/inner_design/asymmetric/refine_dense.py --name mixed2_3 --distance 33/200 --nodes 48
python workstreams/inner_design/asymmetric/refine_dense.py --name mixed2_3 --distance 19/100 --nodes 48
python workstreams/inner_design/asymmetric/generate_kernels.py
```

Witness search is floating-point discovery; exact witness bytes can vary
across environments. The retained covers are the evidence for this run.
Generated kernels go in `../generated/asymmetric_*`; other outputs stay here.

On a Linux host with the repository sources present, build the isolated
targets using `../CMakeLists.txt`, not the supported encoder target:

```text
cmake -S workstreams/inner_design -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --target basis_balanced_masked_bench basis_asymmetric_greedy3_2_sparse_bench basis_asymmetric_greedy3_2_sparse_test basis_asymmetric_greedy3_2_masked_bench basis_asymmetric_greedy3_2_masked_test basis_asymmetric_mixed2_3_sparse_bench basis_asymmetric_mixed2_3_sparse_test basis_asymmetric_mixed2_3_masked_bench basis_asymmetric_mixed2_3_masked_test -j3
bash workstreams/inner_design/asymmetric/run_bench.sh "$PWD"
cmake -S workstreams/inner_design/asymmetric -B inner-build -DCMAKE_BUILD_TYPE=Release
cmake --build inner-build -j3
bash workstreams/inner_design/asymmetric/run_kernel_bench.sh "$PWD"
```

Run the two scripts sequentially, never concurrently. Every benchmark acquires
the shared lock and checks for other benchmark processes. CPU 15 and the
znver4 target are explicit measurement settings; adapt both for another host.
Save the scripts' `asymmetric-*.jsonl`, correctness log, and binary-hash receipts
under `measurements/` when collecting a new run. Do not mix receipts from
different builds. Then run the fast integrity/summary check:

```text
python workstreams/inner_design/asymmetric/summarize.py
python workstreams/inner_design/verify_balanced_artifact.py
```

The first checks source bindings, reconstructs integer map audits, validates
complete cover geometry, and summarizes measurements. It is not an outward
proof verifier. The second confirms the earlier balanced certificate remains
bound to its unchanged implementation and higher-precision replay receipt.
