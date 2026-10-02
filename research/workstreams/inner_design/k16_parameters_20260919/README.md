# Retuning the IMT inner at K = 2^16

The two-round `(64,12)` candidate now has a **full 49.328-bit bound** at
10% relative distance. Its measured 0.342–0.343 ms transposed encoding time
is essentially equal to the old `(128,19)` configuration. See
[FULL49_RESULT.md](FULL49_RESULT.md) for the proof components and performance
report, and [TARGET49.md](TARGET49.md) for the preceding search.

## Previous one-round result

The follow-up expansion-subspace search closed a new **full 45.382-bit bound**
for (t,s)=(64,12) at K=2^16 and 10% relative distance. It is measured about
2–4% faster than the current (128,19) implementation. This is a different
expansion map from the diagnostic (64,12) map tested below. See
[SUBSPACE_RESULT.md](SUBSPACE_RESULT.md) for the exact construction, complete
occupancy coverage, replay checks, and matched timings.

The default and submitted results are unchanged. The rest of this document
records the first screening round and its limitations.

## Initial screening result

The implemented (t,s)=(128,19) inner is sufficient, but not a demonstrated
runtime optimum at this length. The first measured alternatives give a small
speed improvement. None replaces its full 41.818-bit certificate yet.

All comparisons keep the rate-1/2 BCH [256,128] outer, cutoff 13107, independent
row/region permutations, one transvection per update, zero initial state,
output before update, persistent state, and no flush. The new maps are exactly
the deterministic adaptive-length study's map family. Changing t or s changes
the maps; these are not timings of the old map with a different loop size.

| Inner | Outward Q1 margin | Transpose time, ms | Change against matched control |
|---|---:|---:|---:|
| Current (128,19) | Existing full margin: 41.818 | 0.344–0.350 | control |
| (32,15) | 48.714 | 0.352–0.358 | slower |
| (64,19) | 44.833 | 0.358–0.363 | slower |
| (64,15) | 44.538 | 0.335–0.338 | 2.8–3.5% faster |
| (64,12) | 42.628 | 0.327–0.328 | 4.5–5.4% faster |

The Q1 columns for alternatives do not cover multiple active outer rows.
The (64,12) and (64,15) Q1 receipts passed 512-bit replay, using linear region
products instead of repeated squaring. Other new Q1 entries were produced
with outward arithmetic at 256 bits. These values use the general fixed-input
transfer and a six-tilt bank. They differ from the sharper, broader binary64
diagnostic study: 43.537 bits for (64,12) and 45.245 for (64,15).

Times above use unshared emission. Shared-XOR emission was also tested for
every alternative. It used fewer arithmetic operations but did not beat the
unshared implementation in any selected pair. All alternatives reuse the
new four-row BCH kernel and the K16 direct-scatter path.

## What is and is not proved

`model.py` parameterizes the existing outward fixed-input and independent-map
Fourier transfers. It reconstructs both map spectra, the feedback kernel,
fiber caps, and low-weight cancellation counts. It does not modify the frozen
proof-producing modules or reuse their numerical bounds for a new map.

Four tests passed: exact baseline transfer regression, independent binary64
formula comparisons at smaller steps, unchanged historical module globals,
and exact regression of the activation-density refinement. Exact baseline
Q1 comparisons check both ball midpoints and radii.

The multiple-row checks are incomplete:

- (32,15) clears the selected sparse checks q=2,4,8,16, but its current dense
  bound at q=218 and density coordinate 21/128 is larger than one.
- (64,15) clears the same sparse points. A coarse search initially failed at
  q=32; a finer tilt bank recovers a 97.518-bit bound there. At q=63, the
  fine-bank bound still exceeds one. This is a failure to certify, not a
  demonstrated low-distance word.
- (64,12) and (64,15) have strong bounds at the selected q=218 dense point.
  That point alone does not certify its neighborhood or all compositions.
- Reusing the old sparse witnesses without retuning fails. A 180-second
  evaluation of the old dense partition also stopped before evaluating every
  leaf and left weak bounds. No union from that partial evaluation is a
  full certificate. The historical density-subset refinements have not all
  been ported into this initial adapter.

The default therefore remains (128,19). No paper claim or selected certificate
was changed.

## A useful next map question

More updates per region address the short-length retention problem: h rises
from 4 at t128 to 8 at t64. The expansion spectra introduce another tradeoff.
The selected t128 expansion has minimum relative weight 48/128=3/8. The
diagnostic t64 maps contain expansion words of relative weight 16/64=1/4.
For s12 there are just three weight-16 words and one weight-48 word among the
4095 nonzero states. For s15 the corresponding counts are 22 and 18.

These counts suggest a bounded search for a better t64 expansion subspace,
not just more state. Can we suppress these extreme shells while keeping
s around 12–15 and retaining a cheap degree-two circuit? This may help the
retained-state bounds. The present results do not establish that those shells
cause the whole remaining gap. Test both Q1 and intermediate occupancies for
each new map before spending time on a complete proof or performance tuning.

## Performance evidence

Peach Ryzen 7950X, CPU 15, GCC 15.2, Release-style -O3, x86-64-v3 baseline
with znver4 scheduling and separately compiled AVX-512 kernels. Each cell has
three processes, each with 101 in-place calls after three warmups. Setup,
workspace allocation, and input initialization are outside timing. No input
reset/copy is included. Route seeds 1 and 17 were tested; the second repetition
reversed variant order. Benchmarks acquire the two shared locks and reject
other active benchmark executables. All runs were serial.

The first batch measured (32,15), (64,15), and (64,19); the second added (64,12)
with fresh baseline measurements. The local `measurements/remote` and
`measurements/round2` directories preserve these separate batches. Do not use
the second batch's control to calculate the first batch's improvements.

All eight variant builds passed the K16 dense-oracle test, including two
routing layouts, in-place suffix preservation, compaction, a second setup,
and impulses at epoch/region boundaries. Unshared/shared variants have equal
output hashes. The generated map-header hashes match the measured remote
headers. Circuit synthesis also checks every output as an exact binary linear
form. These experimental builds have not received a new sanitizer run.

Retained setup is 532480 bytes for the baseline, 540672 for t64, and 557056
for t32. Workspace remains 3 MiB. The extra setup stores more transvections.

## Files and reproduction

- `model.py`, `test_model.py`: explicit-step outward transfers and regression tests.
- `screen.py`: Q1, sparse points, dense points; `--verify` replays saved checks.
- `coverage.py`: bounded recomputation of historical search witnesses, not a
  claim that they work for the replacement maps.
- `refine.py`, `retune_sparse.py`: activation refinement and finer sparse tilt search.
- `generate_map.py`: degree-two expansion circuit and optional shared emission,
  both symbolically verified against the map columns.
- `overlay.py`: modifies fresh generated experimental build files only.
- `measure.sh`, `download.ps1`, `summary.py`: serial measurement and evidence checks.

Run Python commands from the repository root. Numerical dependencies are the
existing NumPy, SciPy, and python-flint environment. For example:

```text
python -B -m unittest discover -s workstreams/inner_design/k16_parameters_20260919 -v
python -B workstreams/inner_design/k16_parameters_20260919/screen.py --t 64 --s 12 --mode q1 --output <fresh.json>
python -B workstreams/inner_design/k16_parameters_20260919/screen.py --t 64 --s 12 --mode q1 --verify --output <fresh.json>
python -B workstreams/inner_design/k16_parameters_20260919/generate_map.py --t 64 --s 12 --output <map-directory>
```

The Linux timing driver expects the optimized `build-final-default` baseline
and generated map headers under this study's `measurements` directory. Invoke
`bash .../measure.sh <repository-root> "64_12" <fresh-batch-name>`.
It rejects an existing batch directory. Numerical receipts, generated headers,
build products, and raw measurements are ignored; commit the source and this
explanation, not the data. Nothing has been committed or pushed by this study.
