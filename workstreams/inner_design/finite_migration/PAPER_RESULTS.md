# Verified finite IMT results for paper integration

All selected instances use t=128, s=19, and one independently sampled
transvection per epoch. The fixed expansion excludes the all-one word.
Half-rate instances use weight-five feedback; quarter-rate instances use
the separately specified weight-three feedback map `greedy3_2`.
The inner outputs before updating state, starts at zero, persists across
regions, and discards its final state without a flush. No parity fanout is used.

## Complete certificates

The margin is minus the base-two logarithm of a rational upper bound on
setup failure. Failure means that some nonzero message has output weight
at most floor(delta N). These are not estimates of the true failure probability.

| Outer | Rate | log2(K) | delta | Full margin (bits) | Verification receipt |
|---|---|---:|---:|---:|---|
| BCH [256,128,d>=38] | 1/2 | 16 | .10 | 41.8183267960 | `FULL_M16_VERIFIED.json` |
| Same | 1/2 | 18 | .10 | 50.1890763816 | `FULL_M18_VERIFIED.json` |
| Same | 1/2 | 20 | .10 | 50.0620882264 | `../asymmetric/bch256/weight5/FULL_M20_VERIFIED.json` |
| Same | 1/2 | 22 | .10 | 48.3904916829 | `FULL_M22_VERIFIED.json` |
| Same | 1/2 | 24 | .10 | 46.4572549297 | `FULL_M24_VERIFIED.json` |
| BCH [128,32,32] | 1/4 | 20 | .165 | 41.0481676058 | `QUARTER_DEPLOYMENT_VERIFIED.json` |
| Same | 1/4 | 20 | .19 | 30.0334910634 | Same |

The quarter-rate receipt authenticates the certificate, replay, and optimized
implementation. The half-rate proof receipts are complemented by the timing
binding below. No proof is inferred at unlisted lengths.

K16 uses Q1, Q2..63, and Q64..512 contributions. Their margins are
43.5928335021, 60.9887506349, and 42.3171086801 bits. The dense cover has
1,096 leaves and passed 512-bit replay. Its coupled normalizer calculation
is derived in `SHORT_LENGTH_BOUNDS.md`; the actual union is checked by
`verify_budget.py`.

The K16 margin is not on a straight extrapolation of the larger-length
certificates. The paper should plot the certified points without implying
a universal one-bit-per-doubling law or a guarantee between those points.

## Matched transpose performance

Each timing evaluates the complete transposed encoder on 128 parallel binary
instances. Setup and workspace allocation are excluded. Buffers are initialized
once per process, with no copy or reset between calls. Measurements are serial
on Peach CPU 15, a Ryzen 9 7950X, with GCC 15.2 and znver4 compilation.
Each cell is the median of three process medians, each from 101 calls after
three warmups. Quarter-rate timings retain their separately authenticated series.

| Rate | log2(K) | Selected IMT latency (ms) | Retained setup (bytes) | Workspace (bytes) |
|---|---:|---:|---:|---:|
| 1/2 | 16 | 0.523367 | 794632 | 3145728 |
| 1/2 | 18 | 2.153740 | 3178504 | 9437184 |
| 1/2 | 20 | 10.109862 | 12713992 | 41943040 |
| 1/4 | 20 | 15.254337 | 25427976 | 75497472 |

`HALF_LADDER_PERFORMANCE_VERIFIED.json` binds all three half-rate cells to
their exact maps and full distance certificates. `QUARTER_DEPLOYMENT_VERIFIED.json`
binds the quarter-rate cell and both distance thresholds. No K22/K24 timing is
claimed. The half-rate K20 process medians range from 10.095175 to 10.172329 ms.

The fresh matched RM2Sub baseline is 0.562380, 2.305955, and 11.138372 ms
at the three half-rate sizes. These may be presented as an explicitly
historical-inner comparison. They are not IMT measurements. The selected IMT
reductions are 6.94%, 6.60%, and 9.23%, respectively. The older 11.259 ms
headline and its setup-time measurement must not be carried forward as IMT.

## Manuscript integration and reproduction

The manuscript now uses these selected IMT results: the finite construction,
five half-rate margins and engineering curve, quarter-rate theorem and map,
implementation table, external-comparison SPIN row, abstract, and introduction.
The shared transfer argument is reused from the asymptotic IMT appendix;
the finite appendix adds the short-length and quarter-rate reductions.

The selected integration check authenticates 746 files, seven certificate
targets, four timing cells, and 57 map words. It verifies exact half-rate
sums and outward-rounded quarter-rate union ceilings, not new interval
evaluations. Six fast rejection tests cover the new evidence adapter.

```text
python -B paper/check_finite_integration.py
python -B paper/check_imt_integration.py
python -B paper/build_imt_comparison.py --check
python -B -m unittest discover -s paper -p test_imt_results.py
```

The combined comparison generator preserves the external campaign. It
authenticates two subsequently changed shared SPIN files against their
measured bytes in commit `aafb3f59e7c3541e19b8623520c8a407ee2219e8`.
The old generator remains unchanged and is not the current table entry point.

The integrated draft compiles to 60 pages. The updated front matter,
certificate tables, engineering curve, performance/comparison tables, and
finite appendix were rendered and visually inspected. The final log has
no overfull boxes, undefined references, or multiply-defined references.
The selected/asymptotic integration checks, 62 finite-migration tests,
6 evidence-adapter tests, and 8 artifact regression tests pass.
The artifact tests explicitly reject using the old compact RM2Sub receipts
as replacements for absent IMT evidence. No new benchmark or commit was
performed during this integration milestone.

The retained K16 proof and half-rate timing binding can be checked with:

```text
python -B workstreams/inner_design/finite_migration/verify_budget.py --m 16 --q1 workstreams/inner_design/finite_migration/Q1_LOWER.json --sparse workstreams/inner_design/finite_migration/SPARSE_M16_v2.json --dense workstreams/inner_design/finite_migration/DENSE_M16_coupled_v1.json
python -B workstreams/inner_design/finite_migration/half_ladder_binding.py --output <fresh-binding.json>
```

These commands authenticate retained replays and reconstruct exact unions;
they do not rerun the interval search or any benchmark.

The additive `artifact/imt_reproduce.py` entry point now supplies `check`,
`inventory`, and `pack --output <fresh.zip>` for these selected finite results.
The accepted local inventory contains 753 files and 640,779,033 uncompressed
bytes, with zero missing or mismatched files. Eight root receipts contribute
their declared pins; this is not a package of every historical dependency
mentioned inside those inputs. Six packaging tests cover missing roots,
changed files, path escape, conflicting pins, overwrite rejection, and archive
tampering. No production archive was created or published.

The latest regression run passes 66 finite-migration tests, 14 artifact
tests, and six paper evidence tests. Both manuscript integration checks and
the current comparison-table check also pass. The historical artifact
producer retains its original hash. No new encoder benchmark, manuscript
rebuild, or commit was performed during this diagnostic/reproduction step.

## Q1-focused engineering update (2026-09-16)

The user refined the broad figure objective to Q1 sensitivity plus full
certified operating points. The three old RM2Sub plots have now been replaced
by the 130-cell no-constant IMT Q1 grid. The selected BCH-256 engineering curve
shows both its Q1 component and full certificate at each of the five lengths.

| log2 K | Q1 margin | Full margin | Q1 minus full (bits) |
|---:|---:|---:|---:|
| 16 | 43.5928335021 | 41.8183267960 | 1.7745067061 |
| 18 | 50.1891811139 | 50.1890763816 | 0.0001047323 |
| 20 | 50.0622720959 | 50.0620882264 | 0.0001838695 |
| 22 | 48.3904952900 | 48.3904916829 | 0.0000036071 |
| 24 | 46.4572568186 | 46.4572549297 | 0.0000018889 |

Each Q1/full pair uses identical maps and geometry. This small higher-occupancy
loss is established for the selected bounds, not for all diagnostic maps or
against the true failure probability. A final analytic review and public
evidence release remain separate tasks.

The BCH-64/128 parameter slices are a separate study of fixed nested maps.
Their new no-constant IMT Q1 grid covers all 130 geometries. It is explicitly
plotted as Q1, not relabeled as a full margin. `PARAMETER_SLICES.md` records the
refined scope and preserves the earlier full-cover research and diagnostic family.
The old 77/7/46 status counts do not apply to the new maps.

`artifact/imt_reproduce.py figures` regenerates four current inputs; `check`
authenticates the grid and verifies all selected comparisons. Add `--include-q1`
to `inventory` or `pack` for the ten-root, 766-file combined finite/Q1 evidence
set (641,495,560 bytes locally, none missing or mismatched). The default
selected-only package retains its 753-file scope. No data are committed.

Validation of this update: 75 finite-migration tests, 16 paper tests, and
14 artifact tests pass (105 total). Selected finite integration, the 11%
asymptotic binding check, and the comparison-table check pass. The 60-page
PDF was rebuilt and the four plots visually inspected; no overfull boxes,
undefined references, or multiply defined labels were reported. The frozen
`artifact/reproduce.py` hash remains unchanged. No new numerical certificate,
encoder benchmark, production archive, commit, or push was performed.

## Editorial cleanup (2026-09-16)

The rendered manuscript now contains no RM2Sub/MR2Sub names. The prior-inner
timing comparison and historical spectrum-curve detour were removed from the
paper, not from this ledger or the frozen evidence. The small diagnostic outers
are identified as [64,32,12] and [128,64,22]. Section 7 distinguishes K16's
43.593-bit Q1 margin from the additional 1.775-bit higher-occupancy loss, and
states the Q1/full scope without repeating it in every caption.

No selected certificate, numerical grid point, timing, or evidence path changed.
The terminology check allows preserved evidence URLs but rejects the old inner
name in prose or visible link labels. Twenty paper tests and fourteen artifact
tests pass; finite and asymptotic integration checks pass. The 60-page PDF was
rebuilt and the revised figures and performance text visually inspected.

Historical data remain untouched. Generated receipts and raw measurements
are local evidence, not files to add to a source commit. Reproduction must
regenerate or supply their pinned inputs; a source-only checkout must not
claim to contain a complete replay bundle.
