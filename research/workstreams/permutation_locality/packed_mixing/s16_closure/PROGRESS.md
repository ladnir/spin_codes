# Closing the 16-Bit-State Inner

**Closed (2026-09-30):** the t=64, s=16 construction has a complete
greater-than-10% distance certificate with more than 55.08 bits of
whole-code margin at K=2^20, rate one-half. Fresh 384-bit replay covers
all sparse occupancies and 106 dense intervals; its aggregate margin is
55.089909760705 bits. The checked encoder measures 5.3445045 ms, against
a matched 5.3241165 ms retained S19 control. See the
[closure record](CLOSURE.md) for the construction, exact claim, receipt
hashes, implementation checks, and replay command. The following ledger
preserves the unsuccessful t=128 exploration and subsequent progress.

Goal started 2026-09-30: certify at least 10% relative distance and at least
40 bits of whole-code setup-failure margin at K=2^20, rate one-half. The
goal also includes a complete, checked 128-bit-element transposed encoder
and a serial comparison with the retained 19-state implementation.

The [19-state R4 certificate](../R4_CLOSURE.md) remains the baseline:
greater than 10% distance, more than 52.05 bits of margin, and approximately
5.33 ms with the [exact-map GFNI implementation](../../gfni_inner_REPORT.md).
Neither that proof nor its timing transfers to a new 16-state construction.

## Candidate and Completion Criteria

Keep the BCH[256,128] outer, canonical width-eight GL32 mixers, shared
four-row column shuffle, and independent regional packet shuffles.
Start with t=128 and s=16. The expansion is the declared 16-dimensional
subcode in [s16_maps.py](../s16_maps.py). Investigate its fresh weight-five
feedback map and the separate BCH-based feedback map.

A full 16-by-16 update has the same dense GFNI arithmetic cost regardless
of how setup formed it. Candidate setup distributions are products of
eight or sixteen independent transvections and independent uniform
invertible 16-by-16 maps. The uniform case requires its exact refresh law;
it must not be represented by an arbitrarily large finite update count.
All maps are fixed after setup. Output uses the entering state, and
raw-input feedback is added after the sampled state map; there is no flush.

Completion requires all of the following:

1. One explicit construction, fixed maps, sampling distribution, and cutoff.
2. Fresh authenticated outer-count premises and exact map spectra.
3. Outward bounds covering every sparse occupancy and the entire dense
   comparison domain, with no gaps or unproved interpolation.
4. A fresh whole-code replay whose exact aggregate is at most 2^-40.
5. Independent mathematical and implementation review.
6. Full encoder/reference and adjoint checks, then matched serial timing.

## Quantitative Checkpoints

| Checkpoint | Evidence | Status |
|---|---|---|
| Retained S19 baseline | 5.334629 ms; >10% / >52.05-bit certificate | Preserved |
| S16 expansion | Rank 16, exact minimum weight 48 | Passed |
| Fresh weight-five feedback | Rank 16; every packet rank 4; kernel minimum weight 4 | Passed |
| BCH-based feedback | Rank 16; kernel minimum weight 6 | Passed |
| S16 GFNI update core | Physical-basis checks for R4/R8 products, both orientations | Passed |
| Initial weight-five dense screen | At mean .104, R8 best basic log2 proposal +3868.65 | Does not close |
| Explicit-map kernels and replay | All 45 tests in this directory, including independent rounding/overflow checks | Passed |
| Uniform GL16, weight-five feedback, q=1 | Fresh 192-bit upper bound: 62.86-bit margin at 10% | Sparse endpoint closed |
| Uniform GL16, weight-five feedback, q=2 and q=8 | Fresh support covers: 64.41 and 62.09 bits | Closed separately |
| Full weight-five sparse range q=1..32 | Fresh 384-bit replay: 50.15097 bits at 10%; all support ranges | Closed |
| Comparable refined S19/S16 bound | Exact map-specific capped birth density, with adaptive dyadic Walsh enclosures | Implemented; t128 dense failures retained below |
| Complete S16 encoder | Both feedback maps; R8, R16, uniform GL16; two seeds | Independent checks passed |
| Serial K20 screen | Weight-five 5.21–5.23 ms; BCH16 5.26–5.27 ms; S19 control 5.31–5.33 ms | Small preliminary gain |
| t64/S16 sparse and dense complete coverage | Fresh 384-bit replay; 32 sparse occupancies + 106 dense intervals; whole 55.08991 bits | Closed |

The initial positive log2 proposal is an unsuccessful upper-bound search,
not a low-weight counterexample. It lacks refinements used by the S19
certificate. Do not compare that score directly with the final S19 margin.

### First Uniform-GL Checkpoint

The class-only uniform-GL envelope fails at comparison mean .104: the
fresh outward log2 bounds are +4227.09 for weight-five feedback and
+6066.10 for BCH16 feedback. At mean .032 the corresponding bounds are
-436.57 and +233.13. These are selected-point bounds, not code distances.
The outward checks use 192-bit arithmetic and freshly authenticated BCH
count premises. Increasing the output tilt did not rescue the .104 point.

The first proof improvement is to retain a density bound for the state
born from zero. The existing transvection kernel can allocate this birth
between weight classes and a uniform-density envelope; the initial
uniform-GL implementation used only weight classes. Until this option is
restored, an R8/uniform comparison does not isolate the sampling law.
Separately, fixed invertible changes of basis in the feedback map can
change the activation weights without changing its kernel or rank.

The complete encoder was checked at K=2^14 for both maps and all three
update distributions, using seeds 1 and 17, and at K=2^20 for uniform GL16
with BCH16 feedback. Every check includes explicit expansion-map bases,
each epoch's update and transpose, independent inner references, scalar
GL32/production-BCH transposition, a direct forward BCH matrix, full
adjoint equality, and preservation of the unused suffix. Each timing run
also repeats its complete reference check. Setup is outside timing.

Performance evidence is archived under
`tmp/packed-hill/gfni-s16-performance/records.tar`; the remote build is
`/tmp/spin-bch-tune-qIRYRE`. These initial cross-binary screens comprise one
31-call process per candidate and seed, not the final confirmation campaign.
The mathematical setup distribution and the deterministic benchmark seed
are separate scopes, as for the retained implementation.

### Sparse Closure and Feedback-Basis Checkpoint

Restoring the full retained output-tilt menu closes every sparse occupancy
without changing either map or the class-only sparse kernel. The aggregate
bound has approximately 50.1509666 bits of margin at 256-bit precision.
Occupancy q=32 contributes about 50.1913 bits; the next-weakest occupancy,
q=6, contributes about 55.3306 bits. The earlier q32 screens used a coarse
tilt menu and were unsuccessful searches, not a state-size obstruction.
The retained receipt is `tmp/s16-closure/sparse-uniform-weight5-all-p256.json`.
Fresh 384-bit replay confirms all 32 occupancies and the aggregate margin
50.1509666271 bits. The replay is
`tmp/s16-closure/sparse-uniform-weight5-all-replay-p384.json`, SHA256
`cfb49672dca3e84d236b1f8d5fb9e0e87fbbf866e42b5bdbacfa6e0c3b33eb9e`.
It rebuilds exact maps, authenticated outer counts, and numerical
operators; saved endpoints are not premises. This covers q=1..32 only.
Dense coverage is still required for a whole-code claim.

An independent exact scan of sixteen invertible feedback-basis changes
found one weight-five variant whose single-packet activations never enter
the weight-48 expansion class (the original has 3 of 480 labels). Its
birth-weight census is 87 at weight 56, 263 at 64, 123 at 72, and 7 at 80.
The exact matrix and regression census are recorded by
[basis_feedback.py](basis_feedback.py). The feedback kernel and every
column-subset rank are preserved, but the composed inner changes because
the expansion map stays fixed. No full-distance or timing claim transfers
to that candidate; keep it separate from the sparse closure above.

### Dense Refinement Checkpoint

Adaptive-dyadic Walsh inversion restores map-specific birth-density
bounds. Whole-row capped-density allocations are independently reviewed
and checked against exact small-state distributions. At mean .104 the
best tested regional log2 bound improves from +4227.09 to +2101.46; at
.114 it is +3369.89. Both remain failures. At .032 the bound is -446.68.
The fixed feedback-basis candidate changes the .114 proposal by less than
one bit, so that single-packet improvement is not resolving the dense
obstruction. Tracking all 17 expansion-packet profiles instead of five
total-weight classes improves the IID proposals by only about 7.4 bits at
.104 and 6.3 bits at .114. A large profile-based cover is not justified by
those results. With capped births, BCH16 feedback also fails: regional
outward log2 bounds are +4024.99 at .104 and +5306.78 at .114.

The first point receipts are research diagnostics with context schema 1.
Context schema 2 normalizes spectrum dictionary keys before hashing the
map record, so the metadata identity survives a JSON round trip. This
format correction changes no numerical bound; no failed point receipt is
being used as a certified endpoint.

### Next Parameter Change: Shorter Physical Steps

The recommended next control is the existing selected 64-column,
16-state map in `rate_quarter_bch/inner_calibration/maps/t64_s16_selected.json`
(also generated as `Map64S16`). Independent enumeration confirms rank 16
for both maps, C=A^T, CA=0, rank four for every four-column packet, minimum
expansion weight 16 with 20 such states, and feedback-kernel minimum
weight six with 1,984 such words. Its earlier certificates do not apply
to this outer, rate, or setup distribution.

Two physical steps can be analyzed as one 128-bit macro-step, retaining
the same outer and regional geometry. Let L_a be a valid tilted transfer
matrix for one 64-bit step with a active four-bit packets. Conditional on
j active packets in a 128-bit macro-step, the transfer bound is

    sum_a [binom(16,a) binom(16,j-a) / binom(32,j)] L_a L_(j-a).

The sum runs over feasible a. The binomial coefficient ratio is the exact
probability of the support split between the two halves. Conditional on
that split, each half has its own uniform subset and nonzero labels.
Matrix order preserves state across the midpoint: there is no reset or
flush. Each physical step must sample an independent update. For IID
packet activity p, the macro operator is M_64(p)^2.

For the regional argument, condition on the outer support data and the
macro occupancy counts, not the individual slot choices. The independent
regional shuffle then leaves uniform subsets within macros. Conditional
on packet support, the GL32 maps supply independent uniform nonzero
packet labels. The second half's input law is therefore independent of
the weighted state measure produced by the first half. That entering
measure need not be uniform: the physical operators bound arbitrary
nonzero mass, expansion-weight classes, and per-state density. These
nonnegative bounds remain valid under matrix composition.

Equivalently, if F(x)=sum_a binom(16,a) x^a L_a, the macro generating
matrix is F(x)^2. The 64-macro regional placement is F(x)^128, with the
same binom(2048,r) support normalization. This identity uses chronological
matrix multiplication, not commutativity. Capped-density row choices can
give different numerical IID and conditional envelopes; both remain valid
upper bounds, so their numerical equality is not required.

This explicit wrapper keeps 64 macro-steps per region and 16,384 per code,
but performs 32,768 physical inner steps. It must declare both geometries;
it does not invent a 128-column map or inherit a t128 certificate. The
identity above is independently reviewed. The wrapper is now implemented
and tested in [kernel_t64.py](kernel_t64.py); the numerical status is
recorded below. No earlier quarter-rate certificate is reused.

As a separate diagnostic, adding a second independent uniform map after
feedback forces every nonzero outgoing state to be exactly uniform. The
resulting conservative two-state IID bound at mean .114 still gives
positive log2 proposals: about +5445.9 for the weight-five map and +7313.4
for BCH16, versus -4710.7 for the retained 19-state maps. These use the
same comparison mixture and duals. They are floating proposals for a
different recurrence, not certificates or code counterexamples. They
suggest that merely spending more on mixing is less promising than
increasing state size relative to the number of input bits per step.

At the common diagnostic tilt .2948184039 and effective IID packet
activity 304/747, the relevant entries are:

| Fixed maps | Weighted zero-feedback self-loop Z | Uniform-state emission H | Conservative return H/(2^s-1) |
|---|---:|---:|---:|
| S16 weight-five | 7.88429e-8 | 2.54185e-8 | 3.87862e-13 |
| S16 BCH16 | 8.72318e-8 | 2.54185e-8 | 3.87862e-13 |
| Retained S19 | 5.77838e-8 | 2.54916e-8 | 4.86214e-14 |

The emission moment is almost unchanged; feedback cancellations are not.
When the state is zero and CX=0, extra state mixing cannot activate it.
For comparison, averaging C over all independent uniform binary matrices
gives Z=empty+(T-empty)/2^s, where T is the input's tilted weight moment.
The selected weight-five S16 map is only about 0.136% above that benchmark
at this witness. This is neither a universal lower bound nor proof that
no better C exists. It does explain why another generic basis/map search
is less promising than reducing the physical step length.

### T64/S16 Checkpoint

The selected 64-column map, fresh independent uniform GL16 update per
physical step, and explicit two-step wrapper now have nine kernel tests.
They check exact small-state distributions, chronological composition,
the IID/binomial identity, retained state, and rejection of changed maps
or geometry. Eleven sparse-driver tests check scope, coverage, fresh
replay, and rejection of inherited t128 receipts.

At relative distance 10%, a fresh 384-bit replay covers every sparse
occupancy q=1,...,32 and its full support domain. Its aggregate margin is
55.0961813547553899 bits. The weakest occupancy is q=6 (55.18114 bits);
q=32 has 64.46303 bits. This is the sparse contribution to the final
whole-code certificate below.

Fresh outward dense-point checks at means .032, .104, and .114 give
log2 bounds -518.4363, -6299.1476, and -7738.9438, respectively. Thus the
previous t128 dense-point obstruction does not persist at these points.
The complete dense-domain search finished with 106 accepted intervals
and no unresolved cells. Fresh 384-bit replay subsequently checked them
all; their aggregate margin is 62.938762900079 bits. Together with the
sparse contribution, the whole-code margin is 55.089909760705 bits.

The isolated encoder in [gfni_t64_probe.cpp](../../gfni_t64_probe.cpp)
passes independent forward/transpose references, original-matrix update
basis tests, complete-code adjoints, boundary inputs, and suffix checks.
Checks ran at K=2^14 with two setup seeds and K=2^20. Its explicit maps
and generated circuit are bound to the proof maps by
[implementation_t64.py](implementation_t64.py), including the distinct
C++ and proof hash encodings of the same matrices.

Four matched serial 101-call processes per candidate, alternating order
on Peach CPU15, give median-of-process-medians of 5.3445045 ms for t64/S16
and 5.3241165 ms for the retained t128/S19 code (0.38% difference).
These measure precomputed in-place transposed encoding at K=2^20 over
128-bit elements, excluding setup. The t64 update tables occupy 1 MiB
and its packed state occupies 256 bytes. Exact source/binary hashes and
logs remain in ignored `tmp/packed-hill/gfni-t64-performance/`.

Sparse replay receipt: `tmp/packed-hill/sparse-t64-selected-replay-p384.json`,
SHA256 `e74b02217e6daed3516af0d9b0645f5b09bba353c259598d58383bd3b1feb9e0`.
Its search source is `sparse-t64-selected-screen-extend-p256.json` in the
same directory, SHA256
`86dbbfc72505869e76f00bd1df4036e3b5a137a232131d92b6b7c25b54dbc625`.
Dense point and search receipts are in ignored `tmp/s16-closure/`.
The final whole-code receipt is `tmp/s16-closure/t64-whole-p384/whole.json`,
SHA256 `698e849a142641f6af3885f72a4ae5b2f285f0100a17051640c6a97bc0409c3f`.
Both fresh replays and the exact final sum pass. Independent review
confirms their full scope, sum, source hashes, and relation to the
measured construction.

#### Search Continuation

The serial dense search stopped at an intact checkpoint after 145 visited
cells: 40 accepted intervals and 66 unresolved subtrees. The original
receipt is preserved, and its immutable continuation snapshot is
`tmp/s16-closure/t64-dense-snapshot-p256.json`, SHA256
`d7c90c0d2d67d503350db995dbfc438cf68ae43144b3bb8bac7e5fdf95791142`.
Four interleaved searches in [dense_t64_shards.py](dense_t64_shards.py)
continue those exact subtrees. They write `t64-dense-shard0-p256.json`
through `t64-dense-shard3-p256.json` in the same directory. Each worker
rebuilds its authenticated model; no numerical bound is imported as proof.
All four workers finished successfully, accepting 17, 17, 16, and 16
further intervals without additional subdivisions. The merged cover has
106 leaves and no unresolved cells. An independent audit reconstructed
its exact partition from the snapshot and all four shard assignments.
Its receipt is `tmp/s16-closure/t64-dense-merged-p256.json`, SHA256
`55ebfe88c86b908240100801fdde8eaf4393b371daf3809f09ce9e7a272e9589`.

The helper's `merge` mode checked the exact partition and discarded saved
numerical endpoints. The final [whole_t64.py](whole_t64.py) replay uses
that complete cover and the complete sparse search above, at precision
384, target 40 bits, and four dense workers. Its output directory is
`tmp/s16-closure/t64-whole-p384/`. The sparse replay and all four dense
workers completed successfully. Each worker constructed a fresh model
and rechecked its assigned witnesses. Final aggregation uses the exact sum of
new per-cell and per-occupancy dyadics, rounds upward once, and checks the
strict whole-code target again. Source manifests must remain unchanged.

All 129 tests in this closure directory pass, including exact small-state
kernel tests, scope/coverage tampering tests, real archived implementation
binding, safe aggregate serialization, and a tiny fresh-process replay
test. Independent reviews found no composition, encoder, shard-coverage,
or replay-scope blocker. The fresh numerical replay now supplies the
complete certificate; all 521 source/input hashes match their unchanged
before/after manifests. The construction and matched implementation are
summarized in [the closure record](CLOSURE.md).

## Work Split and Progress Rule

Build map-parameterized dense and sparse kernels while implementing the
full encoder separately. Test stronger refresh before holding R4 fixed.
For each experiment, record the construction, exact scope, best bound,
runtime or coverage, and whether the new result improves the frontier.

If repeated parameter or bound changes stop improving the failing region,
record the obstruction and change the construction or proof strategy.
Unchanged failed screens are not progress. Do not start a large coverage
run before representative difficult points have credible bounds.

All numerical receipts remain under ignored `tmp/packed-hill/`; old
certificates and source records remain intact. No raw experiment data is
part of the proposed code changes.
