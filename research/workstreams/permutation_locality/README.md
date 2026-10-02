# Faster permutations for SPIN

Start with the [proof index](PROOF_INDEX.md) for retained constructions,
complete versus partial guarantees, replay entry points, and the local
checksummed preservation snapshot. The [shared-shuffle GF16
experiment](gf16_packets/SHARED_SHUFFLE.md) stopped at a documented plateau,
with q=1–22 covered at the 10% target (restricted margin above 44.12 bits).
The active [pairwise-shuffle follow-on](gf16_packets/PAIRWISE_SHUFFLE.md)
has verified every occupancy q=1–49, with a 48.401-bit restricted margin,
but no whole-code certificate. Occupancies q=50–2048 are not fully covered.
Both are separate from the complete independent-row certificates below.

## Current checkpoint (2026-09-29)

The strongest complete GF16 result is now **distance greater than 10%
with 48.6501 bits of setup-failure margin**, using four inner updates.
The fresh 384-bit full assembly regenerates q=1--48 and replays all
427 dense intervals for q=49--2048, with no unresolved cases. The dense
contribution has 90.9171 bits of margin, so the sparse contribution
controls the complete bound. The dense witnesses previously passed
256-bit checks, and the sparse prefix has separate independent 256-
and 384-bit checks. The two-update design retains its complete >9%
certificate with 46.85 bits, independently replayed at both precisions.
The 10% full replay and exact scope-and-sum audit have passed. The
independently regenerated sparse endpoint matches exactly. Minimum
distance is at least 209716 of N=2097152, except with setup probability
below 2.264e-15. The hill-climb goal is closed; no production default or
paper claim is changed. The unfinished two-update 9.25% assembly and
9.5% search were stopped, preserving their checkpoints and partial logs.

The remaining proof margin is controlled by lightest-support terms in
the q=1 CDF-majorant bound, not the dense range. Support 38 contributes
about 98.46% of that q=1 bound. This is a bound attribution, not an
observed code failure or a claim about the actual BCH spectrum.

A new serial paired benchmark at K=2^20 measures **6.881, 7.184, and
7.484 ms** for two, three, and four updates, respectively. Four updates
cost 8.76% over two. All three share the same routing and GF randomizer.
Regular and sanitizer checks pass. See the
[GF16 record](gf16_packets/README.md#four-update-comparison) for scope,
run order, hashes, and the complete proof record. The next step is to
optimize the certified implementation with its distribution fixed.
All 214 proof-tool tests and fresh normal/sanitizer control checks pass.
The detailed chronology
below records earlier checkpoints, not additional complete claims.

The separate [GF(16) packet randomizer](gf16_packets/README.md) replaces
four-bit lane permutations with independent nonzero field multipliers.
Its two-update version took about 6.86 ms in the first confirmed
implementation and has a **complete 9% distance certificate at K=2^20,
with more than 46.85 bits of setup-failure margin**. The proof retains refreshed-state information,
uses exact GF feedback counts, and couples cancellation to emitted weight.
Exact q=1 placement, fresh sparse covers for occupancies 2--96, and a
complete 116-cell dense cover for 97--2048 were replayed at 256-bit
precision. An independent full 384-bit replay has passed. Earlier
certificates remain available. That gain also retains
the expansion weight of newly activated states.
The replayed sparse and dense contributions have 46.85 and 78.77 bits
of margin, respectively; their sum gives the complete bound above.
An exact regional-placement calculation separately gives 42.12 bits
for all one-active-group messages at 10%. Retaining BCH row parity also
improves the dense comparison and completes the 8.2% result. At 8.5%, a fresh-state density
bound and a variance-sensitive shuffle comparison close all 155 dense
cells for q=97--2048. Its full replay verified every sparse occupancy
except q=4, where the old tilt grid missed its stopping budget. The
8.6% dense partition also closes, with 154 cells; its full replay was
stopped after finding the same q=4 gap. Both logs remain available,
and neither run supplies a whole-code claim. Retaining the expansion weight of newly activated states and
retuning the proof's reference tilt close the 9% dense range in 116 cells;
the fresh full assembly has now completed. The 8.7% dense cover closes in 156
cells. Its sparse q=48 bound closes after expanding the output-tilt grid
and retaining newborn state weights; fresh 256- and 384-bit checks
agree at 175.14 bits. Its first full replay passed every sparse occupancy
except q=4. A fine grid subsequently verifies that occupancy at 9%
with 89.25 bits; the new 9% assembly uses this grid. Retaining regional
packet counts also closes the 9.25% dense range in 172 cells, with
an independent 384-bit aggregate replay giving 76.4285 bits for
q=97--2048. The finer-grid sparse q=96 check also closes, with 50.2158
bits after 247 subdivisions. A fresh full 9.25% assembly was started,
without importing earlier selected-occupancy checks. It was stopped
after the four-update 10% result closed; its partial logs are retained.
The 9.25% partial results do not replace the complete 9% claim above.
The two-update 10% whole-code target remains unresolved. These proof changes do not
change encoder performance.

An actual inner change is now being tested as a separate option: three
updates instead of two. A serial paired benchmark on Peach gives
7.194 ms versus 6.873 ms, about 4.68% extra time. The three-update
regional bound passes two difficult 9.5% intervals with substantial
slack, although the corresponding 10% interval checks still leave gaps.
Its exact q=1 calculation at 10% gives 46.4989 bits, reproduced at 256
and 384 bits. Selected q=2,3,4,8 checks also close at 10%. A full three-update
9.9% dense search used direct regional count masses
to reduce interval losses. Two wider selected intervals pass outward
checks with about 71 and 99 bits; complete coverage remains required.
The all-support check at q=48 passes at 9.9% with 54.31 bits. The first
q=64 search exhausted 200 subdivisions without a certificate. A focused
retry omitted important low tilts. The corrected union-grid run also
finished without a certificate after 400 subdivisions, with floating
log2 sum +32.8086. The q=64 retry with checked joint-return and
lazy-density bounds closes at 117.4341 bits after one subdivision.
An independent 384-bit replay reproduces that margin and verifies
q=96 at 56.0982 bits. The first,
unrefined q=96 check exhausted its budget without a certificate.
An earlier dense handoff is also being tested. At q>=49, three selected
9.9% mean intervals have verified bounds of about 140, 84, and 73 bits.
A full dense search now starts there; if it closes, only q=1--48 will
need the separate sparse proof. That complete q=1--48 prefix has now
been regenerated, with 46.8901 bits of margin. The dense range remains
unresolved, so these results are not yet a whole-code certificate.
An exact small-state sampler test also validates the
one-, two-, and three-update mixture used by the proof.
This is not yet a whole-code three-update certificate; none of the
two-update results is silently reused for it.
Retaining density through lazy updates and using exact three-packet
return counts now closes two troublesome 10% mean intervals with
roughly 220 and 69 bits. Independent 256- and 384-bit checks agree.
These refinements change the proof, not the encoder, and are available
to the resumed 9.9% search. Complete-domain and sparse coverage remain.
They are now optional in the sparse verifier as well, with matching
fresh-regeneration support in the assembler. Three two-update 9.5%
mean intervals now pass outward checks, with bounds of about 887, 498,
and 382 bits, agreeing under independent 256- and 384-bit evaluation.
A full dense cover and a refined q=96 support cover are
running. This could improve the bound without paying for a third update,
but is not yet a whole-code certificate. A floating-search underflow was
fixed without changing the outward verifier; the main three-update
search has resumed from its preserved partition.
The latest diagnostics locate a remaining moderate-activity gap in the
regional comparison. Finer variance parts and the tested output tilts
do not close it. Extending the density allocation also fails there.
A new GF-feedback density bound passes small exact tests but does not
improve the tested actual-size interval, so it remains outside the
certificate path. This motivates retaining mass and density together;
none of these diagnostics is a distance limit or a complete certificate.
The existing constructions below are unchanged controls.

The [two-bit construction](two_bit/FIRST_CLOSURE.md) has a complete proof
at K=2^20: relative distance greater than 9.25%, with setup-failure margin
above 49.11 bits. The [tiled implementation](TWO_BIT_IMPLEMENTATION.md)
preserves that construction and replaces the initial slow direct scatter.
Its confirmed precomputed transpose takes 7.42 ms, compared with 9.65 ms
for the production one-bit encoder on the same machine and workload.
The 9.5% search is paused: sparse coverage reaches 424 active pairs,
but the dense checkpoint still has 199 unresolved cells.
The saved witnesses remain local under ignored `tmp/`; this workstream
tracks source, tests, and proof notes, not generated experiment data.

The [joint experiment](JOINT_TUNING.md) tunes packet size and inner
parameters for end-to-end encoding time. It starts with two- and four-bit packets and
two or three updates, retaining the completed two-bit proof as a fallback.
The earlier four-bit and shared-shuffle work below remains available.
Its partial certificates do not apply to a changed construction automatically.

[The packet-frontier update](PACKET_FRONTIER.md) rejects a slower two-bit
gather schedule and compares two- versus three-update four-bit bounds at
9.25% and 10% distance. The later
[complete four-bit certificate](independent_rows/dense_closure/FIRST_CLOSURE.md)
reaches 0.5% distance with more than 43.744 bits of margin; higher-distance
coverage remains incomplete. That certificate does not transfer to the
new GF(16) construction.

[The shape-potential experiment](independent_rows/SHAPE_POTENTIAL.md)
applies a proof idea from the companion lifting paper without changing
the encoder. It improves selected four-bit bounds by 26--29 bits, but
does not close the middle occupancies. An exact two-epoch return calculation
identifies packet history lost by the current fresh-state envelope; a
small history-aware bound is the next proof experiment.

## Earlier experiments

Two active proof tracks are maintained; see [PROOF_TRACKS.md](PROOF_TRACKS.md).
The shared-shuffle route remains the faster candidate with its existing
partial certificates intact. The separate [independent-row track](independent_rows/README.md)
now has exact averaged outer-support bounds and complete outward covers
for every occupancy from 1 through 58, with a combined margin above
43.74 bits for that restricted message class. See
[LOW_OCCUPANCIES.md](independent_rows/LOW_OCCUPANCIES.md).
It does not replace the shared-route driver or inherit its certificates.

Packet-splitting investigation: [PACKET_SPLIT.md](PACKET_SPLIT.md) tests
independent coordinate permutations for the four adjacent BCH rows. Keeping
the global route cache-line sized and restoring row order locally costs
**6.50 ms versus 5.42 ms** for the shared-permutation two-update baseline
at K = 2^20, about 20% overhead. This changes the ensemble; it does not
inherit the proof coverage below. Production and the paper are unchanged.

For the shared-shuffle route, the two-update replay covers every occupancy from 1 through
14, all ranks and supports, with a conservative combined margin exceeding
41.83 bits. The one-group replay uses 384-bit outward precision; the others
use 192 bits. See [TWO_UPDATES.md](TWO_UPDATES.md). The unrestricted
intermediate and dense occupancies remain open, so this is not a full-code
certificate.

[JOINT_CANCELLATION.md](JOINT_CANCELLATION.md) keeps cancellation and emitted
weight together through three active windows. Independent marginal checks
pass, but the selected 64-group scores improve by only another 3.7 bits.
Retuning the witness has a larger effect: tilt .048 and penalty .95 reduce
the support-160 and support-176 scores to +745.25 and +258.94 bits.
An interior grid improves support 160 further to +555.46; support 176 stays
at +258.94. Removing the penalty (setting it to 1) makes these scores worse.
These selected-point diagnostics still do not close the full two-update proof.

[WINDOW_HISTOGRAM.md](WINDOW_HISTOGRAM.md) replaces a coarse multi-window
moment bound by exact coefficient extraction using the expansion map's 20
window histograms. Outward local coefficients improve, but the selected
64-group scores gain only 15--19 bits and remain open. The complete
fourteen-group replay passes at 192 bits with 152.32 bits of margin.
Production is unchanged.

Latest proof diagnostic: [DUAL_SANDWICH.md](DUAL_SANDWICH.md) records an
exactly checked dual-coordinate reformulation, which did not improve the
current bounds, and a separate counterfactual spectrum sensitivity test.
Binomial-shaped dual counts would close two selected 64-group points but
leave the tested support-160 point open. Those substituted counts are not
certified for BCH; the full two-update proof remains open.

Proof-first update: the two-update inner measures **5.40--5.43 ms**, versus
**5.07--5.08 ms** for one update in the same run. Its complete fourteen-group
cover passes at 192- and 384-bit precision with 152.28 bits of margin. See
[TWO_UPDATES.md](TWO_UPDATES.md). This is a separate ensemble, not an extension
of the one-update coverage below. Neither has a full certificate.

[DENSE_DIAGNOSIS.md](DENSE_DIAGNOSIS.md) records why the current 64-group
screens remain open, the transition sensitivity tests, exact five-/six-window
cancellation counts, and tested total-weight and weighted-fresh refinements.
The elementary F_16 bound is too weak at the intermediate supports; the
tested fresh-state reweighting saves only about 19 bits of the 64-group
diagnostic's roughly 3809-bit gap. No new full-code certificate follows.

[SHORTENING_MOMENTS.md](SHORTENING_MOMENTS.md) gives an exact containment
bound that improves the four-word count by about 21 bits per group at
support 80. It also adds verified positive-polynomial shortening bounds.
These help smaller supports but leave the support-128 diagnostic open.

[FULL_FEEDBACK.md](FULL_FEEDBACK.md) gives exact feedback distributions
through eight active windows. Using their nonzero peaks and expansion
classes reduces the 64-group/support-128 diagnostic from +3809.06 to
+3377.72 bits. This removes real slack but is still far from a certificate.

[WEIGHT_AND_REFRESH.md](WEIGHT_AND_REFRESH.md) records a tested joint
total-weight/all-one-column count and a sensitivity check on the update
count. The first joint-weighting trial does not help the difficult points.
Even the uniform-refresh limit leaves the tested support-128 bound open;
nonzero-state returns and their outer multiplicities remain the priority.

[ODD_COLUMNS.md](ODD_COLUMNS.md) bounds the weight of the XOR of four outer
rows, using pointed bases and shortening. The tested scalar coupling does
not improve the difficult diagnostic. It motivates retaining each group's
complete column-weight histogram and averaging shapes before propagation.
The categorical conditioning identity has exact toy checks; the subsequent
pilot is linked below. The default 14-group certificate was
replayed successfully after the optional weighting changes.

[CATEGORICAL_PILOT.md](CATEGORICAL_PILOT.md) records the implemented
histogram-averaging pilot; it does not improve the balanced intermediate
point. A separate information-set proof certifies dense messages with a
bounded **aggregate histogram-conditioning cost**: 1978--2048 active groups
at 20 bits/group, with over 282 bits of restricted-class margin. This is
not coverage of every message at those occupancies; the atypical complement
remains open.

[EVEN_COLUMNS.md](EVEN_COLUMNS.md) gives a separate restricted dense
certificate for XOR-zero four-row groups with mostly weight-two columns.
Dual-complement shortening tightens its outer count and the general
support counts. With at most eight exceptional columns per active group,
the restricted sum over 2037--2048 groups has over 150 bits of margin.
Neither this result nor the typical-histogram result covers all dense messages.

[DUAL_MOMENTS.md](DUAL_MOMENTS.md) derives positive identities for average
shortening counts in the primal and dual codes. Three exact iterations
improve the support-128 count to 341.78 bits and support-160 to 379.40.
The optional occupancy integration is diagnostic; the full proof remains open.

[FLAG_AND_KERNEL_TESTS.md](FLAG_AND_KERNEL_TESTS.md) records two subsequent
negative tests: averaged extension fibers do not improve the tested
all-one-weighted counts, and the tested weighted dual kernels yield
negligible low-weight gains. Both have exact small-code checks; neither
is enabled by default or supplies new occupancy coverage.

[DUAL_SANDWICH.md](DUAL_SANDWICH.md) applies the audited BCH sandwich LP
to dual shells. Exact residual repair improves the weight-30 cap by
11.71 bits and the four-row support-128 count by a further 4.59 bits.
Numerical infeasibility reports and an unsuccessful bounded SoPlex run
are not treated as proofs. The optional cover integration checks all
accepted bounds rationally; full coverage remains open.

Started 2026-09-23. The leading candidate now uses **one-column bundles** and
takes about **5.08 ms**, or **1.88--1.90x** faster than matched production.
Outward replays cover **all ranks with one through thirteen
active groups**, including every group location, with **48.16 bits** of
combined margin. See [ONE_COLUMN.md](ONE_COLUMN.md) and
[CDF_COVER.md](CDF_COVER.md). Six and nine groups have 192-bit replays;
the other covered occupancies also pass 384-bit replays.
Fourteen or more active groups remain open. The user now accepts roughly
5 ms and prioritizes full proof closure over a strict 2x speedup; see
[PROOF_PRIORITY.md](PROOF_PRIORITY.md). Production and the paper
are unchanged; the candidate does not yet have a full distance certificate.

[MULTI_GROUP.md](MULTI_GROUP.md) records the diagnostics that motivated the
state-memory refinement. Direct support coefficients alone gave only about
32 bits on a truncated three-group cube. Retaining fresh feedback separately
from mature-state density now closes the full three- and four-group support
ranges at 71.90 and 75.43 bits respectively. No full certificate is claimed.

The adaptive extension closes five and six groups at 61.98 and 60.65 bits.
CDF-weighted interval sums close seven through ten groups; combining these
with all-one-column constraints and window averages closes eleven at 73.77
bits. Multi-window averages and joint mixed-support witnesses close twelve
at 57.20 bits in 192- and 384-bit outward runs. The combined upper through
thirteen is below 3.188222e-15. [THIRTEEN_GROUPS.md](THIRTEEN_GROUPS.md)
records the count and cancellation diagnostics that motivated
[MATURE_TAIL.md](MATURE_TAIL.md). Tracking low expansion weights closes
thirteen groups at 75.22 bits in both 192- and 384-bit replays. A further
pair-conditioned tail bound improves a difficult fourteen-group point,
but its full cover has not closed. The revised plan compares a stronger
inner and requires a strategy for the entire remaining occupancy range.

[ALL_ONE_COLUMNS.md](ALL_ONE_COLUMNS.md) removes one source of that slack:
it counts columns where all four grouped rows are nonzero. The sixteen-group
class with every support exactly 80 now passes 192- and 384-bit outward
replays at 53.20 bits. [WINDOW_AVERAGES.md](WINDOW_AVERAGES.md) additionally
verifies the support-84 class at 47.97 bits. Other sixteen-group support
vectors remain uncovered; these points do not extend complete coverage
beyond the full-occupancy results above.

[TWO_COLUMNS.md](TWO_COLUMNS.md) records the preceding two-column family,
which has 47.38-bit coverage for one active group. Each change in bundle
size changes the distribution and requires its own proof replay.

[TWO_GROUPS.md](TWO_GROUPS.md) extends that two-column family's coverage to exactly two active groups
whose row spans both have rank at most two. Their combined contribution has
63.25 bits of margin, including shared-epoch collisions and every pair of
group locations. Higher-rank pairs and larger occupancies remain open for
that preceding family.

[RANK_FOUR_SCREEN.md](RANK_FOUR_SCREEN.md) records the diagnostics leading
to this choice. [PHYSICAL_COLUMNS.md](PHYSICAL_COLUMNS.md) and
[VECTOR_ROUTE.md](VECTOR_ROUTE.md) describe the implementation improvements
that made two-column bundles competitive. The preceding four-column partial
certificate remains in [RANK_THREE.md](RANK_THREE.md).

See [GFNI_BCH.md](GFNI_BCH.md) for the earlier 1.6x canonical-block result,
and [RANDOM_BLOCKS_GFNI.md](RANDOM_BLOCKS_GFNI.md) for the uniform-shuffle
distribution and its preceding 1.48x implementation. See
[CANONICAL_BLOCKS.md](CANONICAL_BLOCKS.md) for the preceding
1.32--1.33x performance candidate, exact canonical-block activation bounds,
and rejected lower-sharing BCH circuits. See
[GROUP_DISTANCE.md](GROUP_DISTANCE.md) for the earlier 16-row candidate's
rank-one output-weight bound, [SHORTENED_SUPPORT.md](SHORTENED_SUPPORT.md) for
the exact counting bound and one-active-group activation result, and
[SHUFFLE_AND_SUPPORT.md](SHUFFLE_AND_SUPPORT.md) for earlier performance
screens and the spectrum-only joint-support bound. See
[WIDE_STREAMING.md](WIDE_STREAMING.md) for the confirmed 16-row result and the
late-activation obstruction, and
[COST_AND_CANCELLATION.md](COST_AND_CANCELLATION.md) for local cancellation checks,
[SECOND_ITERATION.md](SECOND_ITERATION.md) for rectangular
bundles, and [RESULTS.md](RESULTS.md) for the first iteration.
This workstream targets lower full-encoder latency by changing
the permutation and its memory layout. The submitted paper and production
defaults remain unchanged. New timings are experimental; no new full distance
certificate is claimed.

## Objective and baseline

Start with precomputed, half-rate BCH [256,128] SPIN on 128-bit elements,
especially K=2^18 and K=2^20. Hold the outer and IMT parameters fixed initially.
Measure transposed encoding first, then check the forward direction for finalists.
Fresh-code setup and refresh are separate objectives, not substitutes for faster
precomputed encoding. Preserve support for natural lengths in the eventual design;
initial power-of-two experiments need not implement every length.

The production route independently shuffles each outer row's 256 coordinates into
256 regions. Each region independently permutes the row labels. Thus each outer
row contributes exactly one coordinate to each region. The transpose kernel
reverses IMT, routes the values, and applies the transposed BCH circuit.

Re-measure the current public implementation as the control. Historical paper
timings are context, not a matched baseline for new experiments.

## What earlier work already establishes

- `spin/experiments/permutation_bank/TARGET_RESULT.md`: bank-based heuristic
  refresh reached about 10% overhead relative to a precomputed code at K=2^18.
  This was a setup/refresh result, not a speedup over precomputed encoding.
- `spin/experiments/feistel/README.md`: on-demand Feistel evaluation incurred
  dependent table loads. Cheap setup alone did not yield cheap encoding.
- `research/workstreams/inner_design/routing_opt/README.md`: exact-map
  two-level routing and direct gather lost against the existing quarter-rate
  implementation. Avoid repeating these unchanged.
- `research/workstreams/inner_design/redesign/README.md`: confining an entire
  outer row to a local chunk admitted an exact low-distance witness. The short
  outgoing state could be cancelled by a nonzero message from that row.

Those measurements used different implementation revisions and sometimes a
different rate. They motivate experiments but do not predict current timings.

## First design space: globally spread, locally grouped

Keep all 256 regions. Explore small groups of outer rows, with group size
g in {2,4,8,16}; g=1 is the distributional control. For each region, independently
permute the groups and independently permute row order within each group.
Initially retain independent coordinate shuffles for the individual outer rows.

This is a well-defined bijection when g divides the number of outer rows. It
preserves one coordinate per row per region, but correlates positions of rows
in the same group. It does not preserve the existing multi-row routing law.

For a fixed row, its position in each region remains uniform: its group location
and within-group position are uniform. Independent draws across regions, together
with the unchanged row-coordinate shuffle, preserve that row's routing marginal.
This identity alone gives no multi-row certificate.

Grouping row labels is not enough to guarantee contiguous memory accesses.
Independent row-coordinate shuffles generally select different BCH columns.
The existing four-row BCH layout makes equal-column values contiguous, not
arbitrary columns. Test both of the following explicitly:

1. Independent row shuffles with a group-local staging layout. Account for the
   cost of converting that layout into the BCH circuit's input layout.
2. A shared coordinate shuffle within each group, enabling equal-column groups
   and potentially contiguous transfers. Treat this as a stronger change to the
   ensemble, with a greater risk of aligned multi-row supports and cancellations.

The second candidate must pass cancellation diagnostics before expensive tuning.
Neither candidate confines complete rows to a short interval of the recursion;
neither is thereby guaranteed to avoid all locality obstructions.

## Evidence gates

Correctness: enumerate routing addresses; check bijection, inverse consistency,
and one coordinate per row per region. Compare each full encoder to a materialized
reference for the same candidate route. Check the forward/transpose identity.

Performance: benchmark matched controls serially, with identical parameters,
buffer ownership, page policy, compiler, CPU affinity, and timing boundaries.
Measure full encoding, route-only diagnostics, setup, refresh, and retained
memory separately. Route-only timings are explanatory, not additive predictions.
Use fixed-width, unrolled emission paths without encode-time allocation.
Store raw measurements outside version control.

Distance: start with single-row marginals and same-group pairs. Examine support
overlap, activation timing, and cancellation of the outgoing IMT state using
actual BCH messages as well as arbitrary supports. Extend the occupancy analysis
to grouped rows if a candidate survives. A permutation-only statistical screen
does not certify full-code distance or SSD security.

Classification: distinguish exact-map optimizations, changes preserving the
original ensemble, new ensembles with new proofs, and heuristic-only candidates.
Do not attach the existing 40-bit certificate to a changed ensemble by default.

## Goal and first milestone

The active goal is at least 2x faster full precomputed transposed encoding at
K=2^20, retaining 10% relative distance and at least 40 bits of setup-failure
margin through a proof appropriate to the new distribution. Compare against
a matched current baseline, not a historical paper number.

Build an isolated C++ experiment for the grouped families and their materialized
references. Establish whether small groups offer a substantial full-encoder win
at K=2^18 and K=2^20 before undertaking a full proof. An initial engineering
threshold is at least 15% lower latency; this retains promising early candidates
without mistaking that threshold for achievement of the 2x goal.
Reject fast candidates with explicit cancellation witnesses. Only then decide
whether to develop a new proof, retain a heuristic option, or abandon grouping.

The later late-activation obstruction rejects g=4,c=8 and g=32,c=1 even when
their complete zero-state constraint matrices have full rank. Use the current
report's proof filter, not only the initial rank screen.

See [ANALYSIS.md](ANALYSIS.md) for the initial exact feedback checks. The isolated
C++ experiment builds through this directory's CMake project. Modes are `verify`,
`rank`, `census`, `tiled`, `direct`, `line`, and `stream`; arguments are exponent, group size,
shared-shuffle flag, mode, seed, optional timed-call count, and optional column
bundle size. Further kernel and diagnostic modes are listed in the latest report.
`run_screen.sh`
acquires all three shared benchmark locks before running serial measurements.

## Two-bit proof replay

The isolated [two-bit study](two_bit/README.md) pairs independently shuffled
BCH rows while retaining the current IMT maps, two updates, and regional
route. It studies which ideas from the archived complete Riffle g=2 proof
can be reused without substituting its different global shuffle law.
The new complete covers verify every occupancy from 1 through 321 active
pairs, with an aggregate
above 49.11 bits; the dense range remains open. There is no complete
two-bit certificate yet. The four-bit
proof and implementation tracks remain separate and unchanged.
