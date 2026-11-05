# SPIN over characteristic-two fields

This workstream starts a new code family over F = GF(2^m), for
m in {8, 32, 64, 128}. It does not replace the binary library profiles.
The completed RS8 library promotion is commit `3ce07a64`.

We want either less encoding work at 10% distance or substantially greater
distance at a comparable rate. Here K and N count field symbols,
R = K/N is the rate, and weight counts nonzero field symbols, not bits.
For the first track, less expansion means a higher R. The active comparison
therefore seeks higher rate at 10% distance and competitive full encoding
cost, rather than another half-rate scalar-extension guarantee.

The initial rate-1/2 high-distance target is close to 50%.

## Current direction: high distance

### Active experiment: cheap lane mixing, 2026-10-16

The current task returns to the 3.8-ms GF256 baseline and seeks higher
record distance without materially more encoding work. The
[cheap-lane experiment](CHEAP_LANE_MIXING_20261016.md) changes only the
state refresh, using fixed companion multiplication plus one random byte.
Matched K20 measurements show approximately 0--3% overhead. The exact
coefficient law and quotient-support lemmas are recorded, but no new
whole-code distance certificate is claimed. The expensive full-field
ladder below is preserved research, not the active implementation target.
The [active goal ledger](CHEAP_LANE_MIXING_GOAL.md) records audited local
bounds, a directed **40% / more-than-76-bit quiet-class bound**, and a
new **40% / more-than-44-bit one-active-group bound** covering every lane
rank. Nonquiet multi-group classes remain open at 40%. The first raw
scalar-shadow operator has a quantified overcounting obstruction; direct
source-rank accounting avoids that loss in the one-group case.
At 35%, directed checks also cover exactly two through six active source
groups, with more than 450, 250, 100, 50, and 45 bits respectively. A separate
dense-band certificate covers all 7584--8192-group classes jointly with
more than 768 bits. The helpers passed independent proof/source review
and 192-/256-bit directed replay. Later rounds below extend this coverage.
The six-group improvement counts source
constraints at activation; it adds no encoding work. Shared-point
extinction ranks and output-zero constraints are the next proof refinements.
The new output-rank operator has independently audited local algebra and
useful local gains at intermediate-support marks. A bounded 256-/512-group
pilot rejects a coarse half-total-rank shortcut by a wide margin; it does
not reject the code. No implementation change accompanies these proof
updates. A new hyperplane-union lemma preserves exact
scalar source support through 256 groups. The filtered scalar pilot now
uses that lemma, actual output floors, and separate quiet/active coefficient
bounds. Exact local tests and independent review passed, but selected
256-group diagnostics remain far too loose. Exact pair-only regional
coefficients and a joint live-state filter reduce the loss but do not
close it; the repeated-marginal diagnostic is not a feasible source class.
A new audited local hyperplane comparison permits exact auxiliary state
supports with an explicit factor at most 8/7 per packet before alpha.
A cheap block-support implementation adds no meaningful gain. Exact
singleton extinction geometry is now integrated: a one-input active
excursion emits 15 or 16 records, rather than the previous floor 12.
The selected diagnostic improves only slightly. Even an oracle deletion
of all extinctions leaves a gap at that fixed parameter point. A full
nonzero RS block moment passes exact tests and independent review, but
is already essentially covered by the existing local caps. A new audited
operator pools intermediate state masks before the fractional moment.
It improves the selected artificial diagnostic by about 5,785 bits;
four-packet pooling adds another 1,438 bits after route averaging.
Sixteen-packet pooling regresses, so that tile-size scan is stopped.
The expression remains far above zero and is not a certificate.
An audited outer-weight marker is valid, but multiplying live caps by
the worst-case input weight makes the selected diagnostic worse.
Coarse marker and tile-size scans are stopped. The new
[quiet-excursion analysis](CHEAP_LANE_EXCURSIONS.md) instead retains outer
constraints across quiet intervals. It closes a restricted consecutive
256-group/44-packet-window event with an exact 579-bit bound, including
all source intervals and window choices. A separate direct routing bound
covers all allowed-position window concentration events for at least
seven groups with 592 bits. The exact interval calculation now covers
sizes 196--256 jointly with 45 bits under the same path restriction.
For arbitrary size-256 supports, a sharper scalarization and an exact
mixture argument retain randomness discarded by the first coefficient
bound. Directed 192-/256-bit verification closes that entire window class
with more than 76 bits, including all 16640 matching profiles. These are
restricted trajectory results, not whole-code certificates. A new
all-binary-path quiet pilot gives negative floating bounds for four actual
size-256 profiles, but the complementary 17-state pilot remains far too
loose. Controlled comparisons now separate degree conditioning, longer
raw pooling, and actual active-block counts. Each helps, but the exact
regional zero/full-state control still has a positive log2 expression
of +11617 at the selected parameters. No new support class is certified.
Quiet-mark recalibration improves this control by only 25 bits. A new
fixed-path A/B isolates approximately 15278 bits lost solely to averaging
regional routes before the fractional power. Retaining route occupancies
reduces the all-zero/full-state diagnostic to +538; it remains positive.
A weighted-Jensen test retains quiet/birth/return suppression and recovers
the entire diagnosed loss on that fixed path. Pooling partial masks inside
live segments costs only 54 bits in that restricted test. A fixed positive
live-state cone then covers all mask histories without averaging away
quiet routing. Exact incoming supports and RS output floors sharpen its
return cap. After calibrating alpha, the two-state operator at z=.84,
alpha=119/1280 gives a floating log2 bound of -218.67 for all trajectories
with consecutive256-group source support; one-step support memory improves
it to -227.36. A simpler independently valid cone now has directed192-/256-bit
verification: more than200 bits for consecutive size256 supports, at35%
distance. Carry-safe coefficient arithmetic and exact source counting now
extend this to every size256 support: directed384-bit log2 union
-93.5574055068573 over all16640 ordered profiles and all trajectories.
The counts sum to C(8192,256); the largest run class has43 runs, so neither
the consecutive nor fully fragmented case alone controls the union.
A second-precision pass reproduces that bound and closes every size255
support with log2 bound-105.2677756749493. Size-dependent proof parameters
now close the entire7--32-group band jointly with at least178 bits and
every128-group support with more than196 bits, at35% distance. The band
uses directed arithmetic and exact aggregation, not interpolation between
selected sizes. A new audited joint missing-group comparison removes the
old source-padding penalty. Its conservative one-group specialization
closes every257-group support with94.287 bits at35%, by directed
algebraic reweighting of the complete profile sum. The paired-region cell
verifier now extends the contiguous coverage through224 groups, with more
than43 bits jointly at35%. The joint source comparison also extends
the dense certificate to7552--8192 groups with more than781 bits jointly.
The intermediate sizes still need coverage. No encoder work was added.
A fixed large-support coefficient bound has substantial marker-calibration
slack: at the selected512-group profile, its tilted singleton mean is306
per region instead of the conditioned480. Recalibration closes that profile
with107.5447 bits at35%, independently replayed in directed256-/384-bit
arithmetic. The finite-chain derivatives and a marker-reuse formula for
covering profile cells are independently audited. The initial complete FLOAT
size512 screen is positive and identifies clustered supports near
66 runs as its dominant classes. Calibrating that profile recovers only
150 bits. Retaining both regional occupancy residues modulo32 removes
2448.7 bits from the calibrated bound, leaving log2+311.4. A bounded
output-mark adjustment gives+304.4. Convexity rules out alpha-only closure
of that fixed residue witness for every admissible alpha. The next exact
coefficient calculation closes this clustered profile with82.3388 bits
at35%, confirmed at256 and384 arithmetic bits with different route scales.
The remaining gap there was coefficient aliasing, not an encoder change.
Retaining both residues modulo128 then gives a complete FLOAT size512
union of2^-59.9352 over all66048 ordered profiles. A directed hybrid of
exact and cyclic coefficient bands now certifies that entire source size
with60.8664 bits at35%, including every state path. It reproduces the
earlier82.3388-bit exact-profile check. Shared coefficient tables now
certify all supports of sizes508--515 jointly with43.5679 bits at35%.
Including size516 lowers that bound to37.0975 bits, below the40-bit
working target. This does not establish a low-weight codeword.
The next broad-range test proves that a route-independent relaxation
cannot close any size512--7525, although the route-averaged size512 bound
already closes. An exact boundary calculation closes the all-singleton
size1024 class with more than4318 bits. Two clustered classes remain
positive after marker fitting and finer coefficient projections. These
results diagnose proof losses; they do not establish low-distance words.
Tightening one-/two-input moments of the actual return matrix recovers
734/797 bits on those selected profiles without changing the encoder.
Exact three-input moments and a two-packet concavity correction now move
their FLOAT log2 bounds to-1068.237/-1431.422. The correction pools a
binary intermediate state without averaging away occupancy information.
Combining the pair bounds with directed boundary coefficients now gives
a complete FLOAT size1024 sum of2^-688.3404 over all263168 profiles and
state paths, without unresolved fallbacks. A new exact-integer regional
verifier with upward dyadic rounding now certifies that entire source
size with **430.0100 bits at35%**, including every profile and state path.
It takes154.5 seconds, generates no files, and leaves the encoder unchanged.
The directed bound is deliberately reported separately from the stronger
FLOAT estimate. The initial round27 complete FLOAT screens tested actual
sizes2048,4096,6144, including an audited hole-degree transformation
above4096; those initial bounds were positive. New affine-coset and
tensor-row rank bounds tighten the actual
return moment without changing the encoder. With calibrated markers and
proof parameters, a selected size2048 profile improves to log2+419.2413.
Retaining first-packet occupancy information in the pair correction then
gives **FLOAT log2−337.4096** at identical parameters, markers, and counts.
Independent proof/source reviews and exact toy controls pass. This is
neither a directed certificate nor a complete-size bound. The complete
size2048 screen retains all1050624 profiles and reproduces that anchor,
but its bound stays positive: ordinary fallbacks dominate, and the separate
resolved-profile envelope is still+546.18. Round28 calibrates31 regional
witnesses and eliminates every fallback. Coherent proof-parameter changes
and finer pair corrections reduce the complete size2048 FLOAT bound from
+686.44 to+243.59,+51.97,and+13.64. An audited two-sided ratio partition
then treats off-diagonal minima exactly, giving **FLOAT log2-29.9703**
over all1050624 profiles. One explicit midpoint proof parameter improves
this to **FLOAT log2-86.7463**, or-86.7425 with the smaller residue256
grid. Both retain every profile without fallbacks. Round29 reconstructs
the fixed witness with exact integer polynomials and directed rounding.
It certifies **86.7425 bits at35% for every size2048 source support**,
including all1050624 ordered profiles and all state paths. Selected
192-/256-bit controls agree; the complete run takes848.8 seconds and
generates no files. No encoder operation or setup distribution changes.
Round30 reuses each regional table for actual sizes2047--2049, with each
size's own coefficient, normalizer, and source factor. The complete
directed band join certifies **74.2104 bits at35% jointly** for all these
supports and trajectories, including1050625 ordered envelope pairs.
It changes no encoding work. Keeping exact run counts with their own
source size before envelope maxima gives another concrete gain: for one
selected profile across sizes2040--2056, the directed log2 bound improves
from+11.1082 to-65.2459. Round31 completes that wider-band replay:
**60.2491 bits at35% jointly for every support size2040--2056**, all
1058840 ordered envelope pairs, and all state paths. The stronger nested
receipts remain available. No encoder work changes.
Round32 completes the directed size4096 replay: **4621.1552 bits at35%**,
all4198400 ordered profiles, exact indexed counts, and all state paths.
It also closes actual size6144 with **49748.3343 bits at35%**, all1050624
ordered profiles, using ordinary-Cauchy bounds in hole coordinates.
The source factor remains that of actual size6144. This cheaper directed
calculation takes32.3 seconds and needs no regional coefficient tables.
The same shortcut initially fails at4096: its positive log2+262.5035
leader is(2013,36), at coarse marker-cell boundaries. Selecting entrywise
among the same31 ordinary witnesses repairs that loss: the complete
directed bound is **256.2334 bits**, taking85.6 seconds instead of941.6
for the stronger coefficient-based bound. No encoder work changes.
Round33 closes **5632--6528 jointly at35%, with4134.3211 bits**, including
all897 integer sizes. The stronger nested5888--6400 receipt has25386.7993
bits. Normalizing marker-cell assignment by the full-pair fraction fixes
a large witness mismatch without changing the markers or encoder.
The extension to5120 fails with that fixed bank; it adds no coverage.
Round34 fits positive rational markers and independently verifies them
with directed arithmetic. Five new disjoint bands join the earlier
5632--6528 and dense bounds: **every size4864--8192 now has a joint35%
bound exceeding780 bits**. The4352--4863 trial remains positive and adds
no coverage. The encoder, refresh distribution, and its work are unchanged.
Round35 extends that continuous cover to **4544--8192**, still jointly
exceeding780 bits. A new audited count of supports within proof cells
also closes **4352--4383 jointly with221.7125 bits**. It improves the
complete3072-group bound substantially, but that bound remains positive.
Neither result adds encoder work. The whole code remains open. The
[goal ledger](CHEAP_LANE_MIXING_GOAL.md#forty-first-round-reuse-exact-tables-across-wider-size-bands)
records current coverage, failed witnesses, and the bounded next calculation.
Round36 retains run-count/conditional-bound correlation and lowers the
complete3072 log2 bound to+1546.85; it is still positive. A matched directed
test then replaces ordinary evaluation by residue256 coefficient projection
for a leading integer profile. That bound improves from+1468.91 to-3379.95
without changing coefficients or encoding work. The subsequent complete
regional replay closes **every3072-group support at35%, with3356.7491 bits**,
including all2362368 ordered matching profiles and state histories.
Round37 replaces target-specific regional coefficients by verified maxima
over their256-by256 residue tables. These reusable bounds close **every
support size3008--3168 jointly at35%, with more than196 bits**, retaining
each size's own denominators, source factor, and indexed run counts.
A13,754-byte selected witness permits18.4-second band replay instead of
the831-second regional-table generation. It does not replace that
witness's verification. The2944 and3264 endpoints remain unresolved with
this fixed bank even when checked alone. Round38 recalibrates the proof
witnesses at3264 and verifies a new regional-bound bank. It closes
3169--3328 and3329--3392, extending continuous coverage to **3008--3392
at35%, jointly with more than196 bits**. The new compact witness is13,785
bytes and replays the first new band in19.6 seconds. Trying all saved
marker pairs at2944 does not improve its completed bound; the3393--3456
extension also remains positive. Whole-code coverage remains open. Next
extend the verifier below its historical641 wrapper limit and test the
smaller-support gaps. The encoder and its work are unchanged.
Round39 implements and checks that lower-domain adapter. Ordinary bounds
at225 and516 remain positive. Combining the current ordered pair caps
with exact downward coefficient extraction closes225--232 jointly with
334.5700 bits and516 with357.8913 bits at35%. Retuned older witnesses and
the new exact helper extend continuous coverage to **1--279 jointly with
more than43 bits**. The full257--289 trial fails; only its fully computed
258--279 contributions are used in this extension. The complete weight516
result covers67080 ordered profiles and all state paths. It improves both
the local operator and coefficient extraction relative to the historical
37-bit result, so no pure extraction attribution is claimed. Whole-code
coverage remains open, with no new encoder operations. The516--548 extension
initially fails with192-bit coefficient rounding. A matched256-bit replay,
with identical local caps and witnesses, closes **all516--548 jointly with
842.2922 bits**. The failed bound suffered common-scale rounding inflation
on rare coefficients, not a demonstrated code weakness. Together with the
older508--515 receipt, the508--548 band exceeds43 bits jointly.
A new z677/825 witness then closes280--312 with198.4058 bits, extending
the continuous small-source cover to **1--312 jointly with more than43
bits**. No size is inferred by interpolation.
Round40 closes **all195 sizes313--507**, joining the earlier ranges, then
closes **549--581 with86.9979 bits jointly**. This adds228 sizes and yields
continuous coverage **1--581 at35%, jointly with more than42 bits**.
Retuning the proof witness repairs a26.9-bit band bound to440.5 bits,
without removing source sizes or changing the encoder. All currently
certified classes, including true quiet histories, have a combined
ideal-setup failure bound below2^-42. Six support-size gaps remain; this
is not a whole-code certificate.
Round41 closes582--614 and615--647, then tests a new129-size exact wrapper.
The648--776 aggregate is positive, but every completed contribution648--751
is below2^-76; that prefix therefore exceeds69 bits jointly. Continuous
coverage reaches **1--751 at35%, jointly with more than42 bits**.
This adds170 sizes without encoder changes. Strong odd/even jumps in the
failed tail motivate a matched numerical-rounding control before further
extension. The whole-code gaps remain listed in the ledger.
The ledger distinguishes these positive majorants from codeword
counterexamples; the implementation and baseline certificate are unchanged.
The unchanged quiet subcode limits this family to at most
49.8046875% distance; that does not obstruct the next 40% checkpoint.

### Leading research design: 35% distance

On 2026-10-09, we selected the [certified paired-route design](HD_SHIFTED_CLOSED_CERTIFICATE.md)
as the leading research design at 35% distance: K=2^20, rate 1/2,
55.092738508223 bits of ideal-setup failure margin, and 3.798160 ms for
precomputed forward encoding of 128-bit records. The timing excludes setup
and allocation. The guarantee concerns the ideal independent setup
distribution, not a particular benchmark seed.

The complete cover was rechecked with 384-bit directed aggregation of the
pinned, previously verified coefficient bounds: no weights or feasible
profiles are missing. This was not a fresh recomputation of every bound.
No production promotion, paper change, commit, or push accompanies this
selection. The previous certified reference and all partial research remain
preserved.

The approximately 3.4-ms candidates still need substantial proof coverage;
we do not have a short, established route to a complete certificate. The
[wrap-up and restart index](HD_34MS_CLOSURE.md#wrap-up-2026-10-09) records
the stopped scans and the separate evidence for each candidate.

### Next direction: larger alphabets near 50% distance

The user requested a new direction after preserving the 35% design. The
initial working assumptions remain K=2^20, rate 1/2, precomputed encoding
of 128-bit application records, one global routing stage, and at least
40 bits of whole-code ideal-setup failure margin. These are search targets;
neither the existing certificate nor its timing transfers to a larger-field
construction.

The earlier [literal GF(2^128) lift and distance ladder](HD_FIELD128_LADDER.md)
targeted a first closure at 35%, then 40%, 45%, 47%, 48%,
and 49%, while measuring performance. This supersedes starting with GF(2^64).

The literal 35% rung has a complete 108.117-bit ideal-setup certificate.
Its measured encoder is about 22.7 ms, so it does not replace the faster
GF256 leading design. A different inner now gives a
[complete 40% certificate](HD_FIELD128_GAMMA30_CLOSED.md) with
49.171469646557 bits of ideal-setup failure margin at K=2^20 and rate 1/2.
It uses a full RS[64,30] emission map, a 30-symbol persistent state,
and a fixed extension-field feedback scalar. The outer and route are unchanged.

The [version-31 index](theory/hd_field128_gamma30_d40_partial_cover_v31_20261010.json)
covers all 8,192 active-group weights with 274 event bounds. The separate
[completion audit](theory/hd_field128_gamma30_d40_completion_v31_20261010.json)
reconstructs this index, checks complete coverage, and independently sums
the actual event ceilings at 384-bit precision. It inherits the reviewed
coefficient producers and saved-data audits; it does not regenerate coefficients.
The guarantee concerns ideal setup and symbol Hamming distance.
The gamma30 encoder now passes dense whole-encoder checks. It initially
measured 30.35--30.44 ms at K20; equivalent deferred outer reduction and
shared completion products now bring it to about 28 ms. The latest matched
completion change saves 1.3%; a route-load fusion passes correctness but
regresses by 2% and is parked. This remains much slower than
the GF256 reference. State-31/32 screens and bounded retuning do not close
the 45% midpoint weight. The 45% and higher rungs remain open; the
[distance ladder](HD_FIELD128_LADDER.md) records the evidence and next work.
The earlier optimistic no-return diagnostic was misleading: its marked
coefficient bounds were too loose. Adding centered markers changes that
selected-profile score from +3561.91 to -7228.94; independent triangular
diagonal bounds give -9872.89. These optimistic scores are not valid
distance bounds. The actual calculation, with returns retained, remains
at +4069.08 with those marks. Subsequent return analysis is recorded below.
The [version-40 receipt](theory/hd_field128_progress_v40_20261010.json)
records this correction, the unsuccessful input-weight tilt, and a
potentially cost-neutral feedback-basis change.
The [translated-feedback proof](HD_FIELD128_TRANSLATED_FEEDBACK_RETURN.md)
now supplies a conditional return cap. Reading the involutive return map
backward sharpens its rank bound without another setup-failure term.
Joint marker retuning makes one formerly failing profile negative, but
the complete 4,198,401-profile check remains open. Subsequent marker and
witness additions reduce its FLOAT log bound from +4012.90 to +3812.69
at active-group weight 4096. The worst profile is now (1613,1727).
This is not a 45% certificate.
The [version-48 receipt](theory/hd_field128_progress_v48_20261010.json)
preserves the earlier diagnostics and reviews. The
[version-53 receipt](theory/hd_field128_progress_v53_20261010.json)
records the monotone full-profile searches and exact local packet provider.
Its eight tests pass, and all 65 occupancies agree with the FLOAT formulas
at two marks. These local checks are not a directed distance certificate.
A [full-quartet coefficient calculation](HD_FIELD128_FULL_QUARTET.md)
then retains all 65 occupancies and improves the complete midpoint score
to +3712.63. Its chronological products precede the fractional power;
every profile and all earlier witnesses remain present.
The [factorial-moment return refinement](HD_FIELD128_RETURN_TAIL_MOMENTS.md)
tightens 39 local occupancy caps, but improves that complete score by only
0.062 bits, to +3712.57. Local gains alone do not resolve the obstruction.
The [wider-state candidates](HD_FIELD128_WIDE_STATE.md) have fresh
principal-invertibility and return bounds, without self-orthogonality or
reversal assumptions. Separate full-profile screens using the seven saved
witness settings give +3885.27, +3628.55, and +3558.11 for states 32, 36,
and 40. No bound transfers between constructions; all three remain open.
Next retune the state-40 witnesses before paying for a still larger state.
The [version-61 receipt](theory/hd_field128_progress_v61_20261010.json)
indexes the new comparisons, 26 passing tests, and remote preservation.
Directed all-weight coverage remains open.
The [version-37 receipt](theory/hd_field128_progress_v37_20261010.json)
preserves the preceding diagnostics and matched timings.
Bulk outputs and source dependencies remain in a verified remote archive,
with only compact receipts stored locally.

Versions 23 and 24 add six independently audited bands covering 129
weights. The latest four bands were computed remotely and cover
140--164, 180--208, 216--240, and 273--303. These changes improve the proof,
not the construction or encoding cost. Version 25 closes 368--383 by
splitting the widest current deficit interval first. Its 802 evaluations
replace an incomplete 2048-evaluation attempt without changing the cell
bounds. Version 26 adds the three complete bands 385--439, 576--639,
and 640--703, plus weight 913. Version 27 adds the independently audited
bands 209--215 and 1160--1191, preserving all 253 earlier events.
Version 28 adds 165--179, 704--757, 968--999, and 1150--1159.
These four bands pass remote independent count, partition, and saved-bound
audits. Verified durable archives retain their trees and helper dependencies;
only compact audit receipts and archive references are stored locally.
Version 29 adds the independently audited band 920--967, preserving all
259 version-28 events. Version 30 adds independently audited weights
548, 914--916, and 1145, preserving all 260 version-29 events.
Version 31 closes the remaining ranges 521--547, 549--575, 917--919,
and 1146--1149, preserving all 265 version-30 events.
The [profile-wise choice of complete bounds](HD_FIELD128_PROFILE_MIN.md)
gives 464.072 bits at weight 548. The newly covered weight 1145 has
49.694 bits; adding this event lowers the covered union margin without
weakening any earlier bound. All new events are archived remotely.
The two composite bands certify 54.411 bits for 521--547 and 51.148 bits
for 549--575. Both passed separate saved-result audits before admission.
Their complete witnesses include their own outer prefactors before taking
cellwise minima; no mixing of entries across witnesses is used.
Retuning the proof parameters and coefficient markers improves the complete
weight-913 FLOAT log bound to -58.946. A fresh numerical replay agrees,
but directed subdivision initially remains incomplete at 4096 cells.
Direct summation subsequently covers all 208,849 profiles and certifies
58.945660269253 bits at that weight. Independent audits check every
integer count and the saved row sum; the complete result is now indexed.
The incomplete tree remains excluded.
Heavy numerical work uses bounded remote scratch and durable remote archives.
Local additions are limited to compact records, with a 20-GiB free-space reserve.
The [calibration record](HD_FIELD128_LADDER.md#closing-lower-ranges)
separates complete directed bands from numerical screens and failed trials.

The preceding weight-5120 calibration produced an 18-marker witness.
Fresh directed coefficients and six independently audited bands closed
4832--5215 without interpolation. A separate
[empty-pair coefficient identity](HD_FIELD128_EMPTY_PAIR_COEFFICIENTS.md)
also closes the formerly problematic boundary profile at unchanged marks;
that single-profile result is not counted as whole-weight coverage.

A sharper support count closes weight 512 with 75.023 bits of margin,
without changing the encoder. Better coefficient markers and retuned
proof parameters also close weight 768 with 250.943 bits. Both receipts
have independent saved-row audits. The new [midpoint-count lemma](HD_FIELD128_MIDPOINT_COUNTS.md)
gives every profile a cheap joint-count bound within a factor of 4097 of
its exact count. It passes independent proof review and exhaustive small
checks. Directed replays strengthen weight 768 to 299.844 bits and pass
weights 511 and 758 at the 54-bit admission target; both are now indexed.
Fresh-marker directed replays also close weights 128, 384, and 1536
at 757.536, 276.747, and 595.312 bits. A new complete weight-3072
directed replay passes at 2274.503 bits. Its independent count and
saved-row audit passes; v15 adds only that exact weight.
The [midpoint rectangle method](HD_FIELD128_MIDPOINT_BOXES.md) also groups
profiles using exact integer prefix sums. Its first directed cover passes
weight 768 in about 1.8 seconds after table preparation. This is proof
computation, not an encoder speedup or additional integer-weight coverage.

The [explicit weight-band extension](HD_FIELD128_MIDPOINT_WEIGHT_BANDS.md)
now reuses one directed coefficient bank across neighboring integer weights.
Every weight is normalized separately; no interpolation is used. Eight
complete, independently audited bands add 294 weights in v16. Their
profiles are grouped by an exact integer partition before taking positive
matrix bounds. Retuning the proof parameters also resolves weight 1000:
its previous vacuous bound becomes a directed 75.126-bit bound, without
changing the code. The adjacent new bands give 74.771 bits for 1001--1022
and 78.631 bits for 1033--1048. The audits independently recompute counts,
check coverage and provenance, and aggregate saved ceilings at 384 bits;
they do not regenerate coefficients. This round changes proof coverage,
not the implementation, setup distribution, or measured encoding time.

Version 17 adds six independently audited bands covering 257 more weights.
They extend the two intervals to 1000--1144 and 3008--3423, and add 4310.
The last 32-weight band has 60.506 bits of combined margin; relaxing its
per-cell stopping threshold completes the proof without weakening the
54-bit band-admission requirement. A separate dense-range FLOAT scan
identifies the no-empty-pair boundary as the weakest profile under the
old parameters. That scan motivated the fresh dense calibration above.
The [ladder record](HD_FIELD128_LADDER.md#rightward-extensions-and-the-dense-boundary)
also preserves the unsuccessful weight-1280 calibration and its resource stop.

The larger alphabet may reduce the total work needed at a fixed distance.
GF(2^64) remains a possible comparison. The existing entropy model gives rate-half
random-code reference distances of 46.8838%, 48.4386%, and 49.2189% for
32-, 64-, and 128-bit alphabets. These are asymptotic comparison points,
not upper bounds or finite SPIN certificates.

Retain the small RS constituents, local superblock mixing, regional route,
and recursive inner as the starting architecture. Fixed maps may retain
GF256 coefficients in a compatible extension-field basis. Some sampled
coefficients must genuinely leave GF256: scalar extension alone preserves
the existing code's symbol distance.

The first proof question is where full-field randomness must act. In the
current recurrence Q_next=D2 H16 D1(Q+CX), widening only D1 and D2 cannot
remove the event Q=0 and CX=0. First test larger-field randomization in the
outer, before the feedback map C, together with the state refresh. Derive
the actual joint cancellation law before using a 1/(2^m-1) probability.
Restricted fast coefficient families need their own bound and cannot be
assigned the uniform full-field guarantee.

The [earlier larger-field screens](HIGH_DISTANCE_INITIAL_RESULTS.md) and
[arithmetic controls](MEDIUM_ALPHABETS.md) remain useful references. Their
expand-then-puncture construction, sampled-class bounds, and multiplier
timings are not a completed fast encoder. Do not resume its eightfold raw
expansion by default. First establish useful sparse and middle-class bounds
for the current SPIN architecture, then compare full encoders at equal
payload size, including required representation conversions. The literal
lift's certificate and timing above supersede this initial planning stage;
higher-distance candidates still require whole-code certificates and timings.

### Parked research history

The preceding proof checkpoint: the column-shear candidate now passes the exact
selected FLOAT test for all560 three-row sets. Its twelve-row directed
class extends to r16..19, and the all-eleven-row r18 class also closes.
Mixed row counts remain open. An optimistic comparison rules out further
middle-arity searches using the current parity-cover relaxation; the next
test uses a shared-scalar cancellation bound instead. Its first two
four-row examples pass the selected FLOAT test without joint eight-variable
enumeration; complete four-row coverage is next. See
[the scoped proof record](HD_LAMBDA_RENYI.md#shear-continuation-broader-coverage-and-a-limit-of-the-cover-bound).
No native timing or whole-code certificate changed in this checkpoint.

The latest [cheap-mixer and wiring test](HD_SCHEDULED_CHEAP_OUTER.md)
uses one fixed field product per column and no added crossed layer.
Changing both input and output row order between columns clears the
complete selected one- through five-row tests. The two-row bound
improves from +69.27 to -165.56 bits when the output wiring is varied.
These are r18/Fmax<=1 FLOAT bounds, not whole-code margins. Checked native
latency is3.442226ms, versus3.354983ms for the matched H* baseline. The
four-row bound is -574.84 bits across all1,820 sets. A fresh directed
comparison supplies63.942051987187 bits for r18 and total RS arity>=247,
including every pairing profile. Middle arities and broader coverage
remain open. All4,368 five-row sets now give a -168.60-bit FLOAT bound;
combining arities1..5 gives -149.04bits for that selected source class.
This permits different arities across groups, but is not whole-code
coverage. Partial-coordinate bounds also pass representative arities6..11.
A new [Renyi comparison](HD_LAMBDA_RENYI.md) now closes the dense-group
part at r18 with40.82172502955 directed bits, across all99 pairing
profiles. The covered class has every group's arity at least12, or
total arity at least246. Exact third/fourth moments compare the randomized
outer to uniform amplitudes before fractional averaging. Scheduled
products keep each original row set fixed; a positive polynomial counts
mixed arities. Independent256/384-bit replays agree. The remaining r18
class has total arity at most245 and some group of arity at most11.
Other group counts also remain open. This proof-only change does not
alter the3.442226ms native implementation.
The latest extension retains at least nine nonzero columns per active
RS group. All group arities11..16 and total arity198..225 now have a
55.4615-bit directed bound, across all99 profiles. A separate marked
input-weight comparison closes homogeneous arity one across those
profiles. Their union with the preceding dense cover has40.821668515604
bits; it is still only a partial r18 certificate. Paired256/384-bit
replays and `theory/hd_lambda_renyi_activity_verify.py` record the coverage.
The marked comparison does not close arities2..5. An ideal-mixer lower
comparator also shows that more column-local mixing alone cannot repair
the current activity-only uniform comparison. Both follow-up tests are
now complete: the explicit whole-column zero reference gives no gain on
the tested grid, and one cheap transvection slightly worsens the combined
bound. Neither changes the measured encoder. The broader structural
screen now finds a two-stage candidate that adds two independent checks
with seven extra XORs per column, while preserving the local two-input
minimum six. A separate staggered column shear improves exact local
zero-column ranks, at an additional cost to the D1 normalization.
Neither has a full certificate or timing. Fresh selected-class FLOAT tests
give -235.73 for all two-row supports with the two-stage mixer and -395.23
with the shear. Seven of eight prescribed middle-support tests pass for
the two-stage mixer. Exact input-weight/output-support enumeration now
changes the failed fixed six-row case from+529.08 to-3010.59, in3.44 seconds.
This is still one selected FLOAT class, not complete six-row coverage.
Separately, the unchanged encoder's directed cover now includes group
arities10..16 and total arity180..202, at56.1554 bits for that added class.
The expanded partial r18 cover has40.821633581232 bits. The linked Renyi
note records exact scopes, receipts, and the failure of row-floor-only
refinement. A reusable exact catalog now checks eight prescribed six-row
sets together at-2984.44 FLOAT bits; all8,008 sets remain to be covered.
Separately, the shear now has a285.081339551997-bit directed certificate
for r18 with twelve active rows per group, all1,820 row sets, and all99
profiles. Exact independent-set witnesses bound the paired-column fourth
moments; fresh256/384-bit replays agree. Neither modified construction
has whole-code coverage or a native timing. Extend mixed arities and
other group counts without transferring evidence between constructions.
An analytic dense-input check shows the current common measure cannot
close arities12..16, even with the RS support floors restored. This is
a limitation of that positive majorant, not a counterexample to the code.
The dense branch needs another regional measure or the unique-preimage
comparison; the proposed brute-force dense run was not launched.

The preceding [extra-stage experiment](HD_OUTER_PAIR_NATIVE.md) adds crossed
pair mixing inside each outer superblock. Exact local laws tolerate
correlated input amplitudes and improve the troublesome quiet-plane
supports. Checked full-encoder medians are3.37ms unchanged,3.85ms with
one layer, and4.35ms with two. These are new uncertified ensembles.
Those costs motivated cheaper first mixers and changed column wiring,
not further stacking onto the complete H* circuit. The local proof and
its remaining whole-code obligations are in [the pair-layer note](HD_OUTER_PAIR_LAYERS.md).

The [cheaper-first-stage proof test](HD_PARITY_OUTER_PAIR.md#state-aware-low-arity-controls)
now has an exact positive input-weight tail, replacing a finite-marker
screen that still counted the zero word. At the selected source class,
the parity mixer clears singleton inputs but not two-row inputs. The
binary rank-five mixer improves the two-row bound but loses singleton
spreading. Neither has a full certificate or a native timing. The
one-product symmetric pair has a checked local majorant, not a multi-layer
support-only law. The complete certified fallback remains unchanged.

An earlier [sparse-input audit](HD_SPARSE_ACTIVITY.md) retained the H*
outer's nine-active-column condition and cleared several failed classes
numerically. Arity five and six remained open for that construction and
bound. Extra marker tuning did not close them. One arity-five term
allowed330 support patterns where only five exist, then discarded
determined output coordinates. This motivated the exact-fiber work above.
Local output shuffles, stronger syndrome checks,
and exposing the refreshed inner state are separate construction options;
none has a new whole-code certificate or measured3.4-ms timing.

The newest [structural proof pass](HD_34MS_CLOSURE.md#a-broader-proof-decomposition-and-a-zero-product-construction-change)
separates rare bad outer setups from the distance first moment. For the
unchanged H* encoder, the global all-quiet low-arity event is bounded by
2^-50.77. A [fixed nonlinear column schedule](HD_BIT_SYNDROME_COLUMN_PERMUTATIONS.md)
also forces several nonquiet columns without adding field products.
Its separate exact setup bound excludes all-quiet words on at most
eight active rows except with probability2^-105.59 across the full code.
These are new proof tools, not a completed distance certificate. The
next step was a sparse-family bound retaining those nonquiet columns.
Its results and remaining gaps are in the audit above. The scalar-state
obstruction was for the static schedule at the selected mark; it does not
by itself rule out every changed schedule or restricted source class.

The active [3.4-ms closure goal](HD_34MS_CLOSURE.md) keeps K20, rate 1/2,
35% distance, and at least 40 bits of failure margin. It explores a
correlation-preserving proof, a cheap pair-refresh modification, and an
exact-map optimization of the certified checkpoint. The prior failed
envelope fits remain parked; the goal is not yet achieved.
The latest [route-predecode experiment](HD_DECODED_ROUTE.md) preserves the
certified map and improves matched full encoding from about 3.81 to 3.68 ms.
It has correctness and sanitizer checks, but does not reach the 3.4-ms target.
In parallel, a [four-XOR outer modification](HD_BIT_SYNDROME_OUTER.md)
replaces two quiet syndromes by five. Its exact local counts improve, but
some pair inputs spread less. A [reviewed column majorant](HD_BIT_SYNDROME_KERNEL.md)
now includes both quiet and nonquiet words. Its [native implementation](HD_BIT_SYNDROME_NATIVE.md)
passes checks and measures about 3.41 ms in the first noisy matched comparison.
The [full-state gate](HD_BIT_SYNDROME_KERNEL.md#refined-selected-class-result)
now clears the formerly troublesome homogeneous r18/arity8/Fmax1 stratum
in FLOAT: -58.17 bits, or -251.46 after a cheap branch-weight recombination.
Both outer branches and all regional boundary states are included, with
fixed row sets across columns. Fresh [directed replays](HD_BIT_SYNDROME_KERNEL.md#directed-selected-class-replay)
at 256 and 384 bits now agree at -246.69 for the larger F<=1 selected class.
Independent 512-bit final checks also pass. The remaining [source/arity coverage](HD_BIT_SYNDROME_COVERAGE.md)
is still required; the complete 35%/40-bit certificate remains open.
No implementation change was needed for this proof step.
The next [mixed-arity gate](HD_BIT_SYNDROME_KERNEL.md#a-first-mixed-arity-extension)
also clears r18/F<=1 when each group has either eight or nine active rows,
including every mixture. One common measure gives a FLOAT bound near
-202.59, independently reproduced. This extension still needs directed
verification. Two [structural proof repairs](HD_BIT_SYNDROME_KERNEL.md#enforce-the-rs-column-activity-before-combining-groups)
now extend the FLOAT gate to every mixture of arities seven, eight, and
nine, near -76.39 bits: enforce nonzero RS columns and compute the exact
empty-step refresh cap. Neither changes the measured encoder. A separate
[rank-eight coupled-quadratic mixer](HD_BIT_SYNDROME_STAR_OUTER.md) has
checked local algebra at ten extra XORs and no extra field products;
it remains untimed and uncertified. The next gate should test broader
coverage with the proof repairs before adding more isolated class proofs.
The latest [whole-word first-moment control](HD_TWO_SUM_LATE_ALPHA_QUIET.md)
numerically clears the previously troublesome quiet, zero-boundary
component at r=18, arity eight, and source F=1 (log2 term -432.03).
A subsequent [fixed-row-set control](HD_TWO_SUM_GF16_MEAN_COLUMN.md) includes all outer moments
while keeping every regional boundary state zero. It still clears the
target numerically (log2 term -88.13) for that same class. Live boundary
paths and complete class coverage remain open. These are FLOAT diagnostics
for the unchanged encoder, not a new distance certificate.
The ordered full-state control initially isolated repeated activation/return
excursions. A matched two-stream control did not improve the bound. The
current [post-refresh coupled inner](HD_POST_REFRESH_COUPLING.md) has an exact
two-hit first-return lemma and a reviewed eight-family weighted comparison,
but its selected-class bound still fails. Its dominant comparison paths
switch between incompatible outer row sets in successive columns. A
whole-word group Minkowski bound now retains those sets, but the full
selected-class bound remains positive. This variant has not been timed or
certified; more witness fitting is not the current priority.

The structural-search pass retains a 3.398-ms experimental G16 candidate,
but its whole-code proof remains open. The earlier envelope searches are
parked; the active goal tests new bounds and a pair-mixer modification. The
[iteration record](HD_STRUCTURAL_ITERATION.md) separates the certified
checkpoints from scoped proof progress. New
[high-arity receipts](HD_TWO_SUM_GF16_HIGH_ARITY_CERTIFICATE.md) close
three selected classes at 62.25 combined bits; the
[complementary low-arity bounds](HD_TWO_SUM_GF16_LOW_ARITY.md) still fail.
An exact [joint input/output support counter](HD_TWO_SUM_INPUT_WEIGHT.md)
retains the stronger RS row-weight constraints without changing the code.
The [empty-separated state bound](HD_TWO_SUM_GF16_EMPTY_CUT.md) removes
a large regional-envelope loss, but the remaining repeated activation
and return terms still prevent a complete certificate.
The latest [source-stratum control](HD_TWO_SUM_GF16_SOURCE_STRATA.md),
[exact empty-transition constraints](HD_TWO_SUM_GF16_STATE_COVER.md), and
[group-activity constraints](HD_TWO_SUM_GF16_GROUP_SEEN.md) do not close
the remaining classes. The latter two use six frozen controls each;
the larger tested classes show no material improvement. Further fitting
of these same envelopes is deferred. Adding more fixed moment terms is
also not a free repair: the [exact local probe](HD_LOWRANK_OUTER_PROBE.md)
finds singular principal restrictions and weaker two-coordinate expansion.
The [larger-Toeplitz circuit review](HD_LARGER_TOEPLITZ_COST.md) retains
the complete proof interface but finds no lower-cost explicit circuit.
Further work needs a new correlation-preserving proof or a demonstrably
cheaper stronger mixer; the verified research winner remains below.

The latest [certified paired-route checkpoint](HD_SHIFTED_CLOSED_CERTIFICATE.md)
encodes K=2^20 128-bit records at rate 1/2 in about 3.80 ms, with 35%
symbol distance and 55.09 bits of ideal-setup failure margin. A matched
three-way comparison gives 4.05 ms for the retained 72.3975-bit reference
and 3.87 ms for its native-output-layout refinement. The complete candidate
therefore saves 6.19% against the reference, or 1.88% beyond layout alone.
The reference remains preserved; neither the library nor paper is changed.

The outer remains four RS[64,32] constituents with independent invertible
Toeplitz column mixers. The changed route alternates adjacent pair matchings
and independently swaps each routed pair. The inner remains h4 with a
16-symbol state. [Exact cyclic counts](HD_SHIFTED_GENERAL_COUNTS.md),
[sparse coefficients](HD_SHIFTED_SWAP_SPARSE.md), common-state potentials,
and directed box bounds now cover every active-group weight from 1 through
8192. Fresh 256/384-bit replays and independent complete accounting agree.
The [certificate record](HD_SHIFTED_CLOSED_CERTIFICATE.md) pins the proof
artifacts, gives the exact event, and records correctness checks and timings.

The first cheaper-outer [structural iteration](HD_STRUCTURAL_ITERATION.md#next-iteration)
targeted the eight-row RS[32,16] XOR-column outer. It removes about
19.7% of outer GFNI products while keeping the group count and region count
unchanged. The first proof gate is its correlated two-constituent class
at 17 active groups under the paired-and-swapped route. Earlier
[XOR multigroup](HD_XOR_MULTIGROUP.md),
[column-potential](HD_XOR_COLUMN_POTENTIAL.md), and
[joint-column](HD_XOR_JOINT_COLUMNS.md) bounds are partial scalar-route
results, not a certificate for this candidate.

The new [paired XOR2 gate and matched comparison](HD_XOR_PAIRED_GATE.md)
records 3.46 ms against a 3.76-ms paired-route control, an 8.13% time reduction.
Its complete 17-group/two-constituent class has a 50.89-bit directed bound:
fresh 256/384-bit replays agree on all 405 cells and the aggregate, with
independent 512-bit union accounting. Whole-code coverage remains open.
A second [exact column-mask gate](HD_XOR_FIXED_MASK_FLOOR.md) closes the
one-group cases with 3 through 8 active constituents at 164.796 bits.
The [mixed-group gate](HD_XOR_MIXED_MASK_FLOOR.md) now gives a directed bound
of 109.86 bits for one concentrated source class: sixteen two-row groups
with a common row pair, plus one denser group. All 219 dense row sets are
covered using three class-dependent witnesses. Fresh 256/384-bit replays
agree on every class ceiling and the union; separate 512-bit accounting
confirms the rounded-cell union. Other mixed classes remain open;
no row set is resampled between columns.
A separate count transfer initially covered every constituent-pair
assignment for that same source class at only 37.75 bits. Two bounded
worst-row-set refinements now raise this to **61.00 directed bits**:
fresh 256/384-bit replays agree on all 219 classes, and independent
512-bit union accounting includes all `28^16` pair assignments. This
remains one source class, not a whole-code certificate.
The adjacent [one-single-row/sixteen-two-row cohort](HD_XOR_SPARSE_MIXED.md)
now has a 70.35-bit directed bound over all 81 source classes and all
pair assignments. Both precision replays and independent rounded-cell sums
agree; the rest of the candidate still needs coverage.
The [coverage roadmap](HD_XOR_COVERAGE_PLAN.md) also gives
a route-independent terminal argument for all cohort mixtures with at
least 61,191 active constituent rows. Its fresh 256/384-bit scalar bound
gives 165.78 bits. Broader sparse and middle-density coverage remains open.
The new [whole-word Hölder study](HD_XOR_DENSE_HOLDER.md) safely compresses
219 fixed dense row sets into twelve witness/arity types. Its one-dense-group
control now retains **56.18 directed bits**, with all twelve ceilings agreeing
at 256/384 bits and a fresh 512-bit rounded-cell sum. The exact-mask control
has 61.00 directed bits. This is a potential multi-dense-group interface,
not a whole-code certificate. Combining all row sets into one shared type,
or only five witness types, loses too much.
The [bounded sparse-cohort refinement](HD_XOR_SPARSE_COHORT_REFINE.md)
repairs the weak two/three-single-row cells with one new witness. At fixed
r=17, preserving the homogeneous and sibling schedules gives a complete
sparse-family proposal of 50.89 bits. The remaining sixteen sparse cohorts
have a floating-point union below 2^-156.78. Their shared directed replay
is still required; the older two-bank failure remains recorded.
The [iteration log](HD_STRUCTURAL_ITERATION.md) records these distinctions.
A separate [two-sum outer](HD_TWO_SUM_OUTER.md) removes the XOR2 equal-pair
exception: every one- or two-input column has at least seven outputs.
An affine coefficient set also makes all singleton and pair support laws
independent of their row labels. [Matched native tests](HD_TWO_SUM_NATIVE.md)
give 3.537 and 3.574 ms for the two variants, versus certified controls near
3.66 ms and XOR2 near 3.47 ms. The affine variant's complete r=17 sparse-family
calculation covers all 81 source geometries and 18 cohorts at about 247 directed
bits, using one fixed witness and an exact hole mixture. Fresh 256/384-bit
replays and independent rounded-cell sums agree. Neither variant has a
whole-code distance certificate. Both remain default-off research options;
existing XOR2 proof work is preserved.
The [sixteen-row affine proposal](HD_TWO_SUM_GF16_PROPOSAL.md) is the next
structural option: shorter RS constituents and 168 fewer field products per
equal-size outer batch, with a checked label-independent sparse law. It is
now implemented. Its first timing regresses, but the
[native diagnosis and repair](HD_TWO_SUM_GF16_NATIVE.md) remove an outlined
store helper and 64 MiB of intermediate stack traffic. A new matched run
gives 3.398 ms versus 3.616 ms for the eight-row affine candidate, a 6.04%
reduction. Its r=17 sparse-family FLOAT gate clears 511 bits; the dense-input
proof and other source weights remain open. This is not a certified winner.
The [closed affine kernel](HD_TWO_SUM_OUTER.md#closed-kernel-for-the-affine-variants)
now removes labelled-mask enumeration from arity-wise dense bounds. Exact
counts depend only on inside/outside support weights; whole-word Holder
weights depend only on total output weight. Independent algebraic review
and exhaustive Group8 rank-kernel comparisons agree. The resulting
[Group16 dense gate](HD_TWO_SUM_GF16_ARITY_DENSE.md) covers all 65,399
fixed dense row sets for one dense group at source `(8,0)`, r=17, with a
bound below 2^-585.54. Its [directed replay](HD_TWO_SUM_GF16_ARITY_CERTIFICATE.md)
now agrees at 256/384 bits, with independent 512/1024-bit union accounting.
A [broader small-cohort certificate](HD_TWO_SUM_GF16_SMALL_CERTIFICATE.md)
covers all source supports for r=2..17 with zero or one dense group, plus
the two sparse r=1 classes, at 563.33 bits. Its source envelope avoids
matching-profile enumeration. A separate
[single-dense-group certificate](HD_TWO_SUM_GF16_SINGLE_DENSE_CERTIFICATE.md)
closes the remaining r=1 class by retaining the inner state through all
columns. Together these receipts cover all source supports for r=1..17
with at most one dense group, at 53.22 bits. The
[two-dense-group screen](HD_TWO_SUM_GF16_DENSE_TWO.md)
covers all 105 arity pairs for source `(8,0)`, r=17, at 476.63 FLOAT bits;
merging those arities prematurely fails and is preserved as a negative control.
Broader dense-cohort and active-group coverage remains open; the
[coverage plan](HD_TWO_SUM_GF16_COVERAGE_PLAN.md) records the missing ranges
and tests of scalable envelopes. These are not whole-code guarantees.
The subsequent [state-preserving envelope](HD_TWO_SUM_GF16_STATE_ENVELOPE.md)
still fails the broad r=18..32 control. Its negative results motivate
[fractional mask bounds](HD_TWO_SUM_GF16_FRACTIONAL_MASKS.md), which retain
the variation between sampled routes. The initial fractional controls and
a bounded live-state-envelope repair also fail. The same note proves a
stronger uniform comparison indexed only by total active constituent rows,
without a row-set counting penalty. Combining that high-row-count bound
with a marked low-row-count bound is the next numerical obligation.
The [Group8 two-dense-group record](HD_TWO_SUM_AFFINE_DENSE.md) also retains
its complete one-source FLOAT gate and a deterministic obstruction to
halving state without halving width. Those are separate results, not
certificates transferred to Group16.
A proposed larger inner opportunity was a
[rank-one refresh with persistent block-mask analysis](HD_RANKONE_STATE_MASKS.md).
Its first three-component proof has a structural obstruction, not a reason
to keep fitting the same bound.
The subsequent mask bounds also fail the tested medium class. A different
[diagonal-XOR-diagonal refresh](HD_XOR_STATE_REFRESH.md) has a checked exact
support law and restores independent nonzero amplitudes after each packet.
An [exact local kernel](HD_XOR_STATE_LOCAL.md) now covers up to two potential
input positions in one RS block. The [native implementation](HD_XOR_STATE_NATIVE.md)
cuts isolated inner time about 20%, but complete encoding improves only
about 1--5%. It has no distance certificate. A proof-preserving
[two-coset pipeline](HD_COSET_PIPELINE.md) instead regresses to 4.03 ms
against a matched 3.78-ms control; cached producer stores are worse.
Both experiments remain disabled. Prioritize the cheaper outer's proof
over a large new inner proof search or another nearby scheduling sweep.

The iteration log also retains the unsuccessful options: a stronger global
inner costs about 4.22 ms, the proof-compatible RS256 substitution costs
about 4.23 ms, and [halving the outer group](HD_SHIFTED_B128.md) fails its
initial middle-weight proof screen. The older 3.30--3.33 ms XOR2 measurements
used a different route; they are not timings of the next matched candidate.

The active [grouped-outer exploration](HD_GROUPED_OUTER.md) returns to the
binary construction's small RS constituents and joint column randomization.
The grouped outer has the same expected-measure envelope as RS256, so the
existing35% and40% certificates transfer with their respective inners.
Four grouped native variants pass independent checks. The first GF65536
representations did not beat the control. The subsequent GF256 and compact
coefficient variants below do improve full-encoder time. Production
profiles remain unchanged.

The [byte-field continuation](HD_GROUPED_SUBFIELD.md) replaces both outer
and inner random scalars by GF256 coefficients. A new inverse-refresh bound
closes its sparse classes: full256/384-bit certificates give72.3975 bits of
margin at35%. Scalar extension retains this distance over GF65536. Its
single-byte-multiply kernels pass K20 oracle checks. The initial g4 timings
were about 4.8 ms, tied with the old control. Compacting the outer keys gives
about 4.49 ms without changing that map. A new independently sampled
[invertible Toeplitz column mixer](HD_GROUPED_OUTER.md#invertible-toeplitz-column-mixers)
has a stronger atom bound and retains the 35% / 72.3975-bit certificate.
Its first compact g4 implementation took about 4.4--4.5 ms. The subsequent
[exact-map inner optimization](HD_GROUPED_OUTER.md#exact-map-inner-optimization)
reduces this to about 4.0 ms: process short cosets, retain the state in the
striped SIMD layout, and compute feedback from the output using CA=0.
The full-output checksum and 35% / 72.3975-bit certificate are unchanged.
The final separate phases are about 1.87 ms outer and 2.04 ms routed inner,
versus 2.54 ms for the preceding inner. Exact toy-field tests, 30,720 packet
basis cases, native K20 dense-oracle checks, and ASAN/UBSAN checks pass.
All timings are serial CPU 2 runs with an unrelated CPU 0 job sharing L3.
The [outer tile sweep](HD_GROUPED_OUTER.md#outer-tile-sweep-no-material-additional-win)
tests eleven sizes from 4 to 8192 groups. The 128-group option is typically
about 1--2% faster, but reversed runs overlap; retain the 64-group default
and the approximately 4.0 ms research timing. Next investigate further
RS-transform/column-mixer fusion to reduce intermediate scratch traffic.

The active [35% optimization](HD35_OPTIMIZATION.md) preserves the 40% winner
and searches for less encoding work at the same rate and length. A four-check
common-refresh inner with RS128 outer has a complete certificate above
64.77 bits and runs around 5.2 ms. Both 256- and 384-bit replays pass.
A cheaper two-check inner with RS256 outer measures about 4.74 ms at K20.
Its complete certificate gives 76.77 bits of margin at 35% distance;
both 256- and 384-bit replays pass. The 35% study does not replace any
of the 40% checkpoints below.

The user has redirected the search to the
[blockwise SPIN restart](HD_BLOCKWISE_RESTART.md): fast local outer,
one routing stage, and a batched stateful inner, without a pre-accumulator.
The target remains K20, rate 1/2, at least 40% symbol distance with a
40-bit setup-failure margin, and full forward encoding below 5 ms.

The proof target is now closed for a
[sixteen-state refinement](HD_SUMRANK_FEEDBACK_CERTIFICATE.md): 40% distance
at rate 1/2 and more than 59.87 bits of setup-failure margin. It uses
sum-rank-separated fixed maps and independent feedback randomization.
The directed replay covers all 8192 active-row counts. This construction
is preserved as the proof reference. The
[encoding optimization](HD_ENCODING_OPTIMIZATION.md) now has a native
[shared-diagonal RS inner](HD_SHARED_DIAGONAL.md) at about 6.1 ms, with a complete
40%-distance certificate above 49.33 bits. The independent-map
[6.3-ms checkpoint](HD_SYSTEMATIC_MDS.md), with 50.06 bits, is preserved.
The prior 7.4-ms layered LRS
version, its [FFT state maps](HD_FFT_STATE_MAPS.md), and its 48.32-bit
certificate remain preserved. The
sub-5-ms target remains open; this is research code, not a library promotion.

A [block-feedback refinement](HD_BLOCK_FEEDBACK.md) now reaches about
5.5 ms with streaming output, retaining a complete 256-bit directed
certificate above 51.04 bits
at 40% distance. It alternates four/eight-check local feedback, retaining
global RS emission. The 104 rational witnesses are pinned; a separate
384-bit replay confirms all 8192 row counts and 51.04261653 bits of margin.
The streaming implementation computes the same map; the prior
5.7--5.8-ms cached-output checkpoint remains preserved. All 45 native
variants pass Release and ASAN/UBSAN checks.
Pure four-check feedback is about 5.6 ms
but remains uncertified. [Cheaper block emission](HD_BLOCK_EMISSION.md)
now has a certified rank-two alternative, described below. A
[common-map refresh](HD_COMMON_REFRESH.md) reaches about 5.3 ms, but its
local proof has not yet yielded a complete distance certificate.
The [return-bound investigation](HD_RETURN_DESIGN_SCREEN.md) records the
corrected grid-rank bound and failed basis/scalar/state-refresh screens.
A rank-two refresh with eight-check feedback measures about5.6 ms with
cached output and 5.4--5.5 ms with streaming output.
An audit corrected its grouping of complementary transition probabilities.
The corrected 256-bit search now certifies all 8192 classes with more than
50.9966 bits at 40% distance. Its 133 rational witnesses are pinned, and a
separate 384-bit replay confirms the complete certificate.
Earlier, pre-correction partial receipts
remain withdrawn. It is a modest alternative, not a sub-5-ms result.
Exact-map route-layout, address-precomputation, and shared outer-load
controls have not yielded another repeatable improvement. Compacting the
per-step setup also gave no repeatable win. The period-four candidate uses
eight-check feedback only once per four packets and takes about 5.45 ms.
Its complete certificate gives more than 56.6591 bits at 40% distance.
All 119 rational witnesses are pinned; a separate 384-bit replay confirms
all 8192 classes. The simpler schedule has not established a substantial
timing advantage over the other two independently replayed alternatives.
The [optimization checkpoint](HD_OPTIMIZED_CHECKPOINT.md) consolidates the
three certified designs, final K20 validation, reproduction commands, and
the remaining limits. No sub-5-ms or production result is claimed.

The current baseline adds [independent GF65536 outer-output scalars](HD_RANDOMIZED_OUTER.md),
restoring the expected-measure comparison used in the binary packet proof.
The checked prototype costs 7.02ms versus 6.53--6.55ms for its identity
control. Directed [class bounds](HD_RANDOMIZED_OUTER_PROOF.md) close
all low-weight words with at least7548 active rows, and sparse/dense quiet
classes. The [exact sparse ladder](HD_RANDOMIZED_OUTER_SPARSE_LADDER.md)
now covers all words with one through five active rows, including nonquiet
histories, with a union below 2^-56.51. The
[quiet extension](HD_RANDOMIZED_OUTER_QUIET_EXACT.md) covers active-row
counts 1 through 962 and 5865 through 8192 using the combined receipts.
The [occupied-tile state moment](HD_RANDOMIZED_OUTER_FLAT_MOMENT.md) and
[composed local refinements](HD_RANDOMIZED_OUTER_LOCAL_REFINEMENTS.md)
improve the proof bounds without changing the encoder. The directed
[joint-input and coefficient-cap receipts](HD_RANDOMIZED_OUTER_JOINT_INPUT.md)
extend quiet coverage to 1508 through 8192. Quiet counts 963 through 1507
and nonquiet counts 6 through 7547 remain open.
The complete certificate for that four-state native construction remains open.
The latest [tile-conditioned proof](HD_RANDOMIZED_OUTER_QUIET_TILE_CONDITIONING.md)
recovers substantial loss from averaging packet layouts before fractional
powers, but its middle witnesses remain unclosed. A
[support-floor experiment](HD_RANDOMIZED_OUTER_QUIET_SUPPORT.md) identifies
why a coarse input-weight penalty does not resolve them. The algebraically
valid [FFT16 alternative](HD_RANDOMIZED_OUTER_FFT16_PROBE.md) also fails the
tested quiet bounds, so it is not selected or implemented. The next proof
step is to retain state histories longer in the composition. A distinct
[sum-rank-separated inner](HD_SUMRANK_INNER.md) has exact local mask ranks
and input/output moments. Its four-state quiet class is certified with
93.99 bits of margin, but its nonquiet classes remain open. The larger
feedback-randomized refinement above is the fully certified design. The
[closure goal log](PROOF_CLOSURE_GOAL.md) separates accepted coverage from
failed-bound diagnostics. All native timings above are
unchanged; this iteration performed proof computations only.
The [randomizer audit](HD_RANDOMIZED_OUTER_OPTIONS.md) identifies which
weaker distributions can preserve these proof statements. No extra local
shuffle is required at the present routed-symbol granularity.

The following fixed-outer work is preserved as a checkpoint. The
[shared eight-point circuit](HD_BLOCKWISE_SHARED_MAPS.md) adds MDS input
mixing while reusing its transforms for state emission and feedback. Its
exact algebra and structural checks pass. The new
[cube circuit](HD_BLOCKWISE_SHARED_FACTOR.md) shares intermediate XOR sums
across feedthrough, emission, and feedback:59 packed GFNI arithmetic calls
per64 outputs, versus82 for the shared FFT version. The
[native comparison](HD_BLOCKWISE_SHARED_NATIVE.md) now rejects this schedule
as a speed improvement:7.65ms versus6.54ms for the same-map FFT, with a
whole-vector rewrite at7.00ms. The working shift64 map removes a specific
[paired cancellation](HD_BLOCKWISE_SHARED_SUPPORT.md) present at shift8.
The [state-budget bound](HD_BLOCKWISE_SHARED_BUDGET.md) rules out excessively
infrequent state updates, and a [sparse quiet-class proof](HD_BLOCKWISE_SHARED_CUTS.md)
bounds every quiet word using at most256 outer rows below2^-2236.70.
An [all-density type bound](HD_BLOCKWISE_SHARED_QUIET_COMPOSE.md) closes the
high-information quiet class below2^-80.87, and the
[joint-image transfer](HD_BLOCKWISE_SHARED_TRANSFER.md) caps whole-subspace
returns under the actual rank-one mixer. Low-information quiet words and
nonquiet histories remain open for that fixed-outer ensemble. The selected
next step is the randomized-outer proof above, rather than further tuning
of the current cube circuit.

The latest [direct-SPIN native control](HD_BLOCKWISE_NATIVE.md) measures
about6.35ms forward at K20: GRS[256,128], one regional route, and a
width16/state4 inner. It has no full distance certificate and is not promoted.
The first width64/state16 control remains frozen at about9.0ms. Its
[fixed maps and local state bound](HD_BLOCKWISE_MAPS.md) are checked; its
[capped global routing bound](HD_BLOCKWISE_GLOBAL.md) handles one complete-row
trap family. Those global numerical receipts do not automatically transfer
to the smaller-packet construction.

The [quiet-subcode bound](HD_BLOCKWISE_QUIET_SUBCODE.md) now rejects the
width64/state4 and width64/state8 controls for every setup, including with
the larger256-symbol outer. Their timings are cost controls, not viable
40%-distance results. Shrinking packet width together with state avoids
this particular obstruction. The [capped support composition](HD_BLOCKWISE_COMPOSITION.md)
closes all one-active-row messages, while many-row certification remains open.
For the selected t16/s4 geometry, the [interval-span bound](HD_BLOCKWISE_INTERVAL_SPANS.md)
also puts one dimension-forced complete-row trap family below 2^-528.39.
It counts routes once, avoiding a demonstrated uncapped-first-moment failure.
Partial-row shortening and general output histories are not covered.

The [accumulator detour is frozen](ACCUMULATOR_DETOUR_CHECKPOINT.md), with
a verified source archive, measured 8.3-ms control, and all proof results
indexed. It has no complete distance certificate and is not the active
implementation direction. The narrative below records that history.

The following entries describe the preserved pre-restart exploration.
The [high-distance speed goal](HIGH_DISTANCE_SPEED_GOAL.md) targets
40--45% distance below5ms at K20, then about3ms and higher distance, against
the16.6ms forward-RS baseline. Its first exact-map forward control is16.8ms;
no new complete distance certificate is claimed. Construction/proof changes
are needed for a substantial advantage.

The next [two-scalar-scan control](HD_STRIPED_NATIVE.md) reached about 12 ms
but is withdrawn: the [grouped-check bound](HD_BAA_PROOF.md) rules out its
40% target for every setup. Stronger [row-order cuts](HD_REGIONAL_THREE_PROOF.md)
also reject the next three-scan proposal, and a [run cut](HD_RATIONAL_POST_OBSTRUCTION.md)
rejects its rational scalar post-inner. The later
[two-state interval cut](HD_TWO_STATE_CUTS.md) rejects the proposed
[two-additive-state construction](HD_ADDITIVE_TWO_STATE.md) with setup-failure
probability greater than99.8%. Its local laws and
[regional proof machinery](HD_TWO_STATE_REGIONAL.md) remain indexed for reuse.
That search next prioritized fixed affine source-row placement. A directed
spectral bound disables the known quarter-row cut, without proving distance.
Its [forward cost checkpoint](HD_AFFINE_NATIVE.md) reaches 8.327--8.360 ms
with four-row affine groups and streaming stores. Grouping changes the
construction; the [bounded cut audit](HD_AFFINE_CUTS.md) finds no retained
counterexample but supplies no distance certificate. The
[exact source/rank interface](HD_AFFINE_PROOF.md) identifies the remaining
all-support proof requirement.
The [component contraction](HD_COMPONENT_RANK.md) now bounds the entire
permitted-singleton term below 2^-84.93, and the
[permitted-pair bound](HD_COMPONENT_PAIRS.md) below 2^-41.11. These are
positive expansion terms, not a full certificate. The
[dense tail](HD_COMPONENT_DENSE.md) now bounds all classes a>=29943 below
2^-322.58, but also proves that rare routes defeat the uncapped actual-code
first moment in the middle range. Exact MDS filtering cannot fix that issue.
The new [routing-condition gate](HD_ROUTING_CONDITION.md) bounds a family
of these rare structural events below 2^-155.23, uniformly over large row
sets and low-weight outputs with at most 64 runs. For that parked candidate,
further proof work would need to control rank-defect probability, not count
every line in rare large kernels.
The [MDS-aware interface](HD_MDS_COMPONENT.md),
[joint component geometry](HD_COMPONENT_GEOMETRY.md), and
[sparse/offset refinements](HD_SPARSE_REFINEMENTS.md) are preserved as
partial tools; none supplies a full certificate.
The shared upper operator also covers a [three-product post-inner](HD_THREE_PRODUCT.md),
which saves one field product but has not been benchmarked.
The [larger vector-state alternative](HD_VECTOR_OBSERVATION_PROOF.md) has
checked local laws but failed tested middle-class bounds. Its
[rank-refinement review](HD_VECTOR_SOURCE_REFINEMENT.md) rules out one
proposed repair and identifies useful transitions for cheaper inners.
The [composable Cauchy rank refinement](HD_CAUCHY_RANK.md) improves its
saved bounds only modestly; it does not close their large deficits.
Further [grouped cut diagnostics](HD_GROUPED_CUT_FOLLOWUP.md) found no new
counterexample, without bounding finite failure probabilities.
Complete high-distance certification and the sub-5-ms target remain open.

The user accepted the [3.2 ms high-rate winner](GF8_HIGH_RATE_WINNER.md).
Its implementation, proof dependencies, and research notes are preserved in
a verified source-only archive. The 83.33% rate, 10% distance, and 57-bit
K20 certificate are unchanged. Further high-rate tuning is parked.

The active search now targets high symbol distance, ideally close to 50%,
starting at rate 1/2 and K20. Keep 128-bit application records, one regional
routing stage, and a 40-bit whole-code setup-failure target. These are search
criteria, not a new certificate. The coefficient alphabet is open again:
wider records alone do not strengthen a GF256 scalar code's distance.
The [first high-distance screens](HIGH_DISTANCE_RESTART.md) identify a
36% GF256 candidate at rate 1/2 and distinguish it from the larger-alphabet
path toward 50%. Complete certification and native timing remain open.

The [larger-field initial results](HIGH_DISTANCE_INITIAL_RESULTS.md) now
identify promising 42% GF(2^16) and 45% GF(2^32) controls, still without
complete certificates. New validated multiplication kernels make the
larger fields cheaper; the proposed encoder's eightfold raw expansion
remains a major cost concern. The subsequent
[lower-expansion exploration](HIGH_DISTANCE_LOWER_EXPANSION.md) prioritizes
a 4K/[256,64] candidate with positive 40% screens, before completing another
certificate. It also proves a genuine window obstruction for the smaller
outer's direct-2K alternative. The 49% middle bound remains open.

## Preserved GF256 proof and encoding cost

The new [speed-first profile](GF8_SPEED_DESIGN.md) uses GRS[64,40], a
four-input/four-state inner, and three retained outputs. It closes K20 at
83.3333% nominal rate, 10% distance, and 57-bit setup-failure margin.
The [latest native comparison](GF8_SPEED_ROUTED.md) selects relocated
diagonals, sixteen-packet scheduling, a fused outer, and a padded scratch
stride: 3.184--3.201 ms at K20. Restoring input before every timed call
gives 3.332--3.361 ms. This is about20% less time than the matched
[eight-packet checkpoint](GF8_SPEED_NATIVE.md), without changing the code
or its certificate. Production is unchanged.

The [next profiling round](GF8_SPEED_PROFILE.md) splits the remaining work
into about1.25ms inner computation,1.17ms routing, and0.84ms outer work.
Sixteen further exact-map schedules did not produce a material improvement;
the3.2ms checkpoint remains selected.

The [rate frontier](GF8_RATE_FRONTIER.md) now reaches nominal rate 37/44
(84.0909%) at 10% symbol distance with a 183-bit whole-code setup-failure
margin at K20. This is 0.0537 percentage points below the GF256 random-code
benchmark; requested K20 rate is 83.8657% after padding. Four complete new
checkpoints cover fixed and periodic output retention. The strongest proof
uses the [exact MDS spectrum](GF8_EXACT_SPECTRUM_DENSE.md) in the dense class.
There is still one regional route and the same eight-input/eight-state maps.
The [native comparison](GF8_RATE_NATIVE.md) selects the 84% profile:
4.807--4.849 ms for full precomputed K20 transposed encoding of 128-bit
records, with 115-bit certified margin. Fixed-six takes 4.903--4.918 ms;
84.09% takes 5.318--5.350 ms for only 0.3265% fewer retained records.
Folded routing, page-aligned column strides, and fused outer loads/stores
preserve the sampled maps. Production and frozen proofs are unchanged.

The [three-quarter-rate goal](GF8_THREE_QUARTER_GOAL.md) now has a complete
[K20 certificate](GF8_PUNCTURED_WHOLE_CERTIFICATE.md): nominal rate3/4,
10% symbol distance, and116-bit whole-code setup-failure margin. A stronger
GRS[128,84] outer feeds a full-Cauchy8/8 inner; keeping seven of every eight
outputs supplies the rate increase. There is still only one regional route.
The [native research encoder](GF8_PUNCTURED_NATIVE.md) measures about
.240/.938/5.57--5.82ms at K16/K18/K20, with13.38% fewer K20 outputs than
the previous cubic field control at roughly comparable encoding cost.
Only K20 has the new whole-code certificate. No production profile changes.

The earlier unpunctured [full-Cauchy maps](GF8_MDS_PACKET.md),
[partial-live return lemma](PARTIAL_LIVE_RANK.md),
[excursion diagnostic](GF8_EXCURSION_DIAGNOSTIC.md), and [joint source
screens](GF8_JOINT_SOURCE.md) remain preserved. Larger states and wider packets alone did not
close the intermediate-support gap; local puncturing changes that tradeoff.

The [native eight-record implementation](GF8_CUBIC_NATIVE.md) now passes
scalar-reference, large-size, alias, and sanitizer checks. Its first packed
schedule regressed; pairing consecutive packets removes an artificial
state dependency and reaches .237 / .94 / 5.71--5.86 ms at K16/K18/K20.
This is a 14--16% exact-map improvement, but the construction has lower
rate than the previous field profile and remains slower than binary SPIN.
Keep both as research controls; no production change.

The user explicitly wants rates well above 2/3. The [larger-outer
screens](GF8_HIGH_RATE_OUTER.md) show that doubling the outer alone does
not reach 3/4 under the present bound. The GF256 random-code reference at
10% distance is about84.1% rate, not a guarantee for this structured family.
The earlier [inner/return investigation](GF8_HIGH_RATE_INNER.md) tests sharper
rank history, pair mixing, and width8/state6 maps. None of those unpunctured
candidates closes3/4; the newly certified punctured code above is different.

The [new eight-record inner](GF8_PACKET_CERTIFICATE.md) closes a cheaper
arithmetic candidate at K20: GRS[128,83], requested rate 64/99, 10% symbol
distance, and 65-bit whole-code setup-failure margin. A repaired pair of
cubic maps gives a proposed 18 rather than 70 GFNI instructions per sixteen
outputs, at the cost of 4.762% more output records. All 46 new tests and the
complete directed replay pass. Its native result above does not improve the
rate frontier; the certificate and existing profiles remain intact.

The latest [fixed-inner rate sweep](GF8_OUTER_RATE_SWEEP.md) closes a new
GRS[128,88] K18 point: 10% distance, 72-bit whole-code margin, requested
rate 32/47, and 2.083% fewer outputs. K20 dimension88 remains unclosed in
the sparse classes, and dimension89 has middle-class gaps. No new native
timing or production change accompanies this certificate. The timings
below remain for dimension87.

The [GF256 follow-up](GF8_EXPLORATION.md) finds a cheaper outer by changing
the source basis, preserving the entire certified full codeword set.
Fusing the source loads brings full-encoder gains to about 15% at both
K16 and K18 (.244/.956 ms), but both new schedules regress at K20.
The 10% distance, requested rates about 2/3, and
63/79/67-bit certificates are unchanged. This is not yet competitive with
binary half-rate encoding, and does not silently replace its interface.

The parallel [rank-history proof work](GF8_RANK_HISTORY_SCREEN.md) tightens
several bounds but leaves the six-state candidate unclosed. It is a
numerical screen with an independently audited transfer, not a new
whole-code certificate. Existing proofs and production remain intact.

### Preceding medium-alphabet comparison

The latest [GF256/GF32 comparison](MEDIUM_ALPHABETS.md) keeps 128-bit
application records and the useful higher-rate/10% target. Two genuine
GF32 arithmetic implementations pass native and sanitizer checks, but
their isolated coefficient pass is about four times slower at K16/K18.
The [paired proof screens](MEDIUM_FIELD_PROOF.md) do not move the sampled
rate/state frontier: larger coefficients barely help the middle classes.
Keep GF256 as the performance baseline and target the middle-class
analysis or source/packet geometry next. No new full encoder or certificate
is claimed, and no production profile changes.

## Preserved higher-rate checkpoint

**Higher-rate family closed at all three sizes:** GRS[128,87] has complete
10% certificates with 63/79/67 bits of setup-failure margin at K16/K18/K20.
Requested rates are 2/3, 2/3, and 128/189. The new
[K16/K18 replays](FIELD128_SIZE_CERTIFICATES.md) sharpen finite conditioning;
the [K20 certificate](FIELD128_K20_CERTIFICATE.md) remains unchanged.

The [validated native prototype](FIELD128_NATIVE.md) takes approximately
0.288/1.128 ms at K16/K18. Latest K20 matched launches range from 5.73 to
7.46 ms, versus about 3.17 ms for binary SPIN. A compact-route option helps
in direct A/B tests, but the process-level variation is unresolved. The
earlier 6.01--6.10 ms checkpoint is preserved. This is a proved rate benefit,
not a demonstrated encoding-speed win. No production profile changes.

Both exact-map FFT inners were measured and rejected, including a new
lane-local layout: fewer field products do not compensate for their other
costs. The dense inner stays selected. All six current variants pass native
and sanitizer checks. The round is complete for the requested status discussion.
[Design screening](CHEAPER_DESIGN_STATUS.md) records reduced-state and
high-distance gaps, plus actual obstructions for specific cheap mixers.

### Preserved larger-outer checkpoint

**Higher-rate certificate closed:** GRS[256,176], one regional routing
stage, and the symmetric width-sixteen/eight-state inner certify 10%
symbol distance at K20 with at least 119 bits of whole-code setup-failure
margin. The padded rate is 11/16; the requested-input rate is 32/47.
The complete sparse/intermediate/dense replay and scalar-extension claim
are indexed in [the certificate](FIELD256_K20_CERTIFICATE.md).

This is a rate improvement, not a new half-rate/10% result. The
[native implementation](FIELD256_NATIVE.md) passes independent correctness
and sanitizer checks. Exact-map optimizations reach about 0.279/1.23/7.35 ms
at K16/K18/K20, but matched binary SPIN remains faster. Its K20 output is
26.5625% shorter; outer and inner costs consume that benefit. The next step
is a cheaper compatible design, not margin inflation. K16/K18 and width
thirty-two remain open.

### Earlier higher-rate attempts and preserved evidence

The four-state GRS[32,22] profile remains parked. Two
replayed K20 routes force outputs below 10% for every assignment of the
inner coefficients. Unequal regional windows leave shortened source
dimensions 136 and 132 after removing padded rows, but their zero-state
exits impose at most 128 equations. See the [proof and compact
witnesses](WINDOW_DIMENSION_SEARCH.md). This is a realized-code failure,
not merely a loose upper bound. The sampled routes do not by themselves
give a theorem about the ideal ensemble's failure probability.

The replacement had to change outer redundancy, state capacity, or
routing geometry. More randomness in the same four-state inner cannot
repair these routes. GRS[64,44] preserves rate 11/16; GRS[33,22] is a
rate-2/3 alternative. Neither was certified or promoted. The [progress
record](ONE_ROUTE_PROGRESS.md) records this search and preserves all
earlier proof results.

**Objective correction:** the half-rate/10% results below reproduce a
guarantee already available by lifting binary SPIN. They are preserved
controls, not the intended improvement. Further work is gated on higher
rate at 10% distance or substantially higher distance at rate 1/2, with
competitive total cost. See [active targets](ACTIVE_TARGETS.md).

An earlier structural result is a [whole-code obstruction](HIGH_DISTANCE_GLOBAL.md):
a full uniform permutation cannot rescue small outer blocks with a causal
inner. The proposed half-rate/30% instance fails with overwhelming probability
at every requested size. Regional routing avoids that particular obstruction.
Higher-rate regional designs are under proof review. The GRS[32,24]
candidate with state mixing only at region boundaries also fails forty
bits when injection labels repeat every sixteen or thirty-two positions.
Sixty-four labels fix the tested pair mechanism, but a larger shortened
subspace still forces failure probability above 2^-38 at K20. This rules
out the four-state GRS[32,24] profile independently of the injection period;
see [the design review](RATE_DESIGN_REVIEW.md). The next design must change
outer redundancy, state capacity, or routing. The larger GRS[256,176]
checkpoint above supplies the first complete higher-rate certificate.

One reusable positive result is a [proved dense class](REGIONAL_GLOBAL_ROUTE.md)
that holds for any independently chosen invertible inner. The
[capacity comparison](REGIONAL_CAPACITY_COMPARISON.md) now favors [32,22]
and [34,24] four-state outers for the remaining sparse proof. Their rates
are 11/16 and 12/17, with exact one-row and broad window-event bounds, but
no complete multirow coverage. The [cost study](PARITY_STATE_COST.md) compares
their fixed circuits without claiming native timings.

The current [dual-MDS repair](DUAL_MDS_PHASE.md) reverses four state
coordinates before emission, with no extra field products. It strengthens
the period-64 emission geometry and retains the existing conservative
one-row bound. The [sparse proof checkpoint](REGIONAL_SPARSE_COVER.md)
records new joint rank lemmas; complete sparse coverage is still open.
The alternative [fresh random packet maps](RANDOM_PACKET_GATE.md) have
an exact local law but no competitive complete certificate. Their
[cheap source envelope](NONEMPTY_SOURCE_ENVELOPE.md) was tested and parked.

The now-parked candidate uses [Cauchy-form packet maps and fresh
state diagonals](DIAGONAL_REFRESH_PROOF.md), still at rate 11/16 with one
regional routing stage. The rank-moment proof now certifies every inverse
input of weight at least .11N. Separate directed one-row and
[two-row bounds](CAUCHY_Q2_TRANSFER.md) give 107/105/103 and 65/72/77 bits at
K16/K18/K20. The sparse class with at least three rows remains open. Its verified
[four-point circuits](CAUCHY_DIAGONAL_PACKET_COST.md) reduce the estimated
full product volume from 181/88 to 169/88 per source record. This is a
changed construction and an arithmetic count, not a measured speedup or
a complete distance certificate. The prior once-per-region version and
its [.19N-and-above proof](INTERMEDIATE_RETURN_CERTIFICATE.md) are preserved.

The latest [sparse rank transfer](SPARSE_MULTIROW_REVIEW.md) retains exact
regional occupancy and uses the [joint return equations](CAUCHY_JOINT_RETURNS.md)
without counting field-valued messages. Its three/four-row screens are
promising, but larger classes remained unresolved. The evaluated
[actual-state-support transfer](MASK_REGIONAL_SCREEN.md) does not improve
those larger bounds. The window obstruction now takes priority over
further proof refinement for this unchanged profile. No replacement
encoder is promoted.

### Preserved half-rate controls

The archived one-permutation work has a single construction certified at all
three sizes: independent
systematic Toeplitz[64,32] blocks over GF16, one 64-region routing stage, and
a width-16/three-symbol-state inner. It certifies at least 10% symbol distance
with setup failure below 2^-84 / 2^-82 / 2^-80 at K16/K18/K20, respectively.
The proof also removes the input-scaling layer.
See the [construction, complete proof, and replay](K20_RANK_CONTINUATION.md).
The first native timings are about 0.278/1.13/6.81 ms at K16/K18/K20 on
128-bit records. The longer outer is expensive. These certificates close
that ensemble's analysis, but reproduce an existing binary rate--distance
guarantee and are not a useful result for the revised objective.

The earlier complete K16 and K18 certificates are preserved. A GRS[16,8]
outer over GF16, one regional routing stage, and a width-16/four-symbol-state
inner achieve rate 1/2 and minimum symbol distance at least 13108 out of
131072, except with probability below 2^-42 over ideal setup. The inner uses
independent uniform GL4(GF16) refreshes and independent nonzero input scalars.
Scalar extension gives the same distance over all four requested fields.
At K18, six state symbols give minimum distance at least 52429 out of
524288, except with probability below 2^-43. The outer and single regional
route are unchanged; see the [K18 proof and replay](GF16_SIZES.md).
See the [complete proof and certificate](GF16_PROOF.md),
the [progress record](ONE_ROUTE_PROGRESS.md), and
[independent audit](SUBFIELD_AUDIT.md).

The packed four-state precomputed transpose takes about 0.177/0.739/5.28 ms
at K16/K18/K20 on 128-bit records. Only its K16 point has a certificate.
The certified six-state K18 profile measures about 1.014 ms using its exact
FFT alternative, a 3--5.5% gain across matched campaigns. K20 remains
unclosed for this smaller outer. See the [matched campaigns](RESULTS.md) for the
distinction between these timings and the earlier 5.05 ms result. Small sizes use
scatter; K20 uses region-sized reordering and non-temporal stores. A separate
GF256/two-state candidate reaches 4.08--4.10 ms at K20 but remains unproved.
The [implementation](../../../spin/experiments/char2_field_spin/README.md)
keeps both candidates and exact-map A/B schedules separate from production.

## Preserved reference constructions

The first complete construction uses small polynomial coefficients, a GRS
outer, and two input-weighted accumulators. It needs no full-width field
multiplication. The [finite proof](SMALL_COEFFICIENTS.md) certifies 10%
distance at K16/K18/K20 with 41.29/45.30/49.30 bits of setup-failure margin.
A larger [16,8] outer certifies 30.5% distance at K20 with 58.87 bits.
These bounds apply separately over GF(2^8), GF(2^32), GF(2^64), and GF(2^128).

The [implementation](../../../spin/experiments/char2_field_spin/README.md)
contains a scalar reference and optimized precomputed GF(2^128) transpose.
The two global permutations remain expensive. An earlier implementation uses
a GRS[16,8] outer, GF256 coefficients, regional routing, and a two-symbol
batched inner. It measures about 0.133/0.575/6.78 ms at K16/K18/K20 with
128-bit records. Its exact field bases and setup distribution are implemented
and tested; it still has no whole-code certificate.

At K16 and 10% distance, outward bounds cover one and two active outer rows
jointly below 2^-43. A separate bound covers the dense tail below 2^-64.
The middle remains open. An exact window argument shows that counting all
bad words cannot close this construction: rare setups contribute large
subspaces of bad words. The next proof must count independent constraints,
as in the BAA companion. See [proof status](SUBFIELD_PROOF.md), the
[independent audit](SUBFIELD_AUDIT.md), and [measurements](RESULTS.md).

The sections below preserve the initial design analysis. The original proposal
to use full-field GL matrices has been superseded by the user's requirement
to avoid large multiplications. Its algebra remains useful as a reference.

## What the companion contributes

The adjacent [field-universality project](../../../../gen-BAA/field-universality/README.md)
contains three useful ingredients:

1. [Scalar-extension invariance](../../../../gen-BAA/field-universality/paper/main.tex):
   a generator over a subfield has the same minimum symbol distance over
   every extension field, for each fixed setup.
2. [Shortened-space and equation-rank accounting](../../../../gen-BAA/field-universality/paper/distance-proof.tex):
   bound the existence of a bad message by counting support branches and
   their independent zero equations, rather than all field-valued messages.
3. [A convolution rank-mask transfer](../../../../gen-BAA/field-universality/expander-lift/CONVOLUTION_RANK_MASK_TRANSFER.md):
   fresh independent coefficients from a set of Q distinct nonzero elements
   cancel r independent forms with probability at most Q^(-r).

These references describe the local checkout inspected on 2026-10-03.
Its HEAD was `c6613d1fd62f568bfbe94f033fcef339a7196d7b`, with additional
uncommitted manuscript and proof work. That commit alone is not a snapshot
of the inspected files.

The rank and pivot arguments are characteristic-independent. The companion's
integer Newton outer and its opposite-pair refinements are not. For example,
the Newton [16,8] map reduced modulo two sends the eighth basis input to
a weight-two output. Use a genuine characteristic-two GRS constituent
with distinct field evaluation points instead.

Its global permutation formulas also cannot be imported unchanged: a SPIN
region route has different support-placement probabilities. Similarly,
the independent-coefficient pivot bound does not apply directly to the
dependent entries of a uniformly sampled invertible matrix.

## Coefficient field versus storage width

Using the current XOR-only encoder on wider records does not improve distance.
For a matrix G over a subfield E of F, expand each F-valued message in an
E-basis. The support of its image is the union of the supports of its
E-valued components. Consequently d_F(G) = d_E(G).

This also gives a useful implementation option. Every requested field
contains GF(256). One generator with coefficients in GF(256) therefore
has exactly the same symbol distance over all four fields. Multiplication
by its coefficients can act bytewise in an explicitly chosen GF(256)-basis.
An arbitrary consumer's polynomial basis need not have that layout; any
required conversion must be included in timings.

To benefit from the larger alphabet rather than only wider records,
some coefficients must leave that small subfield. One possible split uses
fixed GF(256) circuits for the outer and inner maps, with larger-field
coefficients only in the small persistent state. The current native-field
prototypes instead use low-degree polynomial coefficients throughout.
Those coefficients leave GF256 without requiring large products.

For orientation, the q-ary random-linear-code GV reference gives:

| m | Distance at R = 1/2 | Rate at distance 10% |
|---:|---:|---:|
| 8 | 38.0472% | 0.841446 |
| 32 | 46.8838% | 0.885344 |
| 64 | 48.4386% | 0.892672 |
| 128 | 49.2189% | 0.896336 |

These are asymptotic existence benchmarks, not upper bounds or SPIN
certificates. They use R = 1 - H_q(delta); see the
[q-ary GV derivation](https://www.cs.cmu.edu/~venkatg/teaching/codingtheory/notes/notes2.pdf).
The model recomputes the table with a stable entropy formula.

The first less-expansion screens should use R = 3/4 and 4/5 over GF(256).
At larger m, R = 7/8 is another useful target. These are research targets,
not a change to the public parameter API. Lowering the rate also opens a
small-field high-distance option: at 49% distance the GF256 GV reference
rate is about 0.385382. More output may be worthwhile if it removes enough
expensive full-field arithmetic; that tradeoff needs an end-to-end measurement.

## A cheap baseline and its obstruction

Consider a scalar-state inner. Sample a_i,b_i independently from F minus
zero during setup, then fix them for every message. For input X in F^N,
start Q_0 = 0 and set

    Y_i = X_i + Q_i,
    Q_(i+1) = a_i Q_i + b_i X_i.

The nonzero restriction matters: when X_i = 0, a nonzero state cannot die.
Allowing a_i = 0 would introduce a state-erasure probability of 1/q at
every zero-input step.

There is an equivalent distribution requiring no multiplication on the
loop-carried state. Sample independent nonzero scales delta_i and gamma_i:

    Y_i = X_i + delta_i * sum_(j<i) gamma_j X_j.

The unused delta_0 and gamma_(N-1) may be arbitrary. For i >= 1,
delta_i is the prefix product of a_0 through a_(i-1), and
gamma_j = b_j / delta_(j+1). This change of variables is bijective on the
used nonzero scales, so they can be sampled directly and independently.
The transpose uses a scaled suffix sum:

    (I^T Z)_i = Z_i + gamma_i * sum_(j>i) delta_j Z_j.

For each fixed input, its exact two-category transfer is small. Let z mark
output weight, p = 1/(q-1), and f(z) = p + (1-p)z. Rows and columns index
zero state and uniform nonzero state. Zero and nonzero input steps have
transfer matrices

    T_0(z) = [[1, 0], [0, z]],
    T_1(z) = [[0, z], [p*f(z), (1-p)*f(z)]].

Fresh coefficients make the next-state category independent of the output
category in the nonzero/nonzero case. This is a fixed-input enumerator;
it is not by itself a minimum-distance certificate.

More importantly, a small state allows an adversarial message to silence
long zero-input gaps. This limitation does not vanish as q grows.

**Quiet-gap bound.** Fix a k-dimensional outer constituent whose b outputs
route to positions p_1 < ... < p_b in {0,...,N-1}. Set every other message
block to zero. Suppose the inner is causal and homogeneous linear, with
zero initialization, an s-dimensional F-state for s >= 1,
nonzero scalar direct-feedthrough at every position, and
no output flush. Fix all its coefficients and the route.

Let ell_j = p_(j+1)-p_j-1 for j < b and ell_b = N-p_b-1. Order these b gap
lengths from largest to smallest. With r = floor((k-1)/s),

    d_min <= N - p_1 - sum_(j=1)^r ell_(j).

To prove the bound, require the state to be zero at the start of the r
selected gaps. These are at most rs < k homogeneous equations on the local
message. A nonzero message satisfies them all. Its outputs vanish in the
selected gaps and before p_1. Injectivity ensures that its remaining
output is nonzero. Counting those positions proves the inequality.

For a batched inner, apply this argument to wholly inactive steps and
charge the full output width of every remaining step. Do not apply the
symbolwise formula inside a step that mixes several input symbols.

This bound is a rejection test for actual routes. It is one reason to
retain SPIN's region structure and explore several field-state dimensions,
rather than commit to the cheapest scalar-state inner.

## First batched candidate

Take a genuine GRS[b,k] outer over GF(256), interpreted over F, and route
field symbols with a SPIN-style region permutation. Choose t = 8 or 16
symbols per inner step and begin with s = 2 state symbols. Larger s remains
available if the gap bound or certificate requires it.

Fix A in F^(t by s) and C in F^(s by t), initially over GF(256).
Sample independent M_i uniformly from GL(s,F) during setup. Given
X_i in F^t and Q_i in F^s, use

    Y_i = X_i + A Q_i,
    Q_(i+1) = M_i Q_i + C X_i.

The input/output map is unit block triangular. Hence the inner is
invertible even though its state starts at zero and is not flushed.
When Q_i is nonzero, M_i Q_i is exactly uniform among nonzero state
vectors. Zero input cannot extinguish that state.

Choose A as an injective [t,s] GRS map and C as a parity-check matrix for
an [t,t-s] MDS code. Then nonzero states spread across at least t-s+1
output coordinates, and nonzero inputs of weight at most s cannot lie in ker C.
The joint behavior of A, C, and routed inputs still needs analysis.

At s = 2, the state update uses four full-field multiplications per step,
plus additions. The A/C maps and outer have separate costs; four divided
by t is not the total arithmetic cost per output symbol. Test GFNI-oriented
GF(256) circuits and CLMUL-based larger-field state updates only after
fixing their bases and matching their scalar reference maps.

The formal transpose of one step, given adjoints V_i for Y_i and R_(i+1)
for the next state, is

    U_i = V_i + C^T R_(i+1),
    R_i = A^T V_i + M_i^T R_(i+1).

Start R_final = 0 and scan backwards. These are transposes over F, not
binary bit-matrix adjoints under an undeclared representation.

## Proof work to do

For a fixed source support S, let V_S be the shortened message space and
D_S its dimension. A branch imposes output-zero equations with rank r.
Every branch realized by a nonzero message has r <= D_S-1. Thus for tau >= 1,

    1{a bad message from S exists}
      <= tau^(D_S-1) * sum_(retained bad branches) tau^(-r).

For an MDS[b,k] constituent, shortening to j designated coordinates has
dimension max(0,j-(b-k)). This supplies the local support polynomial

    1 + sum_(j=b-k+1)^b binom(b,j) z^j tau^(j-(b-k)).

Retain its exact support information while deriving the region-placement
factor and inner continuation bound. A fixed-message uniform-state law
does not establish a uniform bound after earlier zero equations.

One useful local identity is available for the GL state maps. Let H be a
fixed s-by-D state-form matrix of rank h >= 1 in the current quotient.
Let a be a fixed nonzero row and let M be fresh uniform GL(s,F). Then aM
is uniform over the q^s-1 nonzero rows. Therefore aMH is zero with
probability (q^(s-h)-1)/(q^s-1), and equals each specified nonzero form
in the row space of H with probability q^(s-h)/(q^s-1).
This follows by counting the fibers of the rank-h map v -> vH.
Joint output rows require a joint-rank calculation, not independent use
of these marginal probabilities.

The immediate next experiment is a rank-aware continuation bound for this
batched inner, together with the quiet-gap screen. Only surviving designs
should receive optimized kernels. Evaluate both outer rate and state
dimension before tuning individual instructions.

## Validation and performance contract

The initial model has 14 tests, including exact GF4/GF256 transition
counts, transpose and inverse identities, and the GL(2,GF4) rank law.
The other directories add arithmetic, proof replay, and independent audits.
Run all three suites explicitly; discovery from the root does not recurse
into these non-package directories:

```sh
python -B -m unittest discover -s research/workstreams/char2_field_spin -p "test_*.py" -v
python -B -m unittest discover -s research/workstreams/char2_field_spin/lane_mixing -p "test_*.py" -v
python -B -m unittest discover -s research/workstreams/char2_field_spin/theory -p "test_*.py" -v
```

Keep precomputed transposed encoding as the first performance target.
Report setup and encoding separately, and compare both symbol count and
physical bytes. A 128-bit record contains sixteen GF256 lanes but only
one GF(2^128) symbol. Packing, basis conversion, routing, and output stores
must remain inside end-to-end encode timings when required by the API.
Run benchmarks serially. The initial model's 14 checks have been supplemented
by field arithmetic, finite-certificate, local-transfer, and native-implementation
tests. Current measurements and exact commands are in [RESULTS.md](RESULTS.md)
and the [implementation guide](../../../spin/experiments/char2_field_spin/README.md).
