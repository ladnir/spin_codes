# Full-alphabet lift and distance ladder

## Objective and preserved reference

The user requested this research goal on 2026-10-09: lift the leading
35%-distance design to GF(2^128), then increase its distance. Distance is
the primary objective. Maintaining or improving full encoding time is a
secondary objective, measured rather than inferred from operation counts.

Preserve the [leading GF256 design](HD_SHIFTED_CLOSED_CERTIFICATE.md):
K=2^20, N=2^21, 128-bit records, 35% symbol distance,
55.092738508223 bits of ideal-setup failure margin, and 3.798160 ms
precomputed forward encoding on Peach CPU 2. Setup and allocation are
excluded from that timing. This project neither promotes a library profile
nor changes the paper.

The first target is a new whole-code certificate at 35% with at least
40 bits of failure margin. The subsequent distance ladder is 40%, 45%,
47%, 48%, and 49%. A passing sampled class is not a completed rung.

## Closed first rung, 2026-10-10

The literal lift now has a complete ideal-setup certificate at 35%
symbol distance, with 108.117061854228 bits of failure margin. It covers
all 8192 possible active-superblock counts. The union has 5854 terms:
one for each weight 1 through 5853, and one terminal event covering
5854 through 8192. No interpolation between sampled weights is used.

The [complete manifest](theory/hd_field128_literal35_closed_manifest_20261010.json)
pins the receipts and reviewed proof programs. Its SHA256 is
`bfdb608e9b1372ffefd5cf3ff06a44568bff72abed462e7726ad245840ab3a49`.
The independently reviewed verifier reconstructs the receipt index,
checks construction consistency and disjoint complete coverage, and
recomputes the probability union at both 256- and 384-bit precision.
It does not rerun every regional coefficient calculation. As a separate
check, a 384-bit native complex-ball replay of the weakest class, weight
197, gives a stronger 112.864923707278-bit bound. That replay is retained
in `theory/hd_field128_r197_d35_acb384_spot_20261010.json`; the complete
manifest conservatively retains the original 108.117-bit union bound.

This guarantee is over ideal independent setup, not the prototype's
pseudorandom seed expansion. The measured literal full encoder is still
about 22.7 ms, versus about 3.6 ms for the matched GF256 control below.
Thus the literal lift is validated, but is not a performance winner.
The gamma30 construction below now closes 40%; higher rungs remain open.

## Closed 40% rung, 2026-10-10

The [gamma30 construction](HD_FIELD128_GAMMA30_CLOSED.md) has a complete
40%-symbol-distance certificate at K=2^20 and rate 1/2, with
49.171469646557 bits of ideal-setup failure margin. Its full RS[64,30]
emission map and fixed extension-field feedback differ from the literal
35% inner. The outer distribution and paired route remain unchanged.

The [version-31 index](theory/hd_field128_gamma30_d40_partial_cover_v31_20261010.json)
covers all 8192 active-group weights using 274 events. The separate
[completion audit](theory/hd_field128_gamma30_d40_completion_v31_20261010.json)
passes: it reconstructs the index, checks the complete nonoverlapping
support partition, counts the terminal suffix once, and independently
recomputes the positive event union at 384-bit precision. The guarantee
is minimum symbol distance at least 838861, over ideal independent setup.
The check inherits reviewed numerical producers and saved-data audits;
it does not replay coefficient generation or certify a benchmark seed.

This closes a proof rung, not a performance target. The gamma30 encoder now
passes independent dense whole-encoder checks at K=8192 and K=16384,
including state continuity, route swaps, aliasing, padding, and full-field
linearity. Its K20 full-encode medians are 30.35--30.44 ms. The matched
literal GF(2^128) control is 22.56 ms; the GF256 control is 3.575 ms.
The [validation and timing receipt](theory/hd_field128_gamma30_encoder_timing_v32_20261010.json)
records the exact scope. The research sampler does not certify its seed.

An equivalent outer implementation reduces each final Toeplitz output
after adding its raw products. It uses the same plan, keys, route, and
inner. Dense same-plan checks pass, and K20 output checksums agree.
Serial matched batches give 28.49--28.71 ms versus 30.47--30.54 ms for the
baseline. The outer alone improves from 16.752 to 14.648 ms. These are
roughly 6% and 12.6% reductions, respectively; component medians are not
additive. The [version-34 progress receipt](theory/hd_field128_progress_v34_20261010.json)
pins the sources and results. This is not yet a competitive full-field lift.

### The next distance rung remains open

Fresh 45% screens retain every matching profile at each selected original
weight, with fresh coefficients and the complete midpoint-count bound.
At weight 4096, retuning the state-30 witness leaves a FLOAT log bound
of +4392.37. With that same witness, states 31 and 32 give +4246.88 and
+4150.91. A bounded 24-trial state-32 search improves the latter to
+4071.24, still vacuous. These upper bounds do not exhibit low-weight words.

State 31 retains a fixed degree-16 feedback scalar. State 32 instead uses
the [random-global-scalar variant](HD_FIELD128_RANDOM_GAMMA32.md).
Its displayed scores are conditional on principal nonsingularity; the
additional setup failure term below 2^-62 must be added once to any final
whole-code union. None of these 45% screens is a directed certificate.
The first state-comparison driver failed while serializing NumPy diagnostics;
version 34 repairs serialization without changing its numerical arithmetic.

Version 35 tests composition before the fractional power across four
packets. It refines 16 selected occupancies at four markers, retaining
every omitted pair through a positive remainder. At the selected profile
(2026,2026), the positive-marker FLOAT log bound improves from +7714.37
to +7666.56. This is not a full 65-occupancy quartet test, an alias replay,
or a certificate; it does not rule out a stronger composition calculation.

The [half-aware return lemma](HD_FIELD128_STATE32_HALF_RETURN.md) removes
some spurious return terms using the fixed map's two Cauchy halves.
Independent source review passes, but its coefficient implementation is
deferred after a more discriminating test. Versions 36 and 37 set every
live-to-zero entry to exactly zero before composition and the fractional
power. This is an optimistic ablation, not a valid probability bound.
With all 79 fixed markers, its positive score is still +7286.57.
The matched alias calculation reproduces +4069.08 before ablation and
gives +3561.91 afterward. All selected alias entries pass the numerical
reliability threshold; arithmetic remains FLOAT, not directed.

Version 40 corrects the interpretation of that diagnostic. The original
markers gave loose coefficients even after increasing the FFT grid.
Nested grids can retain the same distant aliases; stable entries do not
establish accurate coefficients. At the selected profile, the exact
histogram has 2026 empty pairs, 44 single pairs, and 2026 full pairs.
The ordinary multinomial saddle is therefore (44/2026,1).
Adding this saddle, the quiet-weighted saddle, and their neighbors changes
the optimistic no-return score from +3561.91 to -7228.94.

Two independent coefficient bounds give a stronger control. In the
no-return calculation, the regional diagonal entries contain only paths
that remain quiet or remain active. If c=max_j(log2 K_j+128j) and
D=max_j M_j[1,1], their log bounds are respectively
128 alpha c-128 alpha r and 128 alpha log2 D. Every regional route has
total occupancy r, so the histogram coefficient cancels its normalizer.
The live bound also follows from M_j[1,1]<=B*t^64, where
t=z+(1-z)/q. Applying these bounds gives -9872.89 for the optimistic
profile calculation. Neither cap applies unchanged when returns are retained.

The actual full-return score stays at +4069.08 after adding the markers.
Thus the previous conclusion that return refinement had little headroom
was premature. No-return scores are diagnostic, not probability upper
bounds for the code. The [version-40 receipt](theory/hd_field128_progress_v40_20261010.json)
records all controls and their limited scope. Four grid/path tests, five
input-tilt tests, and three saddle-control tests pass remotely.
The simple necessary condition W_in>=33r does not improve this fixed
witness: all seven tested tilts above eta=1 worsen the bound.

A construction candidate changes feedback from A^T to B^T.
Here A interpolates from nodes 128..159, while B uses nodes 64..95;
both evaluate on nodes 0..63. Their shared RS column space gives B^T A=0.
The return matrix AB^T remains symmetric but has dense translation blocks
T_192 and T_224 instead of identity diagonal blocks. Every one-column
return then has all 64 output coordinates nonzero when gamma is outside
GF256. This observation does not prove the needed larger mixed-minor ranks.
The [translated-feedback note](HD_FIELD128_TRANSLATED_FEEDBACK_RETURN.md)
now proves a joint-zero bound and setup analysis. A restriction to one
coordinate half is scaled Cauchy, giving a uniform baseline rank bound.
The return map L=I+gamma AB^T satisfies L^2=I and B^T L=B^T.
Exchanging input and output therefore preserves the constrained solution
space and syndrome. It swaps the active-zero set with the unused set.
Applying the bordered-minor argument in both directions limits the rank
loss by floor(min(c,v,32-h)/16), rather than the forward-only expression.
The input distribution is not changed; all probabilities retain q^j as
their original denominator. Principal failure is below 2^-64 once per setup.

The exact GF256 probe checks every entry and all 3,812,256 disjoint
2-by-2 minors. The new map has 4,480 singular minors, versus 2,336,160
for the old one. These small-minor counts are diagnostics, not larger-rank
assumptions. Independent proof/source review and tiny-field rank tests pass.

At the original profile (2026,2026), the forward-only refinement changes
the FLOAT score from +4069.08 to +4057.94; reversal gives +4009.08.
Retuning z and alpha gives +2987.07. Joint coefficient-marker retuning
then gives -3198.45 for that profile, with actual returns retained.
However, replaying every profile with that last witness gives +32091.67:
its worst profile is (1933,1933). Taking the minimum of three complete
per-profile bounds after each transfer and outer prefix gives +4012.90,
with worst profile (1382,1382). Every one of the 4,198,401 profiles is
retained. Selected-profile success does not close the midpoint weight.

The [version-48 receipt](theory/hd_field128_progress_v48_20261010.json)
records these conditional FLOAT results; none is a directed certificate.
The tested additional-gamma-exclusion counting source is preserved, but
its numerical screen was not run: reversal removes those terms without
paying that extra setup risk. The next step is to optimize the complete
profile envelope while retaining earlier witnesses. At occupations whose
return upper cap already meets its true lower floor, finer local ranks
cannot materially help at the same z; coefficient control and composition
before the fractional power remain separate directions.
Source inspection changes one feedback FFT offset without adding operations.
No encoder change or measured performance claim accompanies the proposal.

The next searches retain every earlier bound. Within a fixed `(z, alpha)`
witness, six coefficient-marker additions reduce the complete profile sum
from +4012.90 to +3896.21. Each completed round replays all 4,198,401 ordered
profiles. Adding a separately recomputed witness gives +3882.89, then
single-marker retuning stalls at the asymmetric profile (1258,1613).
Allowing separate even/odd coefficient markers gives three further complete
improvements, ending at +3812.69. There are seven retained witnesses; their
minimum is taken after each complete transfer and its own outer factor.
The worst remaining profile is (1613,1727), whose term is +3809.9550.
No previous numerical bound is discarded, and incomplete rounds cannot
alter the recorded envelope. These are conditional FLOAT results at one
active-group weight, not a 45% distance certificate.

An exact rational packet provider implements the translated return cap,
the original input normalization, and the unchanged quiet/birth/live bounds.
It keeps the global principal-failure probability separate. Eight tests
pass; an additional comparison covers all 65 occupancies at two marks,
including the saved numerical witness. The 520 packet entries agree with
the FLOAT formulas to at most 2.28e-12 in their finite base-two logarithms.
No directed coefficient bank or all-weight union has been produced for this
construction.

At the remaining worst profile, a four-packet diagnostic separates
coefficient extraction from delayed fractional powering. Four markers
retain the current alias score +3809.9550. Their positive-point pair bound
is +7457.5306. Composing quartets from 16 selected occupancies before alpha,
with every omitted pair retained in a positive remainder, gives +7432.5012:
only 25.0294 bits better. An optimistic boundary-gain diagnostic remains
at +6137.7494. That last number is a lower bound on this selected-quartet
bound, not a probability upper bound and not evidence of a low-weight word.
Thus this limited positive-point refinement does not resolve the profile
at these marks. It does not rule out a full-occupancy calculation or a
quartet coefficient bound that retains the alias improvement.

The [version-53 receipt](theory/hd_field128_progress_v53_20261010.json)
indexes the results, source reviews, tests and remote preservation.
Versions 54--56 now perform the [full quartet calculation](HD_FIELD128_FULL_QUARTET.md).
The point check streams all 17,850,625 ordered occupancy quartets. Its
positive-marker score is +7390.7396, versus +7457.5306 for pairs.
An independent positive coefficient engine groups occupancy multisets,
retains degree 128, and raises each regional quartet matrix to power 32.
Its four point controls agree with the streamed calculation. Replaying
all 4,198,401 profiles improves the complete score by 100.0589 bits,
to +3712.6341. All previous witnesses remain available. This calculation
does not close the midpoint weight and uses FLOAT, not directed arithmetic.

The next [return refinement](HD_FIELD128_RETURN_TAIL_MOMENTS.md) uses
factorial moments of the number of zero output coordinates. Exact zeros
in the moment bounds restrict that number's support. Normalized moment
clipping and a positive tail expansion then supply two additional caps.
The exact scan tightens 39 of 65 occupancies at the saved witness; its
largest return-cap gain is 39.4512 bits, at occupancy 17. Some strictly
positive improvements are too small to resolve in binary64 logarithms.

Version 59 independently reproduces the prior regional hashes and complete
envelope before propagating the new return caps through pairs and quartets.
The complete midpoint score improves by only 0.0621 bits, to +3712.5720.
The local gain therefore does not survive strongly enough to resolve the
45% proof. All 21 small quartet and moment tests pass remotely.
These results motivate a structural test of states larger than 32, with
fresh local bounds; they do not establish a low-weight codeword or rule
out a sharper proof of the existing construction.

The [wider-state proof](HD_FIELD128_WIDE_STATE.md) supplies those fresh caps
for states 32, 36, and 40. Emission interpolates from nodes 128..127+s
to 0..63; feedback is its transpose times a global degree-16 scalar.
The dimension-only return argument needs neither self-orthogonality nor
reversal. Its separate principal-setup risk is at most
`s * 2^64 / (2^128-2^64)`, below 2^-58 for all three states.
These are different constructions from the translated-feedback candidate.

Version 61 recomputes each state's exact local caps and numerical
coefficients using the seven saved witness settings. Every state retains
all 4,198,401 midpoint profiles. The complete scores are +3885.2698,
+3628.5460, and +3558.1095, respectively. State 40's worst profile is
(1252,1772). All five new tiny tests pass, including a nonzero-CA example
and independent return/live enumeration. The full-state comparison takes
141.87 seconds on Peach; this is proof-screen time, not encoder latency.

No cross-construction minimum is taken. None of these conditional FLOAT
screens closes the midpoint weight, and no wider-state encoder was built.
The [version-61 progress receipt](theory/hd_field128_progress_v61_20261010.json)
also preserves the earlier quartet and moment results. Next retune the
state-40 witnesses before widening again; the current settings were
selected for state 32, so this comparison does not establish an optimum.
Preserve the [version-37 receipt](theory/hd_field128_progress_v37_20261010.json)
and every certified rung; none of these diagnostics changes their claims.

For implementation, an exact completion identity can remove duplicate
fixed-map work. For a polynomial p over GF(2^128), if q_j=p(j) for j<30
and deg p<30, set
S=sum_j q_j and T=sum_j j*q_j, with integers denoting the established
GF256 byte labels. Then q_30=T+31*S and q_31=T+30*S. The corresponding
adjoint has the same shared-product simplification. Applying these
identities in emission, feedback, and refresh saves 24 vector GFNI
instructions per packet. The isolated completion variant now passes
exhaustive completion/adjoint basis tests, independent dense full-encoder
checks at K8192/16384, and same-plan controls. Matched K20 full timings are
28.221--28.263 ms versus 28.564--28.646 ms for the deferred predecessor:
about 1.3% improvement. All output checksums agree.

A separate fused route loader passes 4096 focused route tests and the same
dense whole-encoder checks, but is slower. Its full medians are
28.425--28.575 ms versus 27.919--27.978 ms for the matched completion
control; inner-only medians are 13.689 and 13.301 ms. Keep the completion
variant and park this route experiment. Measurements are serial, use
101 samples and three warmups, and exclude setup and allocation.
Arbitrary sampled diagonals still prevent cancelling fixed transforms
across refresh layers.

All new numerical outputs and their source closure are in a
[verified durable remote archive](theory/hd_field128_v34_remote_archive_20261010.json).
Only source and compact receipts are retained locally. The archive does not
replace the separate version-31 distance-certificate archives.
The [version-37 archive receipt](theory/hd_field128_v37_remote_archive_20261010.json)
separately preserves the new diagnostics and exact kernel variants.

## Literal lifted construction

Let F=GF(2^128), with a fixed embedding of AES GF256. Keep K=1048576,
N=2097152, and 8192 outer superblocks of 128 input symbols. Each superblock
uses four fixed RS[64,32] constituents. Their coefficients remain in the
embedded GF256. At every superblock and output column, independently
sample a uniform invertible 4-by-4 Toeplitz matrix over F.

Keep the 256 regions, alternating even/odd adjacent-pair matchings on the
cyclic superblock array, independent uniform permutations of the 4096
pairs in each region, and independent fair swaps within those pairs.

The inner retains 64 input symbols and a 16-symbol state per step. Its
fixed map A consists of four RS[16,4] maps; C=A^T and CA=0. H16 is the
same fixed MDS interpolation matrix, interpreted through the embedding.
Starting from Q=0, process each input X by

    Y = X + A Q,
    Q_next = D2 H16 D1 (Q + C X).

Every diagonal coefficient in D1 and D2 is independently sampled from
F minus zero. All setup choices are mutually independent and are fixed
for every encoded message. State persists across regions. The payload is
now one F-symbol per record, not sixteen independent GF256 coordinates.
Weight counts nonzero records, not nonzero bits.

This changes the sampled code. Merely applying a sampled GF256 code to
wider values preserves its distance and is not the proposed lift.

## Proof and implementation order

1. Audit the Toeplitz counting bound and local inner operators over F.
   Fixed small-field ranks survive the embedding, but numerical bounds,
   source multiplicities, and projective normalization change with q.
2. Build field-parameterized operators without mutating frozen GF256
   helpers. Reproduce the old local calculations at q=256, then screen
   selected source classes at the new field size using bounded runs.
3. For each promising distance, cover every source class and replay the
   full union with directed arithmetic. Explicitly record unresolved
   classes and distinguish a failed bound from an actual low-weight word.
4. Implement the same distribution with an explicit field representation
   and independent arithmetic/encoder oracles. Include necessary basis
   conversions in timings. A restricted fast coefficient family requires
   a new proof and must not replace uniform full-field sampling silently.
5. Run matched complete-encoder timings serially. Only after this literal
   comparison consider removing work using the additional distance margin.

## Initial audit, not a certificate

The outer atom cap generalizes to 1/((q-1)q^3), and the sharpened outer
density is q^160/(q-1)^32. The scalar projective divisor becomes q-1.
The two-diagonal refresh comparison uses (q^16-1)/(q-1)^16. These formulas
must be combined with the actual paired routing and source counts; they
do not themselves establish a distance guarantee.

The literal implementation performs 4,718,592 full-field products in the
outer and 1,048,576 in the refresh diagonals at K20. Nine prepared outer
coefficients per Toeplitz matrix plus the refresh diagonals occupy 88 MiB,
versus 5.5 MiB with byte coefficients. The equivalent compact-seven outer
representation lowers this to 72 MiB without changing the distribution or
the nine-product circuit. These are raw coefficient-storage counts, not
runtime measurements; additional preparation may require more storage.

## Historical first gates, 2026-10-09

The [proof audit](HD_FIELD128_LIFT_PROOF.md) now derives the full-field
outer density, refresh comparison, local transition, and projective
normalization. The sparse evaluator was independently reviewed. A
positive-expression repair prevents tiny positive state births from
disappearing through subtraction of nearly equal floating logarithms.

The new [selected-weight directed verifier](theory/hd_field128_lift_certificate.py)
uses exact witnesses, positive packet expressions, Arb coefficient
arithmetic, and hashes of all recursively imported local helpers.
These are completed selected-weight gates, **not whole-code certificates**:

| Active superblocks | Target symbol distance | Directed failure margin | Precision |
| ---: | ---: | ---: | ---: |
| 1 | 49% | 311.248692113942 bits | 256 and 384 bits |
| 16 | 45% | 411.987106344610 bits | 256 bits |
| 64 | 35% | 521.076451376465 bits | 256 bits |
| 320 | 35% | 889.355605165406 bits | 256 bits |
| 1024 | 35% | 1635.450460201077 bits | 256 bits |

Each row covers all supports and both alternating-matching statistics at
that exact active-superblock weight. Intermediate weights are not implied.
The corresponding receipts are indexed by
`hd_field128_r1_d49_directed{256,384}_20261009.json`,
`hd_field128_r16_d45_directed256_20261009.json`, and
`hd_field128_r64_d35_directed256_20261009.json` in `theory/`.
The weight-320 receipt is
`hd_field128_r320_d35_directed256_20261009.json`; it uses an
outward-rounded interval DFT, not FLOAT FFT coefficients. The same method
gives the weight-1024 receipt
`hd_field128_r1024_d35_directed256_20261009.json`.

FLOAT screens also pass 35% at weights 1, 2, 4, 8, 16, 24, 32, and 48.
The first witness optimization does not close 45% at weight 32, or 40%
at weight 64. These failed bounds do not exhibit low-weight codewords;
the weight-64 optimization also hit its explicit time budget. The
route-independent terminal calculation at 35% starts at weight 5854
and gives a 797.488150197979-bit bound for the entire range through
8192. The 256- and 384-bit directed replays agree at the recorded ceiling
(`hd_field128_terminal35_directed{256,384}_20261009.json`).

At this stage, the full-cover gap was the middle range. Replacing the old
coefficient extraction by a positive-marker envelope loses too much:
at weight 320, the complete old 22-marker starting bank gives a positive
log bound around 2863 even at 35%. A smaller seven-marker bank is much
worse. These runs justify restoring coefficient extraction; they do not
justify changing the construction yet. Sparse and dense endpoint
successes must not be combined into a full-distance claim.
Restoring interval residue aliases closes the selected weight-320 case
as recorded above. FLOAT alias screens also pass weights 2048, 3584,
and 3936 at 35%; their surrounding bands still require fresh directed
coverage. At weight 1024, the conservative FLOAT cancellation floor
prevented the initial proposal from passing; the directed interval
calculation resolves that numerical issue and passes as recorded above.

A first shared-witness sparse sweep evaluated every weight from 1 through
64, but seven weights needed different witnesses: 3, 6, 7, and 11 through
14. A separate directed replay repairs these and strengthens weight 41.
The resulting union covers **every weight 1 through 64**, with
130.42160124837 bits of failure margin at 35% distance. Its receipt is
`hd_field128_sparse35_repaired_directed256_20261009.json`; the original
failed sweep remains preserved. This is a complete sparse-range bound,
not a whole-code bound.

A further 48-evaluation FLOAT optimization at weight 64 and 40% distance
still gives a positive log bound, about 109.096. The previous failure was
therefore not explained merely by its short optimization budget. This
remains a limitation of the present bound, not a low-weight codeword.

### Explicit coverage and the first proof refinement

The first source-checked receipt index covered every weight 1 through 1441
and the terminal range 5854 through 8192. Its directed union over these
events has 108.117061854228 bits of margin. Weights 1442 through 5853
were unresolved in that snapshot; its union is **not** a whole-code
certificate. The index is
`theory/hd_field128_partial35_index_stage19_20261009.json`. It checks
construction metadata, source hashes, row provenance, and scalar
enclosures, but does not independently replay regional coefficients.

In particular, new carry-free coefficient tables close the sparse gaps
65--70, 194--196, and 210--256. Their combined directed margin is
207.847438993285 bits. The middle-range batch continues independently.
The old weight-3936 marker bank cannot simply be extended past 4096:
those histograms admit no empty pairs, whereas the bank was tuned away
from that boundary. Its failed continuation is retained as a failed
bound, not evidence of low code distance.

For the difficult weight-64 class at 40% distance, a proof-only refinement
composes four inner steps before taking the fractional power on selected
arrangements. It partitions all other arrangements into positive terms;
it never subtracts upper bounds. The encoder and its setup distribution
are unchanged. After 48 mark-search evaluations, the FLOAT log bound is
44.387599197648, versus 109.095731846190 for the two-step analysis.
This improves the bound by about 65 bits, but still does not certify even
this single weight class. A separate single-RS[64,16] inner candidate
gives only 84.327137034561 after the same search budget; it changes the
construction and has no implementation or certificate. The proof-only
refinement is the more promising next step.

Nested eight- and sixteen-step refinements retain that same construction.
At weight 64 and 40%, an optimized eight-step FLOAT bound is 16.400048416583;
the sixteen-step bound at those marks is 11.614611833618. Both still fail.
An independently justified regional Jensen bound was also tested: taking
alpha after averaging a whole region loses too much route information here
and does not improve the eight-step result.

### Larger-state candidates

A separate structural screen uses a single RS[64,s] emission map A and
C=A^T, for s=24 or 32. Evaluation on the additive 64-point subspace gives
CA=0 because products of the interpolation polynomials have degree at most
2s-2<=62. The refresh remains D2 H_s D1 with independent full-field
nonzero diagonal entries and a fixed totally nonsingular H_s. The outer
and route are unchanged. These are different constructions with different
state and arithmetic costs, not proof refinements of the literal lift.

At the retained marks z=0.4961394626874333 and
alpha=0.006849071540164115, complete weight-64 FLOAT screens give:

| State symbols | Target symbol distance | Log2 upper-bound screen |
| ---: | ---: | ---: |
| 24 | 40% | -734.998211211324 |
| 32 | 45% | -99.046175324030 |

Each row includes all 1088 feasible matching profiles at this one active
group weight. Neither row covers other weights, uses directed arithmetic,
or measures an encoder. A bounded selected-weight ladder and independent
directed replays are the next checks. No faster-than-RS claim follows.

The first 24-evaluation selected-weight ladder now passes weights
1, 2, 4, 8, 16, 32, and 64 for both candidates. The weakest screened
margin is 236.459269961354 bits for state 24 at 40%, and
97.130881169570 bits for state 32 at 45%; both occur at weight 1.
After mark search, the state-32, weight-64 bound improves to
262.984974430352 bits. These minima describe only the seven sampled
weights, not the intervening integers or the full code. Per-weight
checkpoints are retained in `theory/hd_field128_fullrs_state_ladder_20261010/`.

Directed replays now confirm the state-32, 45% endpoints. At weights 1
and 64, the upward log2 ceilings are -97.130881170152 and
-262.984974429386, respectively. Replays at 256- and 384-bit arithmetic
precision agree on these ceilings. They cover all matching profiles at
their selected weights, not the full weight range. The receipts are in
`theory/hd_field128_fullrs32_d45_directed_20261010/`.

### Wider-inner candidate at 48%

A second construction uses one RS[128,48] emission map, its transpose for
feedback, and the same diagonal--MDS--diagonal refresh on 48 state symbols.
The outer and global route are unchanged. Each region now has 64 inner
steps, each containing 64 pair slots. The additive 128-point evaluation
set again gives CA=0, since 2*48-2<127.

All seven selected sparse weights pass the first FLOAT screen at 48%
symbol distance:

| Active outer groups | Screened margin in bits |
| ---: | ---: |
| 1 | 158.396 |
| 2 | 1125.498 |
| 4 | 1125.366 |
| 8 | 4032.530 |
| 16 | 8468.468 |
| 32 | 823.888 |
| 64 | 157.581 |

These results are witness proposals, not a certificate or an estimate of
the actual minimum distance. Unsampled sparse weights, middle weights,
and dense weights still require bounds for this construction. Results are
in `theory/hd_field128_wide128_state48_d48_20261010/` and
`theory/hd_field128_wide128_state48_d48_fill_20261010/`.
The refresh uses 96 full-field products per 128 outputs: the same product
count per output as the width-64, state-24 candidate. Fixed-map work and
memory traffic have not been measured for either new construction.

Directed 256- and 384-bit replays confirm the weight-1 and weight-64
ceilings as -158.395755031746 and -157.580534145978. These agree across
precisions; they still cover only those two weights. Receipts are in
`theory/hd_field128_wide128_state48_d48_directed_20261010/`.

Failed bounded mark searches must not be interpreted as construction
obstructions. In particular, reseeding the width-128/state-64 weight-64
search at 48% changes its screen from +129.920 to -577.451. Likewise,
the width-64/state-32 weight-1 screen improves from about -24.090 to
-356.742 after a longer, better-seeded search. The wider inner is not
established as necessary. A 32-evaluation width-64/state-32 dense search
at weight 4096 and 45% remains positive (+20586.682), with reliable
selected numerical coefficients; this is an unresolved bound, not a
low-weight codeword. Dense coverage is the next discriminating test.

At those dense marks, width-128 states 48 and 64 also have positive
weight-4096 screens at 48% (+30211.204 and +29601.820). All selected
coefficients pass the FLOAT reliability checks, but these fixed-mark
failures do not rule out better witnesses. In contrast, the width-64,
state-24 candidate passes weight 4096 at 40% with a -12849.597 screen.
Together with its sparse results, that makes 40% the next whole-cover
target while the higher-distance dense bound is investigated.

### The 40% middle range and a fixed feedback scalar

The state-24 candidate does not yet have a complete 40% certificate.
Forty-evaluation searches at weights 256, 1024, and 2048 leave positive
FLOAT log bounds of 109.449, 1964.090, and 3458.675. Their selected
coefficients do not use the cancellation-floor fallback. The fresh
directed terminal bound covers weights 6678 through 8192 at 40%, but
cannot resolve these middle weights. The searches are retained in
`theory/hd_field128_state24_d40_search_r{256,1024,2048}_20261010/`;
the terminal receipt is `theory/hd_field128_fullrs24_d40_terminal256_20261010.json`.

A separate [fixed-scalar feedback construction](HD_FIELD128_GAMMA_RETURN.md)
replaces C by gamma C. Here gamma is the degree-16 tower generator Y.
For the full-RS emission maps with states 16, 24, 28, or 30, a determinant
identity proves that every active principal feedback matrix is invertible.
This permits a tighter joint count of active and inactive output zeros
during a return to zero state. It introduces no bad-setup event.

At the retained marks, gamma alone changes the state-24, weight-256
40% screen only from +109.449 to +107.475. It is not a sufficient fix.
State 30 passes the selected weights 256 and 2048 at 40%, with log bounds
-406.489 and -3278.604, but weight 1024 remains positive at +1246.639.
These are new-construction FLOAT screens, not certificates or timings.

Two proof refinements diagnose the remaining loss without changing this
new encoder. First, retain the chronological regional matrix products
instead of bounding each region by a common-potential contraction.
Second, restore the sparse proof's sharper support count. In particular,
a matching with no singleton active groups can be contracted before
counting cyclic runs. At the original state-30, weight-1024 marks, these
changes reduce the log bound from +1246.639 to +1000.483 and then
+454.950. The bound is still positive. The corresponding receipts are
`theory/hd_field128_gamma30_d40_r1024_{initial,transfer,profile}_20261010.json`.

The [all-pair composition proof](HD_FIELD128_ALL_PAIR_COMPOSITION.md)
justifies composing every chronological packet pair before applying the
fractional power. The earlier good-pair restriction was computational,
not a probabilistic requirement. Combining this refinement with the
sharper profile count gives +429.158 at state 30, weight 1024, and 40%
distance. At state 24 and weight 256, it improves +70.205 to +9.582.
Neither is a passing bound. These fixed-mark screens are retained as
`theory/hd_field128_gamma30_d40_r1024_allpairs_initial_20261010.json` and
`theory/hd_field128_gamma24_d40_r256_allpairs_initial_20261010.json`.
Marks optimized for the common-potential bound can be worse here. A
fresh 40-evaluation search of the refined objective reduces the
state-30, weight-1024 log bound to -63.741904841387. This is a passing
FLOAT proposal for that weight, not a whole-code certificate. The
corresponding state-24, weight-256 search reaches -28.572269756369,
short of the requested 40-bit margin. Both searches retain all 79
source markers and use no numerical fallback in their selected
coefficients. They are saved in
`theory/hd_field128_gamma30_d40_r1024_allpair_search_20261010/` and
`theory/hd_field128_gamma24_d40_r256_allpair_search_20261010/`.
Fresh directed replays of the state-30 proposal at 256 and 384 bits agree
on the log ceiling -86.724953795603. Each replay covers all 263,169
even/odd matching-profile pairs at weight 1024. An independent 384-bit
aggregation of their 513 immutable row checkpoints reproduces the same
ceiling, and both source closures match their pinned bytes. The
[selected-class checkpoint](theory/hd_field128_gamma30_d40_r1024_index_20261010.json)
indexes the witnesses and receipts. This certifies that class at 40%
distance; it does not cover the other weights or provide a timing for
the gamma construction.

The same construction's terminal event is also closed. Fresh 256- and
384-bit inverse-image replays cover weights 6678 through 8192 together,
with log ceiling -2169.486909120572. Identity feedthrough gives a unique
inverse for every fixed gamma/route/inner setup, so this argument uses
the unchanged outer counting law rather than a state-specific local
cap. The terminal term is added once, not once per tail weight. Receipts
are `theory/hd_field128_gamma30_d40_terminal{256,384}_20261010.json`.

Every sparse weight 1 through 64 now has a directed same-construction
bound. Broad bands with a single center's mark initially failed away
from that center. Choosing weight-specific exact marks and recomputing
the coefficients closes those gaps; no bound is interpolated. The
directed union over the selected 64 receipts has 85.087828627317 bits
of margin. The [first partial-cover index](theory/hd_field128_gamma30_d40_partial_cover_20261010.json)
pins those receipts, the weight-1024 receipt, and the single terminal
event. That immutable index remains preserved.

Fresh all-profile computations also close weights 1023 and 1025, using
the same directed coefficient bank and exact marks as weight 1024.
Their respective failure margins are 79.586099545838 and 109.865491904278
bits. An independent audit checked complete profile coverage, source
pins, and the stored row enclosures; it did not replay the coefficients.
The [second partial-cover index](theory/hd_field128_gamma30_d40_partial_cover_v2_20261010.json)
contains 68 events and counts the terminal event once. Its covered
ranges are 1--64, 1023--1025, and 6678--8192. Its 79.544623499091-bit
union over covered events is explicitly not a whole-code margin.

Fresh banks calibrated at weights 256, 2048, and 4096 close additional
contiguous intervals. Each integer weight receives fresh normalized
regional bounds and a common-potential bound for all matching profiles;
the intervals are not inferred from their endpoints.

| New interval | Number of weights | Combined failure margin |
| --- | ---: | ---: |
| 241--272 | 32 | 81.284790953863 bits |
| 1921--2240 | 320 | 871.997857507755 bits |
| 4065--4128 | 64 | 28425.618964266783 bits |

The upper interval imports only proposed marks from an earlier 45%
screen. Its coefficient bank and bounds are recomputed for 40%.
Independent audits check source and artifact hashes, exact proposals,
integer coverage, and directed union arithmetic at 384 bits. They do
not repeat coefficient generation.

The [version-31 index](theory/hd_field128_gamma30_d40_partial_cover_v31_20261010.json)
contains 274 events over nonoverlapping weight ranges and preserves all
265 version-30 entries. It covers every weight from 1 through 8192;
the terminal event is counted once. Its actual union gives
49.171469646557 bits of margin. The separate
[completion manifest](theory/hd_field128_gamma30_d40_completion_v31_20261010.json)
records the whole-code conclusion; the historical partial indexes remain unchanged.
The last 61 weights were covered by composite bands 521--547 and 549--575,
and direct sums at 917--919 and 1146--1149. All nine results passed
independent saved-data audits and have verified durable remote archives.
The new bands give 54.411239507314 and 51.148161164307 bits, respectively.

For historical comparison, version 30 covered 8131 weights and its covered
union gave 49.646477666031 bits; that was not a whole-code certificate.
Version 30 adds weights 548, 914, 915, 916, and 1145. Their respective
directed margins are 464.0724, 82.2270, 98.8402, 126.3709, and 49.6942 bits.
Every new result has a separate count/context/saved-sum audit and a
verified durable remote archive. Weight 1145 uses an explicit 48-bit
local target; the fresh union of actual event ceilings meets the global
40-bit target. Previous event ceilings are unchanged.

Version 30's SHA256 is
`875c572f4ae84b5bc139a4c21be9b0b4e18aed51ab1f288c45c125a260aab88f`.
The [explicit selection](theory/hd_field128_gamma30_d40_v30_selection_20261010.json)
pins the five compact audit/archive pairs.
The [independent version-30 check](theory/hd_field128_gamma30_d40_partial_cover_v30_20261010_independent_check.json)
reconstructs the saved index, preserves all 260 parent events, and
separately recomputes the positive union at 384 bits. It does not replay
the remote coefficient banks or transfers.

Versions 28 and 29 add five independently audited remote bands:

| Added range | Weights | Directed union margin |
| --- | ---: | ---: |
| 165--179 | 15 | 151.623005228399 bits |
| 704--757 | 54 | 65.381520194504 bits |
| 968--999 | 32 | 60.771257274185 bits |
| 1150--1159 | 10 | 58.170620153978 bits |
| 920--967 (version 29) | 48 | 59.996073946746 bits |

Version 29's SHA256 is
`ea1900e2f869a258f16404bd34c3882675a78e960e2703a3754820cc4e63737d`.
The [explicit version-29 selection](theory/hd_field128_gamma30_d40_v29_selection_20261010.json)
pins the five audit receipts and their verified durable remote archives.
All 259 version-28 events remain unchanged. The
[independent version-29 check](theory/hd_field128_gamma30_d40_partial_cover_v29_20261010_independent_check.json)
reconstructs the saved object and separately recomputes its union at 384 bits.
The mathematical tree audits and archive checks were performed remotely;
this local check verifies their pinned compact receipts, not the remote trees.
The remote audits recompute integer counts and check complete partitions,
checkpoint transitions, fresh-bank metadata, and saved directed ceilings.
They do not regenerate coefficients or transfers. Exact seed literals are
the input; historical calibration bounds are neither verified nor reused.

Version 27 adds the independently audited bands 209--215 and 1160--1191,
with respective bounds of 74.478423261615 and 70.941498885492 bits.
Version 27's SHA256 is
`bb5873cc96c861d3c9f08806bdbae63c5778fd9d3d57999a2631f529dd5e3a5f`.
The [explicit version-27 selection](theory/hd_field128_gamma30_d40_v27_selection_20261010.json)
binds those two bands; all 253 version-26 records remain unchanged.
The [independent version-27 check](theory/hd_field128_gamma30_d40_partial_cover_v27_20261010_independent_check.json)
reconstructs the saved index and recomputes its positive union at 384 bits.
It does not repeat coefficient generation.
Version 26 adds three audited bands and the separately audited weight 913.
Its SHA256 is
`45d2854157e47e23145f76df445bc8cae4f318ce828003e8322eb36307ea1f20`.
The [explicit version-26 selection](theory/hd_field128_gamma30_d40_v26_selection_20261010.json)
binds those four events; all 249 version-25 records remain unchanged.
Version 25 adds the audited 368--383 event and leaves the first 248 entries
unchanged. Its SHA256 is
`c92b67290ff9b1a729496598e071a7599a3c76ce2486e6441023ab3dee59aea0`.
The [explicit selection](theory/hd_field128_gamma30_d40_v25_selection_20261010.json)
records all seven additions since version 22.
The [independent version-25 check](theory/hd_field128_gamma30_d40_partial_cover_v25_20261010_independent_check.json)
rebuilds the saved object exactly and rechecks 2036 input pins and 149
helper pins. It checks saved partitions and directed unions, not fresh
coefficient or transfer calculations.
Versions 23 and 24 add six audited bands covering 129 weights to version 22.
The [explicit version-24 selection](theory/hd_field128_gamma30_d40_v24_selection_20261010.json)
distinguishes its four new bands from the two version-23 pilot bands.
The first 244 version-23 entries are unchanged. The version-24 index has SHA256
`703252867f063c8bff2b76f83c582dbb8675577e6a2f56d7e07093b6430d5c45`.
The [independent version-24 check](theory/hd_field128_gamma30_d40_partial_cover_v24_20261010_independent_check.json)
rebuilds the saved object exactly, rechecks the new partitions, and
recomputes the 384-bit union. It verifies 1,225 input pins and 148 helper
source pins without repeating coefficient, regional, or transfer evaluation.

Previously, version 22 retained all 227 version-21 events and added 15 independently
audited bands covering 152 weights. Their combined bound alone has
66.773329912796 bits of margin. The
[explicit selection](theory/hd_field128_gamma30_d40_v22_selection_20261010.json)
admits only complete bands; partial and storage-interrupted runs are excluded.
The [independent version-22 check](theory/hd_field128_gamma30_d40_partial_cover_v22_20261010_independent_check.json)
passes: it checks exact retention of the parent events, the integer weight
partition, source bindings, the 15 new saved subdivision trees, and fresh
384-bit unions. It does not repeat coefficient generation, regional bounds,
or cell-transfer calculations; their reviewed producers and prior audits
remain trusted.

The preceding 37 bands each have an independent count, partition,
provenance, and saved-ceiling audit. Version 21 preserves all 190 version-20 events
and freshly aggregates their union with the new bands. Its 853 new
weights close every gap at or above 1192. The new events alone have a
combined 69.146046963364-bit bound. The
[explicit selection](theory/hd_field128_gamma30_d40_v21_selection_20261010.json)
records all 37 additions and excludes the incomplete attempts. It
distinguishes the new complete 3904--3936 event from its earlier failure.
The [independent v21 check](theory/hd_field128_gamma30_d40_partial_cover_v21_20261010_independent_check.json)
passes. It preserves all 190 parent events, checks the 37 new integer
partitions, freshly aggregates saved leaf and event bounds at 384 bits,
and independently reconstructs the index. It verifies 9,556 file pins
and 134 helper sources, including the separate retry-driver receipts.
It binds the existing independent count audits without repeating them;
it does not replay coefficients, regional bounds, or cell transfers.
The [independent v20 check](theory/hd_field128_gamma30_d40_partial_cover_v20_20261010_independent_check.json)
remains preserved. It checks parent preservation, exact coverage, its 36 new integer
partitions, input and source bindings, and fresh 384-bit unions. It
inherits and binds the independent count audits rather than repeating
their counts. It does not regenerate coefficients or regional transfers.
The [independent v19 check](theory/hd_field128_gamma30_d40_partial_cover_v19_20261010_independent_check.json)
remains preserved. It verifies all 125 inherited events, its 29 new integer
partitions, input and source bindings, and 384-bit leaf and aggregate
unions. It inherits the bands' independent count audits rather than
recomputing those counts again. The check also excludes both the partial
5344--5407 trial and the resource-stopped 4384--4447 trial.
The [independent v18 check](theory/hd_field128_gamma30_d40_partial_cover_v18_20261010_independent_check.json)
remains preserved: it verifies parent preservation, exact coverage, its six new
integer partitions, source bindings, and 384-bit leaf and aggregate unions.
The [independent v17 check](theory/hd_field128_gamma30_d40_partial_cover_v17_20261010_independent_check.json)
remains preserved, including its checks of the 113 version-16 entries.
The earlier
[independent v16 check](theory/hd_field128_gamma30_d40_partial_cover_v16_20261010_independent_check.json)
remains preserved. These checks do not replay coefficient or regional calculations.

Several new bands use an [immutable subset of coefficient markers](theory/hd_field128_gamma_marker_subset.py).
Each retained marker supplies its own valid upper bound; dropping markers
can only weaken the minimum. Retaining 47 of 79 markers closes 2241--2256
and 4129--4192 with respective margins of 1502.571444791055 and
26801.485981291909 bits. The checker preserves all tail markers and
every second earlier marker. Uniform thinning to 21 markers failed the
tested weights. A 31-marker subset passed weights 2048 and 4096 but
failed weight 1922; these overlapping probes are not counted again.
The subset module passes six additional tests, bringing the geometry,
coverage, batch, and subset checks to 58 tests. Independent receipt audits
check lineage, complete weight coverage, source pins, and directed unions;
they do not recompute the coefficient atlases.

Further bands extend the upper covered interval to 3937--4309.
The complete 4257--4320 check fails, even with all 79 markers. Its first
53 rows nevertheless give a valid 323.951548517653-bit bound for
4257--4309. A separate immutable receipt identifies those rows and their
source; the failed full-band event is not included in the index.
Seven complete adaptive covers also add weights 1026--1032, each with
at least 61.957 bits of margin. Those additions lowered the earlier
partial-union margin to about 60.077 bits. An adaptive probe at weight 1016 reached
its 4,096-cell budget without closing and contributes no bound.

### Sharper counts and better coefficient markers

Two separate sources of slack affect the uncovered weights. Near a
matching with few singleton pairs, the previous support count ignores
constraints on the other matching. At weight 768, some coefficient
markers also give loose bounds near the histogram boundary. Neither
issue requires changing the encoder.

The [singleton-deletion count](theory/hd_field128_gamma_near_face_counts.py)
handles the first issue. Fix a binary support on P indexed pairs, with
weight r, f full pairs, and k cyclic runs of ones. There are s=r-2f
singleton pairs. Delete them and contract each remaining equal pair.
The resulting indexed circular word has length P-s, weight f, and
between max(0,k-s) and k runs. The deleted positions, their orientations,
and this residual word reconstruct the original support. Thus an upper
count is C(P,s) times 2^s times the sum of those residual run counts.
We take the minimum of both matchings' bounds and the previous bound.
Constant and empty residual words have zero runs; no rotation factor is
introduced. A singleton cutoff skips optional improvements, not profiles.

With a cutoff of 32, the [directed replay at weight 512](theory/hd_field128_gamma30_d40_r512_near_face_directed256_20261010/summary.json)
covers all 66,049 profiles and gives 75.023496054901 bits of margin.
Its fresh coefficient bank uses 256-bit directed arithmetic. An
independent audit checks every exact support count, all 257 row
receipts, source pins, and a 384-bit aggregation. It does not rebuild
the coefficient bank. The count implementation passes five tests,
including exhaustive small supports; its directed wrapper passes five
additional tests. The same witness fails at weight 504 but passes at
520. A subsequent complete check of 513--519 gives a 94.704316415419-bit
union bound. Together with the separately audited endpoints, it certifies
512--520 without interpolation.

At weight 768, a direct FLOAT extraction of the all-singleton coefficient
reveals up to 10.56 bits of slack per regional matrix entry. Four
optimized positive markers recover almost all of that selected-profile
gain. Adding four markers for the opposite boundary exposes a thin
near-boundary strip as the next obstruction. A further 32 markers reduce
the complete selected-weight FLOAT log bound from +264.51 to +17.80,
using the same z, alpha, construction, and 32-singleton count cutoff.
All 148,225 profiles remain in the sum. These diagnostics are recorded
in the [119-marker witness](theory/hd_field128_gamma30_d40_r768_near_axis119_float_20261010.json);
they are not directed certificates.

Increasing the singleton cutoff to 128 does not close that witness.
Eight more coefficient markers help only slightly. Retuning z and alpha
against the improved bound is decisive: the fixed 127-marker bank reaches
-219.3476 in FLOAT after 40 evaluations. A fresh directed bank uses the
earlier passing witness from evaluation 19, whose FLOAT bound was -215.1790.
The [directed weight-768 receipt](theory/hd_field128_gamma30_d40_r768_near_face_directed256_20261010/summary.json)
then gives 250.943088184407 bits of margin, using a cutoff of 64.
It covers all 148,225 profiles. An independent audit verifies all 385
row receipts, exact marker conversion, source pins, and a fresh 384-bit
aggregation. The audit does not rebuild the regional coefficients.
These changes affect proof evaluation, not the encoder or its runtime.

The neighboring weights require explicit checks. The complete 752--767
band does not meet its target: weights 752--758 fail the 54-bit admission
threshold. A separate saved-row slice gives 54.94509398298 bits for
759--767; the complete 769--784 band gives 274.537452292099 bits.
Independent audits check the linked contexts, exact counts against a
separate formula, all row scopes, and 384-bit probability unions.
The index includes these two disjoint events and the separate weight-768
event, while excluding the failed full-band event. The newly added
54.945-bit suffix accounts for the lower partial-union margin.
Four subsequent complete bands cover every integer weight from 785 to
912. Their respective union margins are 552.902913848983, 862.375027019773,
687.586712126353, and 245.474911195594 bits. Each band contains 32 weights.
The saved-row audits check the original inline proposals against their
immutable proposal files, exact counts, coverage, and fresh 384-bit unions.
The v13 index includes all four bands, not the redundant probes at 832
and 896. No coefficient regeneration occurs in those audits.

The [exact joint-count formula](HD_FIELD128_EXACT_PROFILE_COUNTS.md)
offers a further refinement. It counts indexed supports by the transition
multiplicities between four vertices: even/odd position and zero/one bit.
A positive last-exit-tree formula, summed over feasible flows, gives the
count for each pair of matching statistics. Five unit tests and an
independent review pass. A separate bitmask oracle checks all 65,536
length-16 supports, all 145 feasible-marginal profiles, and their binomial
weight totals. For two previously dominant weight-768 profiles, the old
count upper bounds exceed these exact counts by about 140 and 176 bits.
These are individual-term improvements, not new probability certificates.
The [selected-count evaluator](theory/hd_field128_gamma_exact_refinement.py)
and its directed wrapper now implement that replacement. They preserve
all profiles and pass 11 tests plus independent review. An optional
FLOAT selector ranks canonical pairs by both ordered contributions;
the count identity remains valid regardless of how the selection was made.

The [midpoint-count lemma](HD_FIELD128_MIDPOINT_COUNTS.md) removes the
need to select profiles for most of this gain. After complementing and
exchanging matchings, put r<=P, a>=b, and s=r-2a. The positive fixed-flow
counts are symmetric and nondecreasing up to the midpoint, apart from
the explicitly handled two-alternating-word case. Their sum is at most
(s+1) times the midpoint term. At P=4096 the factor is at most 4097.
The midpoint's four vertex degrees do not depend on a or b. Two
one-dimensional binomial arrays therefore give constant-work integer
profile queries, without a quadratic count table.

Two independent reviews and nine tests cover the lemma, scalar and
vector count interfaces, complete directed transfer, and partial-run
behavior. A further seeded check compares 224 large profiles with the
exact sum, including weights 128 through 6144. All pass; their largest
observed log ratio is 5.130 bits. This observation is not the proof of
the uniform 4097 factor.

Fresh [directed midpoint receipts](theory/hd_field128_gamma30_d40_midpoint_directed256_20261010/weight_0768.json)
use the existing actual coefficient banks. The complete margins are
61.323299692641 bits at weight 511, 25.111068379316 at 756,
50.495929115151 at 757, 74.440374599393 at 758, and
299.84353261793 at 768. Weights 756 and 757 miss the index's 54-bit
event threshold. Version 14 imports weights 511 and 758; it retains
the earlier overlapping weight-768 event rather than counting it twice.
The complete weight-768 FLOAT margin is 299.8426 bits. Thus the new
count removes slack without changing the encoder or coefficient bank.
The [saved-receipt audit](theory/hd_field128_gamma30_d40_midpoint_selected5_audit_20261010.json)
checks every integer count against an independent central-flow formula
and reaggregates every saved row with a fresh outer prefix at 384 bits.
It confirms each target flag; it does not regenerate coefficients.
Taking the union of all five receipts would give only 25.111 bits,
so an audit pass does not mean they all meet the desired margin.
At fixed coefficients, replacing the midpoint caps by exact joint counts
can gain at most log2(r+1) bits. In particular, exact counts alone cannot
raise weight 756 from 25.111 bits to 40 bits; its witness must improve.

Fresh marker calibration is separately important. At weight 384, reusing
the 127 markers from weight 768 gives a vacuous FLOAT bound. Fitting
76 new positive markers at unchanged proposed z and alpha gives 276.678
bits. The [fresh directed replay](theory/hd_field128_gamma30_d40_r384_near_face_directed256_20261010/summary.json)
gives 276.747130952814 bits over all 37,249 profiles. A separate 384-bit
row aggregation reproduces that ceiling; it does not recompute coefficients.
Version 14 includes this selected weight.

The [fresh directed weight-128 replay](theory/hd_field128_gamma30_d40_r128_near_face_directed256_20261010/summary.json)
gives 757.535612086566 bits with 60 markers. Its separate 384-bit
saved-row audit reproduces the ceiling over all 4,225 profiles.
A fresh numerical proposal passes weight 1536 at 595.236 bits with
92 markers. Four rounds of targeted marker enrichment remove the
initially vacuous bound while retaining the same z and alpha. Its
[directed replay](theory/hd_field128_gamma30_d40_r1536_near_face_directed256_20261010/summary.json)
gives 595.311557743521 bits over all 591,361 profiles. The separate
384-bit saved-row audit agrees. Weights 128 and 1536 are now in v14.
At weight 3072, a fresh 76-marker proposal
and a 40-evaluation z/alpha search still fail. The best FLOAT log bound
is +23605.345, with 246 selected entries using the loose positive
fallback. A [four-marker diagnostic](theory/hd_field128_gamma_r3072_f1463_four_marker_diagnostic_20261010.json)
then identifies a gap between targets 1344 and 1520. Four new fits at
the dominant histogram 1463 recover about 150 bits per regional entry.
A fresh directed four-marker calculation agrees with the reliable FLOAT
entries within 3.3e-8 bits. This local diagnostic alone is not a complete
weight bound. The subsequent complete FLOAT replay with 80 markers
gives 2274.503 bits, without changing z or alpha. A fresh
[directed replay](theory/hd_field128_gamma30_d40_r3072_midpoint_directed256_20261010/summary.json)
then gives 2274.502892464027 bits over all 2,362,369 profiles.
The bank and gate take about 250 seconds, including 108 seconds for
the complete transfer. Its [separate audit](theory/hd_field128_gamma30_d40_r3072_midpoint_directed256_20261010_audit.json)
checks all integer profile counts, external row/context links, exact
marker lineage, and a fresh outer prefix with saved-row aggregation
at 384 bits. It reproduces the ceiling without replaying coefficients
or transfer matrices. Version 15 adds only this selected weight to v14.

The [midpoint rectangle method](HD_FIELD128_MIDPOINT_BOXES.md) avoids
enumerating every profile when a coarser enclosure suffices. Six exact
integer prefix arrays sum the raw midpoint count over a rectangle.
The frozen evaluator takes componentwise maxima only after regional
normalization. A new wrapper refines its support count and preserves
the chronological even/odd matrix product. Independent review and seven
tiny and mocked tests check counts, complementation, zero-count events,
batch pin checks, and source-mutation rejection.
Its [first directed cover](theory/hd_field128_gamma30_d40_r768_midpoint_cell_cover_20261010.json)
closes weight 768 at 68.823082024718 bits using 299 evaluated cells.
The cover takes 1.836 seconds after table preparation, versus roughly
7 seconds for exhaustive profile transfer. These are illustrative proof
processing times, not matched encoder benchmarks. Its
[independent audit](theory/hd_field128_gamma30_d40_r768_midpoint_cell_cover_independent_audit_20261010.json)
checks all 299 cell counts and the exact disjoint partition, including
the impossible all-full profile excluded by the existing roots. Its
151 accepted leaves cover every possible support; a fresh 384-bit
union reproduces the ceiling. The audit reuses saved regional bounds
and does not regenerate coefficients. The stronger existing
weight-768 receipt remains in the index. A bounded weight-928 pilot
does not close within 30 seconds; grouped profiles alone do not replace
calibration at a different weight. A second 30-second pilot at weight
1000 also remains partial. The subsequent exhaustive directed midpoint
sum completes in 10.7 seconds with a vacuous log bound of +337.083707264.
Thus this weight needs a better witness or bound, not just more time
subdividing the same cell envelope. Exact-count refinement alone cannot
recover a gap of this size.

### Contiguous bands and the weight-1000 repair

The [weight-band extension](HD_FIELD128_MIDPOINT_WEIGHT_BANDS.md) now
groups both weights and matching profiles. It uses one actual directed
coefficient bank, but normalizes every integer weight separately.
The outer factor is distributed across the 256 regional matrices before
taking componentwise maxima. Exact midpoint-count sums and a three-axis
integer partition then bound each cell. This preserves the chronological
matrix product and includes singleton-free faces; no weight interpolation
or state reset is used.

The new checker and independent auditor each pass seven tests. The
auditor uses a separate five-prefix count formula and an independent
tree traversal. It recomputes every saved cell count, checks proposal and
checkpoint bytes, and reaggregates the accepted ceilings at 384 bits.
It does not regenerate coefficients or regional transfer matrices.

Eight complete bands add 294 weights to version 16:

| Newly covered weights | Combined failure margin |
| --- | ---: |
| 1000 | 75.126 bits |
| 1001--1022 | 74.771 bits |
| 1033--1048 | 78.631 bits |
| 3008--3071 | 110.280 bits |
| 3073--3103 | 156.487 bits |
| 3104--3135 | 451.393 bits |
| 3136--3199 | 650.612 bits |
| 3200--3263 | 153.245 bits |

Every row has its own immutable summary, context, integer-cell checkpoints,
and independent audit, linked by the v16 index. The inherited receipts
for 1023--1032 and 3072 fill the intervening gaps without being counted
again. Thus the index now covers every weight in 1000--1048 and
3008--3263. The new aggregator passes five additional tests and
independent review. All earlier indices and failed computations remain
preserved.

The weight-1000 improvement comes from proof calibration, not a code
change. Marker enrichment alone improves the vacuous bound but does not
close it. Alternating exact-positive marker proposals with bounded
z/alpha searches gives a complete FLOAT log bound of -141.9465.
All 251,001 profiles are included. The final proposal has 127 markers;
its lineage and all 40 final search trials pass an independent audit.
A fresh directed bank and the grouped cover then certify 75.126 bits.
The weaker grouped margin is sufficient for coverage; it is not a claim
that the exhaustive directed sum equals the FLOAT result.

This round is substantive proof progress: 294 additional integer weights
are certified, including a previously vacuous class. It leaves 5,399
weights unresolved and does not establish 40% whole-code distance.
No setup distribution, encoder, or timing changed. The five bands around
3072 take roughly 38--86 seconds each to prepare their normalized tables,
then less than one second each for the grouped cover. These are proof
processing times after coefficient-bank generation, not encoder benchmarks.

### Rightward extensions and the dense boundary

Version 17 adds 257 weights using the same construction and retained
directed coefficient banks. Each weight receives its own histogram
normalization. The six complete events are:

| Newly covered weights | Combined failure margin |
| --- | ---: |
| 1049--1112 | 80.220759 bits |
| 1113--1144 | 60.506005 bits |
| 3264--3327 | 452.479659 bits |
| 3328--3391 | 136.787931 bits |
| 3392--3423 | 90.688487 bits |
| 4310 | 74.626205 bits |

All six independent audits pass. They recompute raw midpoint counts,
reconstruct the integer partitions, verify linked files and source pins,
and reaggregate accepted ceilings at 384 bits. They do not regenerate
regional coefficients or reevaluate the saved transfer bounds. The thin
version-17 aggregator passes five tests and independent source review.

The first 1113--1144 run required 80 bits at every cell. It stopped with
four unresolved singleton cells at weight 1144, whose margins range
from 63.297 to 74.702 bits. A separate run with a 60-bit cell threshold
completes the same domain and gives the combined bound in the table.
The 54-bit band-admission threshold is unchanged. The initial partial
receipt remains saved and is excluded from the index.

The [right-neighbor scan near 3072](theory/hd_field128_gamma30_d40_r3072_right_neighbors_float_20261010.json)
and [scan near 1000](theory/hd_field128_gamma30_d40_r1000_far_right_neighbors_float_20261010.json)
guided these extensions. Their sampled outcomes do not certify intervening
weights. The directed bands, rather than those samples, provide coverage.

The [dense-neighbor diagnostic](theory/hd_field128_gamma30_d40_r4096_dense_neighbors_float_20261010.json)
evaluates all profiles at original weights 4310, 4352, 4480, 4608, 4864,
5120, 5632, and 6144. The same old witness passes only 4310 in this scan;
its FLOAT log bounds grow from -73.757 there to +444046.201 at 6144.
Every selected coefficient entry is reliable under the screening policy,
with no selected fallback. At each sampled weight r, the largest profile
has a=b=r-4096: neither matching has empty pairs. This identifies a
target for fresh coefficient markers, not a low-weight codeword.
Only support counts are complemented above 4096; the matrices,
normalization, and outer factor retain the original weight.

That diagnostic motivated explicit boundary calibration with finite
positive marker seeds. At weight 5120, the boundary histogram has
3072 singleton pairs and 1024 full pairs. The follow-up below fits
this histogram and interior histograms before optimizing the complete
midpoint sum. It does not import literal-lift numerical bounds.

A separate [weight-1280 attempt](theory/hd_field128_gamma30_d40_r1280_midpoint_mark_search_20261010/recovery_summary.json)
starts from the weight-1536 witness. Eight marker-enrichment rounds and
28 completed z/alpha trials do not close it: the best FLOAT log bound
is +7919.325. The next trial hit a memory-allocation error. All completed
trials and exact proposals remain preserved; the incomplete trial is not
counted. Releasing completed normalized-table caches recovered about
5 GiB of committed memory while retaining the actual coefficient banks.
Fresh weight-1280 markers or a different seed are needed before directed replay.

The version-17 goal turn made proof progress: six audited additions cover 257 new
integer weights. The 40% whole-code certificate remains open. No encoder
benchmark, construction change, production promotion, commit, or push
accompanies this update.

### Fresh dense markers and exact boundary extraction

The dense follow-up changes the proof calculation, not the code.
[Exact empty-pair coefficient extraction](HD_FIELD128_EMPTY_PAIR_COEFFICIENTS.md)
reduces the no-empty regional boundary to one variable. At the original
weight-4096 marks, direct extraction proves a 264.166-bit bound for the
weight-5120 matching profile with 1024 full pairs in each region.
Fresh 256- and 384-bit coefficient replays agree. The independent audit
recomputes the support count, prefix, and final transfer from saved
regional caps. This isolated result is not added as a whole-weight event.

For complete weight coverage, the
[fresh dense calibration](theory/hd_field128_gamma30_d40_r5120_dense_calibrated_float_20261010.json)
adds 32 positive markers, reducing the full FLOAT log bound from
+195,673.94 to +67,906.23. A bounded
[mark search](theory/hd_field128_gamma30_d40_r5120_dense_mark_search_20261010/summary.json)
completes 31 trials and reaches +26,384.16. Its next trial times out
without a completed objective. All trials keep the original weight
in their coefficients and outer factor.

The [adaptive follow-up](theory/hd_field128_gamma30_d40_r5120_dense_adaptive_float_20261010.json)
then fits the current worst profile after each marker addition. Six
additions produce a complete FLOAT log bound of -12,000.77997.
Every completed step retains all 2,362,369 profiles. Selecting 18 markers
that preserve every regional entrywise minimum and rebuilding from their
exact literals reproduces the complete score and regional-array hash.
The [independent consistency audit](theory/hd_field128_gamma30_d40_r5120_winner18_replay_20261010/independent_consistency_audit.json)
checks that provenance and equality; it does not rerun the numerical screen.

A fresh directed 18-marker bank supplies the subsequent band proofs.
Every original weight receives its own coefficient normalization and
outer factor. These calculations use the existing reviewed integer-band
checker, not interpolation from sampled endpoints.

| Newly covered weights | Combined failure margin |
| --- | ---: |
| 4832--4895 | 12948.867947 bits |
| 4896--4959 | 12672.594272 bits |
| 4960--5023 | 11837.706145 bits |
| 5024--5087 | 10536.030837 bits |
| 5088--5151 | 8760.479455 bits |
| 5152--5215 | 5973.512585 bits |

Each band closes at its two parity roots, covering every feasible matching
profile for all 64 original weights. Independent audits recompute counts,
check partitions and source bindings, and aggregate saved ceilings at
384 bits. They do not regenerate regional coefficients. Version 18 adds
these six events over nonoverlapping weight ranges, preserving all prior
coverage. Its adapter
passes six tests and independent source review. No construction change,
encoder benchmark, promotion, commit, or push accompanies this update.

The [neighbor screen](theory/hd_field128_gamma30_d40_r5120_winner18_neighbors_float_20261010.json)
passes weights 4864, 5056, and 5184, but fails at 4480, 4608, 5376,
5632, and 6144. Those failures are positive upper bounds, not witnesses
against the construction. They delimit where fresh marker calibration
is useful; they do not provide coverage between sampled weights.

### Adaptive calibration beyond the first dense band

The next round applies the same strategy to previously unresolved weights
1280, 4352, 4480, and 5888. These are proof-parameter changes only. The
gamma30 construction, setup distribution, distance target, and source
multiplicities are unchanged.

The fitter evaluates every matching profile, selects a largest term, and
adds a positive coefficient marker fitted to one of that term's matrix
entries. It then evaluates the complete profile sum again. Each added
marker is stored as an exact rational; subsequent directed runs regenerate
coefficients rather than importing FLOAT bounds. A selected subset of
markers is accepted for further experiments only after a fresh complete
replay reproduces its parent score and regional-array hash.

| Original active-group weight | Initial FLOAT log bound | Final FLOAT log bound | Retained markers |
| ---: | ---: | ---: | ---: |
| 1280 | +8594.2144 | -364.9246 | 74 |
| 4352 | +10833.7485 | -7425.2414 | 44 |
| 4480 | +34806.2856 | -6290.0341 | 35 |
| 5888 | +127609.6634 | -4534.6161 | 17 |

These rows are complete numerical screens at the listed weights, not
directed certificates or interpolated ranges. The lower-half adaptive
adapter adds a finite positive chart for zero full-pair count and passes
seven tests. At weight 1280, seven adaptive additions close the FLOAT
screen without changing its initial z or alpha. At weights 4352 and
4480, an initial fixed set of target profiles is followed by one adaptive
addition. At weight 5888, 17 additions are followed by a bounded search
over z and alpha; all 40 completed trials retain the full profile domain.

The independent saved-result audits for
[4352](theory/hd_field128_gamma30_d40_r4352_winners_replay_20261010/independent_consistency_audit.json),
[4480](theory/hd_field128_gamma30_d40_r4480_winners_replay_20261010/independent_consistency_audit.json),
and [5888](theory/hd_field128_gamma30_d40_r5888_winner_replay_20261010/independent_consistency_audit.json)
check exact marker lineage, source pins, original-weight domains, and
the recorded score arithmetic. They do not repeat numerical objectives.
The [1280 lineage check](theory/hd_field128_gamma30_d40_r1280_adaptive_lineage_check_20261010.json)
separately records its exact proposal conversion and fresh replay.

Fresh directed coefficient banks turn these numerical proposals into
29 complete, independently audited weight bands. The additions cover
1268--1472, 4311--4831, 5216--5279, and 5345--5983: 1,429 previously
unresolved weights. Each integer weight is evaluated explicitly. The
audits check exact counts, integer partitions, source bindings, and saved
ceilings; they do not rerun coefficient generation or regional transfers.
The [dense-band manifest](theory/hd_field128_gamma30_d40_r5888_winner17_directed256_20261010/audited_band_manifest.json)
records the ten bands covering 5345--5983 and excludes the failed trial.

The neighboring weight 5344 exposed a limitation of the 5888 witness.
One unresolved profile had 1249 full pairs, 2846
singleton pairs, and one empty pair in each alternating matching.
Its directed log ceiling is -24.264403800332, short of the desired
40-bit whole-code margin even before combining other events. This is
a limitation of the current bound, not a low-weight codeword. The
[5344--5407 trial](theory/hd_field128_gamma30_d40_r5888_winner17_directed256_20261010/band_5344_5407/summary.json)
remains partial and contributes no indexed event. A separate complete
[5345--5407 run](theory/hd_field128_gamma30_d40_r5888_winner17_directed256_20261010/band_5345_5407/independent_audit.json)
closes that suffix with 92.315734102991 bits of combined margin.
Fresh calibration at weight 5312 subsequently closes this gap, as
recorded below. The failed 5888-witness trial remains preserved.

A [bounded weight-1200 attempt](theory/hd_field128_gamma30_d40_r1200_from1280_adaptive_float_20261010.json)
retains the weight-1280 z and alpha but recomputes every numerical bound.
Its 32 marker additions improve the complete FLOAT log bound from
+2335.46845 to +1213.29136 without closing it. Each of the 33 completed
sums retains all 361,201 profiles; the final regional entries use no
positive fallback. This attempt contributes no certificate. The subsequent
weight-1200 calibration changes z and alpha and closes neighboring bands.

A failed 4384--4447 trial using the 4480 witness exhausted disk space
while writing intermediate checkpoints. Its
[resource-stop record](theory/hd_field128_gamma30_d40_weights4384_4447_midpoint_band_20261010/resource_stop.json),
context, and first three checkpoints are retained. Only that run's
21,196 incomplete, reproducible checkpoint files were deleted, freeing
1,764,090,003 bytes. No certificate or prior research was removed. This
run contributes no event to the index. Subsequent covers cap the number
of evaluated cells at 512 or 2048 so failed trials remain bounded.

### Closing gaps with fresh proof parameters

The next round retains the same gamma30 encoder and setup distribution.
It changes only coefficient markers and the scalar proof parameters
z and alpha. Each numerical screen still includes every feasible pair
of matching profiles at its original active-group weight.

At weight 5312, two marker additions turn a positive FLOAT log bound
of 5095.1452 into -9615.5753. Fifteen markers suffice to reproduce the
complete score and regional-array hash. A fresh directed bank closes
[5280--5343](theory/hd_field128_gamma30_d40_weights5280_5343_midpoint_band_20261010/independent_audit.json)
and [5344](theory/hd_field128_gamma30_d40_weights5344_5344_midpoint_band_20261010/independent_audit.json)
with respective margins of 5106.1486 and 12650.5443 bits. These results
include the formerly troublesome one-empty-pair profiles. They resolve
the gap without increasing state or adding encoding operations.

At weight 1200, a fresh geometric proposal interpolates positive proof
parameters from the 1000 and 1280 witnesses, not their bounds. Seven
adaptive additions then give a complete FLOAT log bound of -273.2088.
A fresh replay retains 71 markers. Eight complete directed bands cover
1192--1267; the [progress manifest](theory/hd_field128_gamma30_d40_r1200_directed_progress_20261010.json)
pins their independent audits. Partial wider trials are excluded.
Some bands stop splitting cells at 72 bits rather than 80; their
combined event must still meet the unchanged 54-bit admission threshold.

The [weight-3680 proposal](theory/hd_field128_gamma30_d40_r3680_winners_replay_20261010/independent_consistency_audit.json)
uses fresh geometric parameters from the 3072 and 4096 witnesses.
One adaptive addition improves its complete FLOAT bound from
+26924.9114 to -7453.2360. Sixty markers reproduce the score exactly.
The [weight-3456 calibration](theory/hd_field128_gamma30_d40_r3456_calibration_20261010/adaptive.json)
reuses these positive parameters but recomputes every bound at 3456.
Three additions reach -6207.4955; a fresh replay retains 30 markers.
These numerical successes select witnesses for explicit directed bands;
they do not certify intervening weights by interpolation.
Fresh directed banks now close 3424--3903 in eight bands, each with
an independent saved-data audit. The new lower and middle additions
1192--1267 and 2432--2847 appear in their
[separate manifest](theory/hd_field128_gamma30_d40_kernel_v20_band_manifest_20261010.json).

The dense gap needs separate witnesses near its edges. At weight 6400,
retuning z and alpha followed by adaptive markers gives -5549.7272 in
FLOAT. Twelve retained markers yield five complete directed bands over
6176--6495. Retargeting these parameters to 5984 and adding four markers
gives -4374.2366. A fresh 13-marker subset yields three bands over
5984--6175.

At 6677, 17 marker additions initially reach only -519.8006. Fifteen
more additions improve the log bound by less than 1.5 bits. A bounded
search over z and alpha, using the earlier 19-marker subset, reaches
-19354.9338. All 40 trials retain all 574,564 matching profiles. Fresh
directed coefficients then close three bands over 6496--6677. Together,
the [11 audited bands](theory/hd_field128_gamma30_d40_upper_gap_5984_6677_manifest_20261010.json)
cover every integer in 5984--6677. Intermediate witnesses
at 6080 and 6560 remain numerical proposals; they are not additional
certificate events. This comparison again favors retuning proof
parameters before changing the encoder.

The [original 3904--3936 trial](theory/hd_field128_gamma30_d40_weights3904_3936_midpoint_band_20261010/summary.json)
does not close with the 3680 witness. Its tree contains 454 evaluated
nodes and 47 unresolved leaves. Some unresolved leaves are individual
profiles, so further subdivision alone cannot repair this bound. For
example, at weight 3930 with 130 full pairs in each matching, the
directed log ceiling is -29.638876702061. A stronger bound is needed;
fresh proof parameters are one option. The partial run
contributes no probability event.

The two adaptive-marker helpers and the v20 index adapter pass all 19
targeted unit tests. Independent source review also checks the explicit
36-event assembly and its exclusion of partial runs. No encoder benchmark
was run in this round.

### Extending the recalibrated cover

The next round resolves the 3904--3936 gap without changing the encoder.
Recomputing all 3,876,961 matching profiles at weight 3936 and adding
one coefficient marker changes the complete numerical log bound from
+3604.9977 to -9243.6007. A fresh replay retains 35 markers with exactly
the same numerical regional bounds and score. Their directed bank closes
[3904--3936](theory/hd_field128_gamma30_d40_weights3904_3936_r3936_band_v21_20261010/independent_audit.json)
in two cells with a 5524.1740-bit combined margin. The original failed
trial remains preserved and contributes no event.

Weight 112 illustrates a different source of slack. Starting from the
weight-128 witness, 17 marker fits add 11 distinct markers but barely
change the positive numerical bound, +73.6936. A separate search changes
z and alpha while retaining the original 42 markers. Its 40 trials each
include all 3,249 feasible matching profiles; the best reaches
-1018.3548. A fresh replay retains 20 markers. Thus the successful chain
is 42 fixed markers, retuned parameters, then 20 retained markers; it
does not inherit the separate failed 53-marker attempt. The
[lineage audit](theory/hd_field128_gamma30_d40_r112_calibration_20261010/independent_consistency_audit.json)
checks these saved inputs and results without repeating their numerical
calculations. Fresh directed evaluation then closes
[104--120](theory/hd_field128_gamma30_d40_weights104_120_r112_band_v21_20261010/independent_audit.json)
with a 130.6521-bit combined margin.
Three further directed bands close 101--103, 121--127, and 129--139.
They use the fresh 20-marker bank and a separately rebuilt 42-marker
weight-128 bank. Together with the existing weight-128 event, they cover
101--139; the neighboring failed numerical probes are not included.

Two further calibrations close
[1473--1535 and 1537--1920](theory/hd_field128_gamma30_d40_kernel_v21_band_manifest_20261010.json)
in 26 independently audited bands, preserving the existing weight-1536
event. The weight-1473 bank retains 88 markers after fresh fitting and
replay. The weight-1792 bank starts with geometrically interpolated
positive proof parameters, then fits and retains 99 markers. No bound
is interpolated. Three intervals within 1701--1728 exceed the initial
512-node evaluation budget. New attempts using the same bank close them
in 610, 864, and 838 evaluated nodes, within a 2048-node cap. They retain
the same 72-bit cell and 54-bit band admission targets. All ten partial
attempts remain excluded.

Fresh banks near weights 2336 and 2944 close the remaining intervals
[2257--2431 and 2848--3007](theory/hd_field128_gamma30_d40_gaps2257_2431_2848_3007_manifest_20261010.json).
Their six independently audited bands cover 335 weights. The weakest
combined band margin is 87.0147 bits. The witnesses retain 71 and 68
markers, respectively; each bank recomputes the actual directed
coefficients before evaluating every requested weight.

A disk-space guard interrupted one further band before it wrote any
cells. Lossless NTFS compression of the proof directory recovered about
1 GB without deleting files or changing their contents. The failed
attempt is retained separately from subsequent retries. Unused in-memory
banks were also released. No encoder benchmark accompanies these proof
runs.
The two adaptive-marker helpers and the v21 index adapter pass all
19 targeted unit tests. All in-memory banks from this round were released
after their checks; saved proposals, successes, and failures remain.

### Closing lower ranges

The next round adds 152 weights in 15 independently audited bands:
72--88, 304--367, and 440--510. These are fresh directed calculations
for every included weight, not interpolation between successful probes.
The encoder, setup distribution, and 40% target are unchanged.

At weight 80, reusing a small marker set from weight 112 leaves the
numerical bound positive even after retuning z and alpha. Fresh occupancy
fits provide 76 markers; another parameter search reaches -1279.3134.
A fresh replay retains 38 markers with identical numerical regional
bounds and score. Its directed bank closes 72--79 and 80--88 with
combined margins of 302.1724 and 208.4989 bits, respectively. The
[lineage check](theory/hd_field128_gamma30_d40_r80_fresh_v22_independent_consistency_audit_20261010.json)
distinguishes this successful chain from earlier unsuccessful fits.

Fresh calibration at weight 320 gives a 45-marker proposal. Its rebuilt
directed bank closes 304--367 in eight explicit eight-weight bands.
The weakest band margin is 68.2859 bits. A separate 51-marker bank
closes [440--510 in five bands](theory/hd_field128_gamma30_d40_kernel_v22_band_manifest_20261010.json).
Neither result transfers its bound to neighboring weights.

The independent checker needed a provenance correction: some witness
ancestors link to Python source files, not JSON, and some group helper
hashes in nested maps. Separately named checker versions now distinguish
those formats, rehash every linked source, and retain the original
count, partition, enclosure, and union arithmetic. Old checkers and
receipts are unchanged. These audits do not regenerate coefficients,
regional matrices, or cell transfers.

At the version-22 checkpoint, several further proposals passed complete
numerical screens but had no new directed certificate. Replayed proposals at weights
68, 96, 152, 192, 224, 288, and 376 are retained in the theory directory;
their log bounds range from about -1128 to -297. The kernel manifest
also preserves the passing 400, 576, and 640 proposals. The upper-gap
work retains passing proposals at 960 and 1168; a fresh weight-913
screen remained positive. None of those uncompleted ranges was admitted then.

Low disk space interrupted further checkpoint writes. Those historical
producers used a 750-MiB free-space guard, and final bookkeeping used a
600-MiB reserve. Neither guard is sufficient for new numerical work:
the [current storage policy](AGENTS.md) requires 20 GiB of free local space,
including planned output, and bounded remote scratch for bulk work.
All owned in-memory banks were released; no proof files were deleted.
No encoder benchmark was run during that stop.

### Checks during the storage stop

A [saved-proposal audit](theory/hd_field128_gamma30_d40_root_edge_v22_independent_consistency_audit_20261010.json)
checks the six calibration chains at weights 68, 96, 152, 224, 288, and
376. Their exact positive witnesses, source bindings, saved numerical
replays, and neighbor-screen metadata agree. The audit does not regenerate
coefficients or evaluate the numerical objective. These inputs are ready
for directed runs over 65--71, 89--100, 140--164, 216--240, 273--303,
and 368--383. All 116 weights remained uncertified at that checkpoint;
the subsequent remote bands below cover 100 of them. Only 368--383 remains open
among those six proposed ranges.

The interval 165--179 needs fresh probes or calibration. In particular,
the saved weight-192 witness gives a positive numerical bound at 168,
although its probes at 180, 192, and 208 pass. Neither its successful
center nor the weight-152 endpoint probe at 164 certifies this bridge.

The [weight-913 replay](theory/hd_field128_gamma30_d40_r913_gap_mark_search_20261010/replay.json)
had a different limitation. Its largest term came from 410 full pairs
in each matching: each matching has 3593 empty pairs, 93 singleton pairs,
and 410 full pairs. This term contributes about 93.4% of the numerical
sum and has a log bound of +11.3198 after the outer factor. The full
bound is +11.4184. These are failed upper bounds, not low-weight words.
No selected regional entry uses the numerical fallback. The next step
is to refit markers at the best saved z and alpha, beginning with this
profile, then reevaluate the complete objective. Subdivision alone
cannot repair this individual-profile bound.
Exact support counts alone cannot repair it either with these regional
caps held fixed. The midpoint count exceeds the exact count by at most
a factor of 94 for this profile, or about 6.555 bits. Removing that
entire loss would still leave this term's bound positive, near +4.765.

The partial directed trees near weights 640 and 1168 instead retain
coarse matching-profile cells. Relaxing the 72-bit cell target does not
close them from saved receipts. In the
[interrupted weight-640 retry](theory/hd_field128_gamma30_d40_weights640_647_midpoint_band_nodes2048_v22_20261010),
even the best saved partition has an individual cell bound above
2^430. Thus no union of that partition can meet the whole-code target.
The retry reached only 710 evaluations before the storage stop, not
its 2048-evaluation limit.

For the [1145--1191 retry](theory/hd_field128_gamma30_d40_r1168_gap_directed256_20261010/retry2048/summary.json),
a bottom-up check takes the smaller of each saved parent bound and
the sum of its two child bounds. Missing child receipts contribute an
unbounded value. The best complete partition still has log bound
+959.1261. Its dominant cell groups odd weights 1145, 1147, and 1149,
with deficit ranges 38--74 and 0--74; both children lack evaluations.
These checks reuse saved directed ceilings but regenerate no regional
coefficients. They admit no new event.

Completed sibling splits in both trees reduce the envelope much more
when splitting a deficit axis than when splitting the weight axis:
typical saved improvements are hundreds of bits versus a few bits.
This supports prioritizing deficit refinement in the next bounded runs,
with separate calibration for weak lower endpoints. It does not prove
that the remaining cells will close. No construction change is justified
by these partial trees alone.

### Bounded remote continuation, versions 23 and 24

Version 23 admits the completed pilot bands 65--71 and 89--100. Its
independent check preserves all 242 version-22 entries and reaggregates
their union with those two events at 384 bits.

Version 24 adds four fresh directed bands. Each uses the original weight
in its regional normalization and outer factor; no interpolation is used.

| Newly covered range | Weights | Directed union margin |
| --- | ---: | ---: |
| 140--164 | 25 | 92.814830949126 bits |
| 180--208 | 29 | 87.196379333826 bits |
| 216--240 | 25 | 84.007480278738 bits |
| 273--303 | 31 | 69.783621943004 bits |

Each saved tree passes an independent typed audit of source bindings,
integer counts, complete partitions, checkpoint transitions, and the
384-bit union of saved leaf ceilings. These audits do not regenerate
coefficients, normalized regional matrices, or cell transfers. The
four events cover 110 new weights; the last band required 528 cell
evaluations after its earlier 512-node attempt stopped incomplete.

The version-24 368--383 attempt stopped incomplete after 2048 evaluations, with
896 accepted leaves and 258 pending leaves. It contributes no event to
the index. Its normalized-width splitting rule spends early work on
the weight axis despite much wider deficit intervals. A different
partition strategy may reduce this enclosure loss; this is not evidence
of an actual low-weight codeword or a failing individual-profile bound.

At weight 913, fifteen additional marker fits initially improve the complete
FLOAT log bound from +11.418359 to -31.043655. The selected witness
retains 79 markers. A fresh no-fit replay matches its numerical score
and regional-table hash over all 208,849 profiles. The largest term now
has 357 full pairs in each matching, rather than 410. This remains a
numerical screen, short of the 54-bit band-admission target; it is not
a directed certificate. Retuning the global positive proof parameters
is the next step before another directed attempt.

Completed results were imported in a 1.28-MiB archive; incomplete trees
remain remote. The [archive receipt](theory/hd_field128_gamma_remote_archive_v24_20261010.json)
identifies a verified durable copy of the remote execution snapshot,
including its sources and partial attempts. Older local ancestry remains
required for complete provenance. New runs must respect the storage
policy rather than relying on historical helpers' smaller guards.

### Compact subdivision and parameter retuning, version 25

The new controller changes the subdivision order, not the numerical bounds.
It splits the interval with the largest current absolute width, preferring
the two deficit axes on ties. The fresh 368--383 run completes after
802 evaluations with 402 accepted leaves. Its directed union margin is
70.81199658695 bits. The independent typed audit checks the complete tree,
803 checkpoint transitions, source bindings, and a fresh 384-bit union.
It does not regenerate coefficients or cell transfers.

Checkpoints no longer repeat the pending frontier; the final summary
retains the entire tree and frontier. This run produces 3,347,197 bytes,
versus 16,692,884 bytes for the previous incomplete attempt. Both the
different subdivision order and compact checkpoint format contribute to
the reduction. Neither change alters the encoder.

At weight 913, retuning z and alpha improves the full FLOAT sum from
-31.043655 to -55.977219. Eight additional marker fits reach -58.945542.
The pruned witness contains 83 markers; a fresh no-fit replay agrees over
all 208,849 matching profiles. These are numerical screens, not certificates.
The subsequent fresh directed tree stops at 4096 evaluations with 1133
accepted leaves and 1831 pending leaves. No evaluated singleton failed;
the remaining intervals have widths at most three on either deficit axis.
This establishes a resource-limited incomplete replay, not an adverse
individual-profile bound. The subsequent direct enumeration below resolves
that question without generating a larger tree.

The [completed-results archive](theory/hd_field128_gamma_remote_archive_v25_20261010.json)
preserves the compact band and FLOAT retuning results. The partial
weight-913 tree stays remote in a
[verified durable archive](theory/hd_field128_gamma_r913_partial_archive_v25_20261010.json).
It consumes no additional local tree storage and contributes no indexed event.

### Direct profile sum and three additional bands, version 26

The weight-913 replay uses the same 83 exact markers but replaces cell
envelopes by the existing directed sum over all 208,849 matching profiles.
It certifies 58.945660269253 bits. The independent saved-data audit checks
every integer midpoint and final count, all 457 external row checkpoints,
the reconstructed metadata context, and the 384-bit row/prefix sum.
The audit does not regenerate coefficients or transfer matrices.
Portable proposal links are resolved from the theory directory, exactly
as in the producer; the original receipt bytes are not rewritten.

The same widest-interval subdivision strategy also closes three bands
using existing exact witnesses. Every original weight is normalized
separately; bounds are not interpolated between sampled weights.

| Newly covered event | Weights | Directed margin | Evaluated cells |
| --- | ---: | ---: | ---: |
| 385--439 | 55 | 70.952761046342 bits | 306 |
| 576--639 | 64 | 70.725681362784 bits | 624 |
| 640--703 | 64 | 68.879257798808 bits | 1296 |
| 913 | 1 | 58.945660269253 bits | Direct profile sum |

Each band passes an independent typed audit of its integer partition,
checkpoint transitions, source bindings, and positive directed union.
The new index adds these 184 weights to version 25, retaining every prior
event. Including 368--383, this continuation adds 200 weights to version 24.
This is a proof improvement for the unchanged construction, not an
encoding optimization or a complete 40% certificate.

Numerical jobs ran serially on Peach and exited normally. The direct
weight-913 output is 1,183,012 bytes; the earlier incomplete subdivision
tree is 20,080,344 bytes and remains remote. The four returned result
archives plus uncompressed contents total 15,442,580 logical bytes before
deduplicating source files already present. Local audit and index output
is separately bounded. Durable copies are recorded in the
[weight-913 archive receipt](theory/hd_field128_gamma_r913_midpoint_archive_v25_20261010.json)
and [three-band archive receipt](theory/hd_field128_gamma_completed_band_archives_v26_20261010.json).
No new encoder timing, paper change, library promotion, commit, or push
is part of this continuation.

### Other retained proof refinements

A [state-31 variant](HD_FIELD128_GAMMA31.md) has a reviewed algebraic
description, but no directed certificate or encoder timing. It greatly
improves the old weight-512 FLOAT bound, yet barely changes weight 768.
The count and marker refinements therefore take priority over adding
state. The state-30 certificates remain attached to their original
construction.

The FLOAT and directed weight-1024 bounds differ by about 22.983 bits.
An independent diagnostic checks all 65 local caps, selected alias
entries, normalization, matrix order, profile counts, and source pins
without finding a mismatch. The difference is not yet fully localized;
rejection of small FLOAT alias residues is a possibility, not an
established explanation. The directed 256- and 384-bit results and their
checkpoint reaggregations agree. Evidence is retained in
`theory/hd_field128_gamma30_r1024_float_directed_differential_20261010.json`.

An independently reviewed [two-half Gram-rank refinement](HD_FIELD128_HALF_RANK_RETURN.md)
can further tighten the same construction's return cap. It accounts for
the actual cross-rank between inactive output constraints and remaining
active inputs. It is an algebraic result only: no new coefficient table
or numerical certificate uses it yet. In particular, the current maps
do not have generic cross-rank; same-half submatrices can have rank at
most four at state 30.

The independently reviewed [adaptive-cell coverage lemma](HD_FIELD128_ADAPTIVE_CELL_COVER.md)
provides a way to avoid a quadratic matching-profile loop at
every middle weight. It combines directed numerator maxima, convex
normalization bounds at every cell vertex, and exact support-count
projections. Singleton-free matching faces are kept separate to retain
their contracted-cycle counts. The [exact geometry](theory/hd_field128_gamma_cell_geometry.py),
[support counts](theory/hd_field128_gamma_cell_counts.py), and
[directed cell checker](theory/hd_field128_gamma_cell_certificate.py)
and the [bounded cover driver](theory/hd_field128_gamma_cell_cover.py)
are now implemented and independently reviewed. All 30 supplied tests
pass. A separate exhaustive audit checked 8,178 small boxes against
actual supports. These checks validate the machinery, not a completed
middle-weight cover.

The first fresh directed probe rebuilt the all-pair bank at 256 bits.
At weight 1024, the profile `(a,b)=(310,511)` receives a cell bound of
`2^-97.943290231451`. The box with `a=300..320` and the same `b` is
uninformative, as is the box with weights 1023--1025 and fixed `(a,b)`.
These are failures of coarse envelopes, not low-distance witnesses.
The six [probe receipts](theory/hd_field128_gamma30_d40_cell_probe_20261010/bank.json)
retain the rebuilt bank metadata and each cell's exact scope. They
overlap and are not aggregated into a coverage claim. The selected
profile is not established as the worst profile under the directed bound.

The narrow singleton count near this difficult profile motivates a
different coordinate choice for future weight bands. Fix the weight
parity `epsilon`, and write `r=2u+epsilon`, `a=u-d`, `b=u-e`.
Then the singleton counts are `2d+epsilon` and `2e+epsilon`, independent
of `u`. This aligns a cell with the near-boundary profiles instead of
allowing the singleton count to vary as the weight changes. The affine
change preserves the coverage lemma's convexity argument. Its separate
[geometry and counts](theory/hd_field128_gamma_deficit_cells.py) pass seven
tests and an independent audit of 5,916 small boxes. No numeric coverage
adapter or actual distance certificate uses these coordinates yet.
One proposed use streams the exact factorial normalization along each
fixed-deficit slice, looking up the actual alias numerator at every integer
point. A shortlist of markers selected at an anchor remains valid
elsewhere, although it can weaken the bound. Maxima would be taken only
after multiplying the numerator by its normalizer. This band method is
not implemented or used in any coverage claim.

The initial fixed-marker cover at weight 1023 stopped at its 300-second
budget after 993 cell evaluations: 35 accepted leaves and 924 pending
leaves. Its [partial tree](theory/hd_field128_gamma30_d40_r1023_cell_cover_20261010/summary.json)
passes a fresh partition audit but provides no selected-weight bound.
This method is not yet an efficient replacement for exact profile
summation. A refinement [normalizes regional bounds before grouping
profiles](HD_FIELD128_REGIONAL_TABLE_COVER.md), so the enclosure need not
separately maximize a marked coefficient and its inverse normalizer.
Its first bounded pilot evaluated 621 cells in 180 seconds: 156 accepted
leaves and 310 pending leaves. This is again only a partial cover.
The [guarded batch evaluator](theory/hd_field128_gamma_cell_pinned_batch.py)
now checks the complete evaluator and coefficient-source closures before
and after a cover run, instead of repeating those checks inside every
cell evaluation. Standalone cell checks remain unchanged. Intermediate
checkpoints are conditional on successful final source verification.
The geometry, counting, numeric, coverage, and batch modules pass all
52 supplied tests; the batch change does not alter their arithmetic.

The guarded replay then completed weight 1023 in 9.0697 seconds,
using 1,287 cell evaluations and 644 accepted leaves. Its directed
failure margin is 61.180906663845 bits. The independent
[coverage audit](theory/hd_field128_gamma30_d40_r1023_batched_regional_cover_audit_20261010.json)
checked the disjoint partition of all 262,144 profiles, exact support
counts, source pins, and the final union at 384 bits. The audit reused
the saved coefficient-derived bounds; it did not recompute the coefficient
bank. The stronger exhaustive 79.586-bit receipt remains in the partial
index. These are proof-computation times, not encoder timings.

The fixed multiplication by Y passes standalone scalar and SIMD checks
on Peach. It uses a byte shift and sparse reduction, not a general field
product. The [native check receipt](theory/hd_field128_mul_y_native_checks_20261010.json)
pins matching local and remote sources. This isolated check predates the
gamma30 whole encoder and timing above; no timing follows from the instruction count.

An [optional random-global-gamma variant](HD_FIELD128_RANDOM_GAMMA32.md)
extends the principal-matrix argument to state 32 with a separate setup
failure term below 2^-62. It has no numerical distance certificate and
loses the fixed generator's cheap multiplication rule.

## Arithmetic prototype

The isolated prototype lives in
`spin/experiments/char2_field_spin/field128_lift/`. Its initial native
polynomial basis uses X^128+X^7+X^2+X+1 and an explicitly verified AES
subfield embedding. Portable and four-record CLMUL checks pass on Peach
with GCC 15.2 and `-O3 -march=znver4`. Checks include all 16,384 binary
basis products, all 65,536 embedded-byte products, basis round trips,
independent SIMD coefficients, aliasing, alignment, and buffer guards.
This is arithmetic validation, not a full encoder or timing.

A second representation has passed native scalar and SIMD checks: GF256[Y] modulo the
irreducible polynomial Y^16+Y^5+Y^3+0x32. This permits the existing byte
GFNI fixed maps without runtime basis conversions; random coefficients
still range over the whole 128-bit field. It changes representation,
not the sampled construction or distance proof. Its vector kernel was
checked against an independent polynomial oracle, including every binary
basis product and full-field linearity. The full encoder now passes a
dense, independently implemented scalar oracle at K=8192 and K=16384,
including cross-tile and cross-step boundaries, buffer guards, unaligned
access, overlapping input/output buffers, and both output-store policies.
Its restriction to GF256 coefficients agrees with the retained byte
encoder. A non-subfield scalar test verifies full-field linearity.
The lifted proof tools also pass 53 local unit tests, and the portable
field arithmetic passes 13 tests.

## First matched timing

The first literal lift is substantially slower. On Peach CPU 2, GCC 15.2
with `-O3 -march=znver4`, serial 101-sample batches use the same payload,
store policy, prefetch policy, warmup count, and measurement harness.
Setup and allocation are excluded; output is hashed once after each
batch, not between timed calls.

| Measured work | GF256 control | Literal GF(2^128) |
| --- | ---: | ---: |
| Full encode, first batch | 3.604889 ms | 22.641665 ms |
| Full encode, reverse-order repeat | 3.550417 ms | 22.831610 ms |
| Outer alone | 1.431953 ms | 16.750518 ms |
| Routing and inner alone | 1.731603 ms | 5.697395 ms |

Standalone components have their own warm cache state, so their medians
must not be added to reconstruct a full-encode time. The historical
3.798160-ms reference is preserved rather than overwritten by this new
harness's control. The full-field coefficient payload is 72 MiB, versus
5.5 MiB for the byte control; both have a 2-MiB packet route. The timing
receipt, commands, and measured-source hashes are in
`theory/hd_field128_first_matched_timing_20261009.json`.

The byte-compatible full-field multiplication is now the principal
implementation target, especially in the outer. These measurements do
not establish a distance improvement or justify promoting this design.

A separately compiled GFNI-only convolution control passes the same
whole-encoder checks but is slower: 23.840774 ms versus a fresh CLMUL
control of 22.719691 ms. Its isolated outer and inner medians are
17.819374 ms and 5.759832 ms. These are again serial 101-sample batches.
Keep the CLMUL implementation as the arithmetic control; replacing
carry-less multiplication with byte multiplication alone is not a win.
The measured CLMUL outer loop also spills temporary vectors, so its cost
is not explained by key traffic alone.

A source-level follow-up confirms that this path performs direct arithmetic
in the sparse byte-tower representation. It does not repeatedly convert
between field bases. Each four-by-four Toeplitz map uses nine complete
field products. One representation-preserving candidate would combine
unreduced products before reducing the four final outputs. Its nominal
GFNI count falls from 54 to 24 per four parallel maps, but wider
temporaries could worsen register spills. The version-34 experiment above
now implements and measures this candidate. A native-polynomial-basis
alternative also needs an efficient
isomorphism for the actual sparse tower; the existing byte-basis helper
uses different extension coordinates and cannot be substituted directly.

The [first-principles diagnosis](HD_FIELD128_LIFT_DIAGNOSIS.md) separates
this implementation cost from the distance obstruction. In the present
fractional-moment bound, rescaling alpha by log(q) leaves the leading
rank and route contributions nearly unchanged. A candidate stronger
emission/feedback map can change those ranks and return weights without
increasing the state dimension. That is a proposed construction change,
not an inherited certificate.

The final gaps were closed using bounded remote runs.
The new weight-172 and weight-730 calibrations have already closed their
bands in version 28; they change proof parameters, not the encoder.
At weight 548, two separate witnesses give complete FLOAT log bounds
of +210.186 and +222.507. Taking their minimum separately for each
matching profile reduces the summed log bound to +75.931, still vacuous.
Each witness must finish its own chronological transfer and outer factor
before taking this minimum. Mixing matrix entries from different marked
parameters is not justified. A third witness tuned to the residual
profile sum initially improves the bound to +61.747 but does not close it.
Restarting that search from the other baseline yields -464.070466.
A fresh numerical replay agrees. Only two witnesses contribute to the
minimum; each alone still fails. The
[profile-wise argument](HD_FIELD128_PROFILE_MIN.md) explains why these
complete bounds can be combined without changing the code.
A fresh directed calculation covers all 75,625 profiles at weight 548
and gives 464.072437820792 bits. Its independent audit passes and version 30
includes this weight; nearby numerical successes are not interval certificates.

A separate retuning from the weight-913 witness gives a complete
weight-936 FLOAT log bound of -616.363. Its fresh directed band
920--967 passes with 59.996073946746 bits in 1568 cell evaluations.
The independent remote audit passes; its
[compact audit](theory/remote_receipts_v29/audit_920_967.json) and
[durable archive reference](theory/remote_receipts_v29/archive_920_967.json)
are retained, and version 29 admits this band.
The adjacent 914--919 directed run remains incomplete at 4096 cells,
despite passing numerical screens at 914, 917, and 919. Direct profile
summation subsequently closes weight 914 with 82.227040360459 bits over
all 209,764 profiles. Its independent audit passes and version 30 includes it.
At weight 1145, additional markers improve the complete FLOAT screen from
47.031 to 49.694 bits. The directed replay certifies 49.694183364555 bits;
version 30 includes its independently audited result. The
whole-code target is 40 bits; 54 bits per event is a bookkeeping choice,
not a required property of the construction. Any smaller admitted budget
must be included explicitly in a fresh union over all covered events.
The [remote diagnostic archives](theory/remote_receipts_v29/float_archives_20261010.json)
preserve the exact proposals and numerical history without local bulk copies.
Prefer compact direct profile sums when subdivision leaves nearly
singleton cells; use grouped cells when large regions close together.
Retune these parameters when additional markers stop improving a complete
numerical bound; subdivision cannot repair a failing individual profile.

Keep the original weight in every regional coefficient and outer factor;
only support counts may use complementation. Retain every feasible
profile and distinguish complete directed bands from numerical screens.
Version 31 and its separate completion audit now close the 40% cover.
The next proof target is 45%. Performance work must address full-field outer products;
the gamma30 encoder is validated but remains much slower than the GF256 reference. Preserve the certified literal
lift and the faster GF256 reference as separate controls.
