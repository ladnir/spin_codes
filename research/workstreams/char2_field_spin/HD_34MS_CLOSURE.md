# The 3.4-ms high-distance search (parked)

## Wrap-up, 2026-10-09

The user accepted the certified 3.8-ms design and requested that we wrap up
unless a faster design could be closed quickly. We are therefore parking
this search, not declaring its original target achieved. The selected
research deliverable is the [paired-route checkpoint](HD_SHIFTED_CLOSED_CERTIFICATE.md):
K=2^20, rate 1/2, 35% symbol distance, 55.092738508223 bits of ideal-setup
failure margin, and a measured median of 3.798160 ms.

The complete cover checker passed again with 384-bit directed aggregation
of the pinned, previously verified bounds. It covers every active-group
weight from 1 through 8192 and all 11,465,830,401 feasible matching profiles.
This check authenticates the retained evidence and recomputes coverage and
aggregation; it does not freshly recompute every coefficient bound.

There is no established quick closure for the cheaper designs. Mixed row
counts, additional source-pairing patterns, and other active-group counts
remain unresolved. The shear also lacks a native timing. Its predecessor's
3.442226-ms timing must not be assigned to the modified encoder.

### Preserved partial work

- The [scheduled cheap outer](HD_SCHEDULED_CHEAP_OUTER.md) and
  [Renyi proof record](HD_LAMBDA_RENYI.md) retain the unchanged Hlambda
  implementation, double-transvection candidate, and column-shear candidate.
  Their certificates cover different ensembles and must not be combined.
- For the shear, the shared-scalar cancellation lemma is implemented in
  `theory/hd_lambda_shear_shared_ratio.py`. It bounds marked expressions
  under a positive common measure, not the entire code's failure probability
  by itself. All 560 three-row choices were checked with the earlier exact
  method. Both results retain their recorded source-class restrictions.
- The four-row shared-scalar scan stopped after 1,152 of 1,820 labelled sets.
  The three preserved receipts are
  `theory/hd_lambda_shear_shared_a4_batch0000_0384_float.json`,
  `theory/hd_lambda_shear_shared_a4_batch0384_0768_float.json`, and
  `theory/hd_lambda_shear_shared_a4_batch0768_1152_float.json`.
  Each batch is complete, but the complete four-row class is not covered.
  These are FLOAT checks at r=18 and Fmax<=1, not directed certificates.
- `theory/hd_lambda_shear_shared_ratio_extended.py` and
  `theory/hd_lambda_shear_shared_extended_control.py` preserve the reviewed
  exact single-column machinery through six rows and paired degree twelve.
  Fourteen tiny/mock tests passed; no production a5/a6 gate was run.
- `theory/hd_lambda_shear_profile_measures.py` preserves an adapter for
  all source-pairing profiles, including feasible active-subset constraints.
  Eight tiny/synthetic tests passed. It has not had independent review or
  an actual regional-bank run and supplies no new certified coverage.

If the search is reopened, first review and test the full-profile regional
measure before resuming row-set enumeration. Then close mixed arities and
other active-group counts, verify the complete union with directed
arithmetic, and measure the actual modified encoder. These are research
steps, not a final verification pass. For the accepted design, the next
useful task is packaging the existing implementation and certificate for
reproduction, without changing the construction.

## Target and preserved checkpoints

The original goal was a complete precomputed forward encoder at K=2^20,
rate 1/2, over 128-bit records, taking approximately 3.4 ms on Peach CPU 2.
The required guarantee is at least 35% symbol distance with at least
40 bits of failure margin over the ideal setup distribution.

The [paired-route checkpoint](HD_SHIFTED_CLOSED_CERTIFICATE.md) remains
the certified reference: 3.798160 ms and 55.092738508223 bits. Its original
72.3975-bit reference is also preserved. The
[G16 two-sum implementation](HD_TWO_SUM_GF16_NATIVE.md) takes 3.397943 ms,
but its whole-code certificate is open. Timings come from different
matched series; they do not establish a new matched speedup ratio.

This search followed the previous
[structural-search pass](HD_STRUCTURAL_ITERATION.md). It does not reopen
the same failed envelope fits. A partial certificate, a cheaper operation
count, or a weaker distance target does not satisfy this goal.

## Three parallel approaches

1. Preserve correlations across each complete RS group when counting the
   two cancellation sums. A new syndrome-rank lemma can bound exceptional
   setup choices and determine fiber sizes. It is not a distance bound:
   message coordinates conditioned on a fiber remain correlated.
2. Strengthen the G16 outer cheaply. The first candidate replaces its
   final scalar diagonal by independently sampled invertible Toeplitz2
   maps on pairs, using a fixed schedule of different pairings. Prove a
   useful local comparison before implementing it. The schedule must
   remain fixed across columns in the subsequent proof.
3. Reduce implementation cost without changing the certified encoder.
   A compact representation of the Toeplitz4 keys stores seven prepared
   coefficients instead of nine. Two linear identities reconstruct the
   omitted coefficients. Test the exact map before measuring complete
   encoding against the unchanged implementation.

## Validation and stopping rules

Keep exact finite-field tests separate from numerical distance bounds.
Before trusting a new bound, identify what is fixed, what remains random,
and whether a comparison is pointwise or averaged. In particular,
concavity does not permit replacing a powered fiber expectation by the
expectation of the powered integrand in the upper-bound direction.

Native candidates are opt-in research variants. Run correctness checks
against an independent encoder and matched full-encoding timings, one
benchmark at a time. Retain source and binary hashes, parameters, and
checksums. Shared-host load may affect timings.

Completion requires whole-code coverage, independently checked directed
aggregation, and matched performance measurements. No library promotion,
paper edit, commit, or push is part of this goal.

## Initial status

The [syndrome-rank lemma](HD_TWO_SUM_SYNDROME_RANK.md) gives a simultaneous
exceptional-setup bound below 2^-46.6757. Its exact tests pass and its rank,
gauge, and row-Hölder arguments have been independently reviewed. This is
not a distance certificate. Six bounded quiet-branch comparisons show that
the main opportunity is delaying mask fractionalization, not replacing an
old approximate fiber size by its exact size. The fiber bound gains 62--77
local bits per group over the column-first control, but remains 15.5--26
bits worse than the old mean-first control. These are weighted counting
diagnostics, without the full distance-bound factors.

The [pair-refresh law](HD_TWO_SUM_TOEPLITZ_PAIRS.md) has exact small-field
tests and a positive fixed-set kernel. It disperses an isolated quiet mask,
but the bounded one-column averaged gate does not establish a useful
whole-code improvement. A [pre-mixer variant](HD_TWO_SUM_PRE_TOEPLITZ_PAIRS.md)
has a valid comparison invariant
under independent input-diagonal scalings. That permits the RS measure
comparison without adding those diagonal multiplications back to the
implementation. Its ordinary scalar-moment gate also gives little gain;
this does not settle the support-resolved fractional comparison.

The [compact-key implementation](HD_TOEPLITZ_COMPACT7.md) passes native
basis, full-encoder, and K20 oracle checks. It reduces hot keys by 1 MiB,
but three matched runs give 3.678 ms versus 3.686 ms: effectively a tie.
It remains default-off. No new whole-code certificate or meaningful
timing improvement is claimed.

## Next decisive gates

First evaluate the original G16 construction with an entire outer column's
mask sum retained before the fractional power. The syndrome note now states
the coefficient-level interface. Any geometry partition must describe the
placement of the entire fixed source set, common to every mask. Moving a
power through an expectation in the wrong direction, or changing the route
event for each mask, would invalidate the comparison.
Independent review confirms the coefficient/state expansion and the
independent, identically distributed route-block step. Fixed-route caps must
be constructed before route averaging; the existing averaged regional banks
may enter the Jensen mean, not the fixed-route fractional power.

In parallel, evaluate the proposed pre-H pair layer's complete labelled-mask
count before taking its fractional power. A positive pair polynomial reduces
this to joint input/output pair types. Retain the fixed constituent set and
matching schedule across all sixteen columns. Ordinary weight moments alone
are not the relevant acceptance gate.

Run only a bounded sparse/dense diagnostic first. Build a full directed
cover, and implement/time the changed construction, only if that diagnostic
shows a credible path to the target. The exact-map Compact7 optimization is
not the primary speed route after its effectively tied timings.

Root independently reran 33 new exact/scalar tests across the compact-key,
rank/fiber, and pair-law helpers; all pass. Those checks are distinct from
whole-code directed certificate validation. The goal remains active.

## First correlation and pre-mixer results

The [column-late-alpha control](HD_TWO_SUM_LATE_ALPHA_QUIET.md) uses one
fresh unpowered regional bank at the frozen 18-group, arity-eight anchor.
Summing the complete column's quiet masks before the fractional power
reduces the selected majorant term from log2 2813.42 to 1203.04.
This is a substantial reduction, but the remaining term still prevents
this majorant from closing. It is neither an actual failure lower bound
nor a bound covering the complete quiet event. The next proof-only test
conditions on route events common to every mask of the fixed source set;
the [common-source partition](HD_COMMON_SOURCE_ROUTE_PARTITION.md) retains
the correct geometry of pairs spanning two emission blocks.

The [pre-H pair control](HD_TWO_SUM_PRE_TOEPLITZ_FRACTIONAL.md) gives a
different opportunity. For one fixed eight-row set, a matching entirely
across its boundary reduces the powered mass of output weights 1--4 by
52.56 bits and weights 5--8 by 96.34 bits. An entirely internal matching
gives no improvement. These are local counting comparisons, not margin
bits. The fixed sixteen-column matching schedule includes both cases.

An incidence count over matching orbits now supplies the missing
fixed-target sum over all row sets. It permits the same row set to be
retained across columns through Hölder, without substituting the earlier
fixed-row-set diagnostic for that sum. The next construction test combines
these weights with the existing fractional regional operators and keeps
the matching schedule in chronological order. It must first pass a bounded
homogeneous-cohort test before a complete cover or native implementation
is warranted.

Root independently reran ten tests for the pre-H fractional helper and
thirteen tests for the late-alpha and common-source helpers; all pass.
Both bounded numerical controls completed. Neither supplies a new
whole-code certificate, and no new timing claim is made.

## Source geometry changes the next proof step

The fixed-schedule pre-H comparison improves its selected full-state
majorant by about 830 log2 units, but remains far above the target. A
separate positive-component calculation restores the exact empty-state
branch. Its all-zero-boundary contribution still fails at each frozen
witness. Thus a much larger state-path calculation with the same bounds
is not justified. The pair modification is retained as a candidate, not
as a proved or timed improvement.

For the unchanged G16 code, retaining more of the unpowered occupancy
curve gives a smaller improvement. A positive two-atom moment envelope
reduces the original selected quiet contribution from log2 1203.04 to
1026.17. Its coefficientwise concentration argument applies to every
fixed eight-row-set vector; it does not independently choose row sets
in each column.

The more consequential change is restricting source geometry before
building that envelope. Let F be the larger full-pair count in the two
cyclic source matchings. For F=1, every regional subset has at most one
full pair. A fresh unpowered bank gives log2 R00(18)=-184.87 instead of
the unrestricted -132.37. With the exact F=1 source count, the late-column
control drops to log2 440.85 and its moment refinement to 295.01.
These remain positive contributions of an upper-bound calculation, not
failure probabilities or evidence of a low-distance codeword.

The [common-route diagnostic](HD_COMMON_SOURCE_ANALYTIC_GATE.md) supplies
a complementary improvement. At regional occupancy nine, conditioning
on no emission block receiving more than four potential source positions
gains 3.81 log2 units over unconditional Jensen. The event concerns the
entire source set and is shared by all masks. The next test combines this
conditioning with the F=1 moment bound, keeping the complete column sum
inside the fractional power. It remains scoped to the quiet, zero-boundary
contribution until the other syndromes and state paths are included.

All new numerical gates are bounded FLOAT diagnostics with exact small
tests and independent source review. They are not directed certificates.
No native benchmark ran during this proof pass.

## Whole-word ordering clears the selected obstruction numerically

A cheaper ordering change supersedes the proposed combined-flag test.
Fix the source set G, but retain all of its fixed row-set assignments
inside the conditional first moment. Sum the complete sixteen-column
input-weight tail before applying a fractional power. The same moment
polynomial bounds each fixed row-set vector, so its row-set multiplicity
belongs inside this power. Only the source-stratum count remains outside.
This ordering does not resample either the row sets or the route.

The saved F=1 control reproduces its old log2 term, +295.013579.
Delaying the whole-word tail and the row-set sum reduces that term to
+54.944642 at alpha=1/5. More importantly, the regional provider is already
unpowered: the direct first-moment endpoint alpha=1 gives
**-432.03056185400425**, including every source in the F=1 stratum and
every fixed arity-eight row-set vector for r=18.

This is a successful FLOAT check of one positive component of the
distance-bound majorant. It covers zero outer moments and zero state
at every regional boundary. It does not cover nonquiet moments, live
boundary-state paths, other arities, or other source strata. It is not
yet a directed receipt or a whole-code distance certificate.

Six exact or mocked order-control tests pass. Two independent reviewers
agree on the ordering, normalizations, and scope. The run reuses the
pinned F=1 bank and moment envelope; it builds no provider and fits no
new parameter. See [the ordering derivation](HD_TWO_SUM_LATE_ALPHA_QUIET.md).

The next priority is a complete arity-eight/F=1 class bound that includes
nonquiet moments and live regional boundary states. Keep the successful
whole-word ordering, then extend class coverage and perform directed
verification. The additional pair-mixer implementation, combined-flag
gate, and large labelled-floor calculation are deferred until that test.
The existing 3.398-ms encoder is unchanged.

## Nonquiet outer moments clear the zero-boundary control

The next control retains the same r=18, arity-eight, F=1 class, but includes
every outer-mixer syndrome. For a fixed eight-row set, a quiet column has
equal input and output weights, at most eight. Every nonquiet column has
output weight at least nine, whereas its input weight is at most eight.
Thus min(output weight,8) is a valid upper marker for input weight in both
cases. The necessary whole-word floor remains 1296.

Holder over the fixed row sets supplies the full output-mask weights.
This comparison counts each row-set choice once across sixteen columns;
it does not independently choose row sets by column. Under the same
R00 moment law, the quiet-only and full-outer controls both give FLOAT
log2 term **-88.13020516064898**. Their equality at displayed precision
does not mean that nonquiet terms vanish. This is a weaker quiet bound
than the preceding -432.03 control, but it extends to every outer moment
and still clears the requested 40-bit target for this positive component.

The saved control is `theory/hd_two_sum_gf16_mean_zero_gate_float.json`,
SHA256 `610eb126c6dd2c34bf8ad8a36ec447d372d5fed3321d05b216d2b54003e57cd2`.
Eleven exact or mocked helper/driver tests pass independently at the root.
The calculation took 0.182 seconds and reused the saved regional bank.

Live boundary-state paths remain unresolved. A coarse entrywise comparison
R(m)<=M(m)C loses 14699.75 log2 units when C is raised to the 256th power.
Its worst-case ratios combine the steep R00 decay with much flatter live
entries; this numerical failure is not evidence of a low-distance word.
The next calculation keeps separate positive moment measures for the four
state transitions and composes their chronological column matrices. It
also retains the requirement that a state-changing region be nonempty.
No native code, timing, or certified reference changed in this pass.

## Ordered paths isolate repeated activation and return

The separate-entry moment calculation now includes every regional state
path and every outer syndrome in the selected homogeneous class. Its
complete positive polynomial tail gives FLOAT log2 comparison
**+3497.359683472526**, so it does not close the class. The quiet-only
calculation agrees at displayed precision. The nonquiet remainder is not
the dominant issue in this comparison.

A positive-term diagnostic locates the largest tail coefficient at the
floor 1296. One largest column word has nine columns of degree 144 and
seven empty columns. In a full degree-144 column, almost all of the bound
comes from paths with eight empty 00 regions, four 01 transitions, four
10 transitions, and no 11 transitions. Each of the eight active regions
receives all eighteen groups. These are repeated activation/return
excursions, not sustained live-state propagation.

This diagnosis concerns terms in an upper bound, not actual low-distance
messages. In particular, the outer alone may make some such patterns
rare. A separate degreewise outer-count cap can test that loss without
changing the encoder.

Enforcing that every state-changing region is nonempty gives no measurable
improvement under the matched scalar comparison. A global convex choice
of the input-weight marker gives log2 X=1.109035132306877 and comparison
+3543.61424510912, with or without that restriction. The earlier adjacent-
coefficient marker fails because the polynomial has strongly separated
degrees; it is preserved as a negative control, not a useful witness.

The [moment and coverage record](HD_TWO_SUM_GF16_MEAN_COLUMN.md) indexes
these checks. They rule out empty-transition slack as the current repair.
The next structural candidate is to interleave two independent inner
state streams, retaining the arithmetic per packet but requiring two
returns when both streams activate. Packet placements between streams
are correlated and must be counted jointly. Neither its distance nor its
complete encoding time is established yet.

A complementary [common-source route lemma](HD_COMMON_SOURCE_ROUTE_PARTITION.md)
is now exact: with probability greater than 1-2^-63, every r=18, F=1
source set has at most twelve regions containing an emission block with
more than four potential source positions. The proof uses exact rational
coupon probabilities and counts source sets once; it assumes independence
of regional routes for a fixed source, not independence between sources.
This is a setup lemma, not a distance certificate. The repeated excursions
can still occur on routes satisfying that local occupancy constraint.

The degreewise outer-only cap has also completed. It counts fixed row sets
once in an unmarked column polynomial, without the marked calculation's
Holder relaxation. It improves the quiet full-state comparison from
+3497.36 to +2943.12, with dominant degree 1830 instead of 1296. Thus rare
minimum-weight outer patterns explain some loss, but removing that loss
does not close the bound. See the
[cap derivation and receipt](HD_TWO_SUM_LATE_ALPHA_QUIET.md).

The [two-stream proposal](HD_TWO_STREAM_INNER_PROPOSAL.md) is the recommended
next proof-first gate. It includes a joint paired-slot polynomial interface
and explicitly retains one-stream-only placements. No native speed claim
is made for it.

The root independently reran all 53 tests for the new live-state, full-outer,
matrix-moment, path-diagnostic, coverage, rare-route, and degree-cap helpers;
all pass.

## A matched two-stream control and a coupled alternative

The joint paired-route provider for two interleaved states reproduces the
saved one-state regional bank exactly at FLOAT precision. It retains one
global routing denominator and extracts the shared occupancy coefficient
before normalizing. Its exact three-state quotient also agrees with the
four-state product comparison. Fourteen small tests pass at the root, and
both the provider and the whole-class gate received independent review.

The frozen scalar-potential control does not improve the selected class.
Its log2 comparison is +4462.31 for one stream and +8109.76 for two streams.
These are failed upper-bound calculations, not failure lower bounds.
Incoming one-live states determine most of the two-stream occupancy curve;
the calculation does not force both streams to pay a return cost. The
[matched gate record](HD_TWO_STREAM_MEAN_GATE.md) preserves the parameters,
source pins, complete result, and scope. No native variant was run.

The next candidate couples the two states instead of merely interleaving
them. For each input packet X, emit X+A U and update
(U,V) to (V,M(U+V+C X)). It retains one emission, projection, and refresh
per packet, adding four ZMM XORs. The arithmetic count does not establish
its runtime; register pressure and the new dependence need measurement.

The [five-cone derivation](HD_COUPLED_TWO_STREAM_CONES.md) identifies a
necessary proof distinction. A cancellation can leave a small point mass
in one state. Treating that mass as uniform would lose its probability or
give an invalid return bound. An additional cone retains its total mass
until the next packet. The proposed comparison also uses a fresh bound
for two live states; the old inverse-refresh sharpening does not transfer
to that transition.

The next decisive test is one raw coupled-state regional bank followed by
a complete arity-eight/F1 class control. A positive result would justify
expanding class coverage and implementing a native candidate. It would
not itself close the whole-code certificate. The goal and both certified
reference constructions remain unchanged.

## Post-refresh coupling and fixed row-set chronology

Moving the companion-state XOR after the refresh gives
`(U,V) -> (V,M(U+CX)+V)`, still emitting `X+AU`. For a fixed input word
starting at zero, a first nontrivial return requires two independent
refresh hits. The exact three-packet condition and a longer first-return
bound are recorded in the [post-refresh note](HD_POST_REFRESH_COUPLING.md).
These unweighted statements are not substituted for a weighted distance
bound. Eight comparison families retain point spikes and conditional
flatness explicitly.

The one new regional bank contains all 36 histograms through occupancy 18
and one full pair, at the unchanged mark 691/700. Its fixed live-reference
compression passes the local tests and independent review. The complete
homogeneous arity-eight/F1 FLOAT gate is still positive, at
**+2282.476801465376**, compared with +2906.605466 for pre-refresh coupling.
Neither value is a failure lower bound or a distance certificate.

A saved-polynomial diagnostic reproduces the new tail exactly. Its largest
column word contains nine degree-144 columns and seven empty columns. Each
full column is dominated by eight empty zero-state regions, one cross-state
transition, and seven sustained-live regions. Alternating activation and
return columns put the eight active regions at opposite ends of the column.
For the dominant quiet weight-eight outer branch, the output row set equals
the fixed active constituent set. That set cannot switch in this way across
columns. The existing Holder comparison is sound but forgets this constraint.

The next gate therefore changes the proof contraction, not the encoder.
A whole-word group Minkowski inequality permits summing bounds with one
fixed row set across every column. It applies at a common scalar input-weight
marker; it is not a coefficientwise inequality after taking roots. The
quiet part will retain the locations of empty-region gaps under the same
positive comparison measure, while a separate flagged calculation bounds
the nonquiet remainder without subtracting two rounded totals. In particular,
the homogeneous row-set terms cannot each substitute a sharper empty-region
operator: that would change the common measure required by Minkowski and
could discard heterogeneous-row-set contributions. Whole-class coverage
and directed validation remain future obligations even if this gate passes.

No native code, performance measurement, library profile, paper, or certified
reference changed in this pass. The 3.398-ms timing still belongs only to
the original experimental G16 encoder, not to either coupled variant.

The positive nonquiet calculation subsequently completed. It reproduces
the old total but still gives **+2226.121421174038** for the nonquiet
component. Its small relative size is not a useful absolute bound, so it
cannot close the complement of a new quiet proof.

The full fixed-row-set helper instead uses a three-way positive Young
bound inside each column: quiet masks, full outside support, and one
outside hole. Full support uses the common moment matrix evaluated at
occupancy 18. A positive degree-one flag inserts the empty moment matrix
at exactly one outside position for the hole branch. All three column
matrices are added before the sixteen-column power and row-set Minkowski
sum. Exact small controls and independent review cover this construction.

The first marker from the old construction was poor for the new bound.
Three output marks were screened with the original Holder contraction:
99/100, 497/500, and 499/500 gave +1403.86, +2282.09, and +7601.90,
respectively. These are witness probes, not new configurations or claims.
A full fixed-row-set run at mark 99/100 and log2 activity marker 2.5
completed in 33.74 seconds, but still failed at +2292.32; its quiet
component was +2288.38. The unshifted full/hole matrices are cached so
subsequent witness choices do not repeat those dynamic programs.

At mark 497/500 and activity marker one, the next full fixed-row-set
calculation finished in 33.74 seconds. The quiet component is +517.63.
Four shared branch-weight choices give full bounds +1825.86, +1601.02,
+1503.52, and +2024.65. Thus retaining fixed row sets removes an identified
relaxation but has not closed even this selected class. The full/hole
Young domination now also has visible cost; the construction has not
been disproved. Both branch banks and all controls are preserved.

The final focused regression passes 54 packet, region, compression,
return-identity, moment, diagnostic, positive-remainder, and coherent-row-set
tests. The certified cover manifest is unchanged. Next work should assess
the remaining joint-branch/marked-emission loss or a performance-friendly
construction change, rather than treating these witness probes as a
whole-code certificate. No further envelope sweep is queued.

## Exact-map performance and a new outer candidate

The next pass prioritizes two concrete changes rather than more fits of
the same unsuccessful envelope. Setup-time route-mask decoding preserves
the certified construction and gives about 3.68 ms against a matched
3.81-ms control. Independent review, basis and dense native checks, and
ASAN/UBSAN support this implementation change. The existing complete
55.092738508223-bit certificate transfers by exact-map equality.
The [experiment record](HD_DECODED_ROUTE.md) also retains negative source
fusion and address-scaling results. The 3.4-ms goal remains open.

The new [bit-syndrome outer](HD_BIT_SYNDROME_OUTER.md) adds four XORs to the
G16 mixer without adding field products, keys, or passes. Its quiet kernel
has five independent affine-coordinate syndromes instead of two. Exact
counts improve the troublesome arity-eight quiet branch, but some pair
inputs now leave six holes. Old pair and nonquiet-support bounds therefore
do not transfer. The next proof step is a positive column bound covering
both branches, followed by a selected full-state gate. Local improvements
alone do not establish a whole-code margin or a native runtime.

The subsequent [native experiment](HD_BIT_SYNDROME_NATIVE.md) passes
basis, dense K20, and sanitizer checks. Its three full run medians are
3.335718, 3.410037, and 3.621712 ms, with the original one-state inner and
predecoded swap masks. A matched old-G16 control gives 3.350616, 3.284903,
and 3.478614 ms. This noisy comparison makes the stronger outer worth
testing at the target cost; it is not a distance claim.

The [positive kernel](HD_BIT_SYNDROME_KERNEL.md) now bounds every input
column for each fixed arity-eight set. Exact fiber counts retain outside
outputs and an information set of the quiet kernel. Omitted coordinates
have explicit max(1,t) factors, so the bound also covers moment atoms
above one. Exact singleton and pair terms can be restored through integer
count subtraction before adding their true masks. Four affine orbits
supply universal coefficient caps only; they do not collapse the final
ordered sum over labelled source sets. Fourteen tests and independent
review pass. No numerical contraction has yet used this new kernel.

The next full-state gate is now implemented and independently reviewed.
Its initial universal coefficient bound fails, dominated by nonquiet
words. Separating four exact coefficient types retains all labelled
row sets and avoids assigning impossible low-weight patterns to most of
them. One prescribed weaker output marker clears the quiet branch but
leaves long-lived state paths dominant. A single intermediate marker
z=993/1000, with X=1, clears the complete r18/arity8/Fmax1 stratum in FLOAT:
the bound is -58.17 bits with balanced branch weights. Recombining the
same cache with weights1/16,15/16 improves it to -251.46 bits.

This [selected-class result](HD_BIT_SYNDROME_KERNEL.md#refined-selected-class-result)
uses the unchanged measured 3.41-ms bit-syndrome encoder and fresh
original one-state regional curves. It includes both outer branches,
all boundary-state paths, and every fixed arity-eight row-set assignment.
It is not yet directed and does not cover the whole code.

An exact singleton-outside fiber calculation is also preserved. It forces
at least four inside nonzeros universally, or five outside32 exceptional
sets, improving a relaxation in the current bound. It was not needed for
the successful selected-class FLOAT result.

Fresh directed selected-class replays at 256 and 384 bits now agree:
r18, homogeneous arity eight, F<=1 has log2 upper bound -246.68569098238194.
The replay includes both outer branches and all boundary-state paths.
It uses newly computed regional bounds, checked exact dyadic atoms,
upward positive matrix arithmetic, and a high-precision final sum.
The measured implementation has not changed. Independent 512-bit final
checks pass for both receipts: all 12,870 cached word bounds dominate a
separate successive-row-product calculation, and the source counts and
global sum agree. Whole-code coverage remains open.

The first mixed-arity extension also makes progress. An actual-fiber
kernel handles all arity-nine sets, including the 240 sets whose
zero-outside fiber is larger than their quiet kernel. Its initial bound
fails through long-lived boundary-state paths. One stronger marker then
helps arity nine but hurts arity eight. An interior common marker
z=9927/10000 clears the complete mixed a8/a9 class at r18,F<=1 in FLOAT:
log2 bound -202.5934919020, independently reproduced. This is wider
coverage for the same measured encoder, not yet a directed mixed bound.

The next arity-seven test exposed an artificial zero-column floor.
Enforcing the RS requirement of at least nine nonzero columns per active
group changes its FLOAT bound from +124.12 to -71.27 bits, without native
changes. An independent high-precision recombination confirms the gain.

A separate proof improvement computes the actual state-weight law of an
empty inner step. It removes almost all of the repeatedly charged refresh
density penalty while retaining the same uniform-live reference. The
empty-region improvement is at least 11.564160923761 bits, checked with
Arb256/384. Rebuilding the complete mixed a7/a8/a9 selected gate gives
-76.39182947 bits in FLOAT, including all regional boundary paths.
The arity-nine component improves from -347.50 to -1238.99 bits; arity
seven now limits this selected union. This is not whole-code coverage.

The creative construction branch also produced a
[rank-eight coupled-quadratic mixer](HD_BIT_SYNDROME_STAR_OUTER.md).
Ten additional XORs, with no extra field products, reduce its quiet
dimension from eleven to eight. Exact local and pair checks pass; the
coupling removes the new eight-hole pair cases of a simpler candidate.
It still has neither a measured native runtime nor a distance certificate.

Next: use the activity and exact-empty-step repairs to test a broader,
simpler all-arity proof for the unchanged 3.41-ms encoder. Do not invest
in more one-class arithmetic before testing that coverage strategy.
Keep the coupled-quadratic mixer available if an actual structural
cancellation obstruction remains. The active goal is still unfulfilled.

## A broader proof decomposition and a zero-product construction change

The next pass tested the proposed simplifications before investing in
their complete geometry. The result is not a whole-code certificate.
It identifies a stronger route than continuing one-arity envelope fits.

The exact scalar helper now obtains every zero-outside fiber enumerator
from the outside image of its complement, using integer MacWilliams
inversion. This avoids a separate exponential inside-support calculation.
However, collapsing the inner to one scalar potential loses too much.
At the saved r18/F<=1 mark, put c=log2 R11(18)=-84.62262821330349 and
b=log2 t<=0. The 140 quiet affine-plane row sets alone force that proof
majorant to be at least

\[
 7932.104519645711+256\max(0,c-18b)
 +18\log_2\!\left(140\binom{16}{9}255^{-23}\right)+648b.
\]

Its minimum is approximately +1947.03 bits, at b=c/18. This rejects the
scalar comparison at this mark, not the construction. No full scalar
geometry or optimizer run was needed.

A separate conditional test keeps both boundary states but assumes a
uniform information-set cover for every column subcode. Under that
unproved hypothesis, the cached bounds for arities seven through eleven
are negative. Arity seven reaches -111.70 bits. The same relaxation fails
for sparse arities and for arities twelve through sixteen. Limited exact
basis checks found no counterexample, but do not establish the hypothesis
universally. This is a possible component of a hybrid proof, not a
replacement certificate. The isolated control is
`theory/hd_uniform_cover_conditional.py`.

Two more structural ideas are now supported by exact calculations:

1. [Separate rare outer setups](HD_BIT_SYNDROME_QUIET_EXPURGATION.md).
   For the unchanged measured H* design, the probability that any of the
   8192 groups has an all-quiet nonzero word on at most seven active RS
   rows is at most 2^-50.771255305346. This global bound uses no output
   mark. It permits a failure decomposition into that bad-setup event
   and an unconditional first moment restricted to words with a
   nonquiet column in each low-arity group. Do not reuse independence
   after conditioning, or multiply this global charge by a source count.
2. [Vary the input order by column](HD_BIT_SYNDROME_COLUMN_PERMUTATIONS.md).
   Use P_j(i)=inv_GF16(i XOR j), with inv(0)=0. The schedule is fixed;
   only the existing diagonals remain random. Every four-row group then
   has at least five nonquiet columns, and every five-row group has at
   least two. The four-row statement also has a short quartic-derivative
   proof. This changes no field-product count. Whether the loads can
   absorb it without a latency penalty remains unmeasured.

The scheduled cached a7/a8/a9 comparison retains every original labelled
row set and multiplies its sixteen column matrices in order. Its mixed
FLOAT bound is -76.39183098, essentially unchanged; the a9 component
improves from -1238.99 to -1920.64. This test does not yet exploit the
new sparse nonquiet-column floors. The exact singleton kernel also gives
a cheap independent control: preserving simultaneous activity in all
sixteen regions yields -4003.46 bits for r18/a1/F<=1 in FLOAT.

The scheduled quiet-space calculation is stronger than the initial
four-row argument. Across all 12,870 eight-row sets, each can become an
affine hyperplane in at most one of the sixteen columns: 12,390 sets
never do, and 480 do once. Their total quiet dimension across columns
is therefore at most49, rather than the static worst case64. The exact
dimension helper covers all49,946 labelled sets of arities four through
nine; it retains each original row set across columns.

The corresponding exact global setup-event unions are now computed.
For the inverse schedule, all-quiet words on at most eight active RS
rows occur with probability at most 2^-105.587491005579 across the entire
8192-group code. Extending the event through arity nine gives at most
2^-45.106001459030. Arb256 and384 ceilings agree. These exclude only
all-quiet outer words, not low-distance codewords in general. The receipt
is `theory/hd_bit_syndrome_scheduled_expurgation_le9_exact.json`.

A first positive nonquiet-flag control is also implemented. Expanding
at the first nonquiet column before the group power gives an ordered
quiet-prefix/nonquiet/full-suffix sum. Including its68-bit convexity
cost, the scheduled homogeneous-a8 restricted first moment is -243.05
bits at r18/F<=1 in FLOAT. The unscheduled restricted control is -213.76,
but no small static-a8 bad-setup charge has been established. These are
restricted counts, not unconditional distance certificates. The helper
is `theory/hd_nonquiet_word_control.py`; it never subtracts upper matrices.

An [affine-partition Toeplitz proposal](HD_AFFINE_BLOCK_TOEPLITZ.md) is
retained as a different construction option. Four small randomized
blocks replace H* and both diagonals, with a whole-group pointwise
counting bound. Its estimated arithmetic is comparable, but its latency
and distance are untested. It is not the immediate implementation choice.

Next: integrate the nonquiet-column requirement into the two-state
sparse-family comparison. Test it first using the unchanged encoder's
global setup-event bound, then exploit the stronger deterministic floors
of the inverse schedule. Neither native implementation nor another
large regional scan is justified before that bounded proof test. The
certified fallback and its manifest remain unchanged.

Validation for this pass: the168-test existing-plus-new suite passes,
as do six additional first-nonquiet/dimension tests and five scheduled
setup-event tests. Independent review checked the column schedule,
incidence argument, and positive first-nonquiet matrix recurrence.
No benchmark ran, and no production or paper file changed.

## Sparse audit and distinct construction alternatives

The [next audit](HD_SPARSE_ACTIVITY.md) completed the bounded arity2..7
test at r18/F<=1. Retaining nine active columns clears arities two, three,
four, and seven in FLOAT. Arity five remains at+423.63 and arity six
at+105.53. Adding the stronger nonquiet setup exclusions and input-weight
markers barely helps; further tuning of those markers is parked.

The failure has a concrete explanation worth testing before changing the
encoder. On affine-basis five-row inputs, an outside-weight-four term
replaces five actual support patterns by330 and discards the five inside
outputs that are already determined. Exact projective rays now recover
those outputs. Both affine-rank-four and rank-five contributions need
repair; neither alone explains the complete remaining class.

The creative search also produced an independently shuffled local output
ensemble, a stronger nonbinary syndrome mixer, and an inner that emits
its refreshed state directly. Their algebra, implementation costs, and
specific risks are indexed in the audit. None is silently substituted
for the measured design or the certified fallback. The active goal still
requires complete directed coverage and a matched measured implementation.

## Extra stages inside the superblock

The [crossed-pair implementation and cost test](HD_OUTER_PAIR_NATIVE.md)
now compares one and two additional local stages. Their exact support
laws compose without new intermediate diagonals. Varying the row
differences removes the common four-row invariant planes of the simple
constant-wiring proposal. This is a useful local structural improvement,
not a new full distance certificate.

Matched K20 full encoding measures3.37ms for H*,3.85ms for one added
stage, and4.35ms for two. All variants pass the independent scalar
reference checks. Therefore additive stages alone miss the current
performance target. The next candidate should replace some existing
outer mixing with one crossed stage, with its sparse quiet class tested
first. The certified fallback remains unchanged and the goal remains open.

The [replacement screen](HD_PARITY_OUTER_PAIR.md#state-aware-low-arity-controls)
is now implemented with exact positive input-weight tails. An initial
finite-marker bound was dominated by the excluded zero word and is
superseded, not used to reject the kernel. With that loss removed, H0
plus one pair layer clears the selected singleton class but fails the
two-row bound. Hbin improves the two-row bound to+104.39 but fails the
singleton bound at+926.38. These FLOAT tests cover only r18/Fmax<=1;
they are not whole-code certificates or counterexamples to distance.
The one-product symmetric pair gives essentially the same selected
bounds. Thirty-nine local tests and independent algebra/controller
review pass. No new native experiment is justified by these results yet.

## Cheap mixing and varied column wiring

The [next bounded pass](HD_SCHEDULED_CHEAP_OUTER.md) changes the mixer and
the local input/output row order, rather than stacking more stages.
H_lambda uses 57 XORs and one fixed GF(256) product per column. Choosing
lambda in the GF(4) subfield gives exact small relation-kernel counts.
The prescribed inverse input schedule plus output-XOR schedule clears
all selected one- through four-row classes without a crossed layer.
The pair bound is -165.56 bits, the triple bound is -1934.58 bits, and
the complete four-row bound is -574.84 bits across all 1,820 sets.
They refer only to r18/Fmax<=1 at the pinned mark, not total failure.

This pass also removes an avoidable proof relaxation: independent
columns are bounded separately before chronological matrix composition.
All original row sets and their input-weight floors remain in the bound.
Independent review checked the new field embedding, support counts,
state order, route interface, and normalization. The certified fallback
manifest is unchanged. Middle/high arities, broader source coverage,
directed verification, and a matched native timing remain
required before closing the active goal.

The subsequent dense diagnostic gives a useful stopping condition for
this particular bound. One positive live-state path forces the current
common-measure comparison above +312 bits already at arity12, even
after all original RS support floors are restored. This does not refute
the construction. It rules out extending the same sparse measure to
those dense classes, so the new dual-code controller was not launched.
Use a different regional measure or the unique-preimage comparison for
the dense branch. The sparse results and 66 passing local tests remain
recorded in the linked construction note.

The prescribed H_lambda construction now has a checked isolated native
implementation. Matched serial K20 medians are3.442226ms for the candidate
and3.354983ms for the old H* control. Independent scalar, basis, gauge,
aliasing, and natural-output checks pass. Thus implementation cost is
close to the target, but is not a new certified performance result.

Fresh directed arithmetic supplies a63.942051987187-bit bound for the
class r18 with total active RS arity at least247, over all pairing profiles.
The low-total-arity complement and other group counts remain open.
Retaining quiet input activity now closes all4,368 selected five-row
sets in FLOAT at -168.60bits; the same-measure mixture of arities1..5
gives -149.04bits. Partial-coordinate accounting extends to larger
quiet fibers and passes representative arities6..11. The arity12
representative remains at+180.60, consistent with the dense-measure
obstruction. Three fresh directed uniform-comparison controls also fail
at total arity216; the best tested witness still starts at247 for60 bits.
The next test needs a support-resolved regional interface before
fractional averaging; the old dual formula already retains codimension.
Test its dominant live-state path before enumerating all remaining row
sets. The construction note indexes exact
scopes, pinned receipts, and native logs. The original certified fallback
and production library remain unchanged.

The next round closes the former arity12 gap through an exact Renyi
comparison, with no implementation change. All1,820 mapped complements
give maxI3 approximately1.0078591195655953. Holder before fractionalization
lets the unchanged uniform provider certify69.240061061278bits for r18,
exactly twelve active rows in every group, and all99 pairing profiles.
Fresh256/384-bit replays agree, including a512-bit rounded profile union.
The construction note and HD_LAMBDA_RENYI.md record the proof and hashes.
This is genuine new class coverage, not whole-code closure; mixed row
arities and the remaining group counts still require certificates.

The same round extends to a genuinely mixed dense-group class. Fourth
moments and actual scheduled column products certify total arities216..245
when every group has arity12..16. Union with the prior total-arity>=246
comparison gives40.82172502955bits for r18 across all99 profiles. Exact
orthogonal-symmetry enumeration retains all2,517 row sets and agrees with
the independent full1,820-set third-moment sweep. The remaining r18 gap
is total arity<=245 with a group of arity<=11. Bothprecision replays agree,
60tests pass, and no native, production, paper, commit, or push action
was taken. Extend the middle-arity moment catalog and its mixed-class
partition next, before enlarging the active-group range.

The next round adds all4,368 eleven-row local moment cases and retains
the RS-implied nine-active-column floor. Directed256/384-bit replays
close r18/all99profiles with all group arities11..16 and total198..225
at55.461496706851bits. Exact original-input-weight moment polynomials
also close homogeneous arity one at975.628296489068bits. The recorded
union with the preceding dense cover is40.821668515604bits, still partial.
Seven intermediate moment orders improve lower totals but not cutoff225.
The full marked arities1..5 uniform comparison fails, dominated by arity
two; this does not invalidate the earlier selected-profile sparse bounds.

The fixed activity-only comparison fails even when every nonzero column
is ideally mixed: its all-arity r18 expression is at least2^380.0995.
This identifies a proof relaxation that another column-local stage
cannot fix at the current mark and tested moment orders. Preserve the
construction and test a reference with an explicit whole-column zero
branch. In parallel, an isolated four-coordinate transvection can target
actual nonzero-moment spikes at a small XOR cost; it has no transferred
certificate or native timing yet. HD_LAMBDA_RENYI.md records all scopes.

Those two controls are now complete. The one-transvection variant keeps
the minimum two-source output weight six, but changes the troublesome
high-arity subspaces rather than removing them. Its directed log2 bound
for the tested total198..245 class is+166.2600, versus+163.2424 before
the change. The reference that explicitly accounts for whole zero columns
also gives no improvement on its tested grid. Both256/384-bit replays
agree on every completed bound. These are negative proof diagnostics,
not distance counterexamples. The next structural screen considers
several independent cheap transvections and, separately, local wiring
that couples different columns. No native code or fallback was changed.

The bounded structural search is complete. Of5,040 eligible two-stage
configurations,1,290 preserve the local two-input minimum six. The cheap
supports0123 and0167 also avoid the single-stage variant's enlarged
proportional-column class; shared circuit sums reduce their added cost
to seven XORs per column. This is a promising candidate, not a certificate.
A staggered adjacent-column XOR shear is a separate option: exact
uniform-comparison zero ranks double for every singleton and pair support
with even/odd-row wiring. It breaks the independent-column proof and
requires restoring eight D1 products per group relative to the measured
gauge-normalized kernel. HD_LAMBDA_RENYI.md records the complete search
scope, cost accounting, and next joint-distribution test.

The next bounded proof round extends the unchanged encoder's directed
r18 cover to the class with all group arities10..16 and total180..202.
That class has56.155410727944 bits; its union with the earlier partial
cover has40.821633581232 bits. Fresh256/384-bit ceilings agree and an
independent audit clears the normalization and final verifier.

Both modified mixers now pass fresh complete two-row FLOAT tests at the
selected r18/Fmax<=1 measure: -235.73 for the cheap double transvection,
and -395.23 for the joint-column shear. Seven of eight prescribed middle
supports pass for the double transvection. The failing six-row set has
an exact quiet-weight calculation, but its aggregate RS floor hides that
two original rows cannot use some quiet columns. The next gate enforces
one of those rows' nine-column floors separately. These tests neither
transfer the unchanged encoder's certificate nor supply a new native
timing. The original certified fallback and production code are unchanged.

The targeted row floor gives no useful gain. A recovered term of the
same proof expression has log2 value+380.17 despite satisfying all six
surrogate row floors. The loss is therefore in the nonquiet surrogate,
not merely aggregate activity accounting. Exact GF4 relation-kernel
enumeration removes that relaxation: the fixed six-row diagnostic changes
from+529.08 to-3010.59 with its true input-weight/output-support counts.
The complete calculation takes3.44 seconds. Independent review clears
the new counts, normalization, and coefficientwise chronology. This is
progress on the cheap double-transvection candidate, not a new whole-code
certificate. Reuse the local categories to cover broader labelled sets
next. A separate ideal paired-column comparison also leaves room for the
shear's dense branch; an actual moment upper bound is still needed.

Both follow-ups made progress. The reusable a6 table compresses564,357
relation kernels to19,414 exact categories, preserving total input weight.
It reproduces the original histogram and bound, and the union over eight
prescribed row sets gives-2984.44 FLOAT bits. This is not all8,008 sets.

For the shear, small parity witnesses give an actual fourth-moment upper
bound. All1,820 twelve-row sets and14,560 paired-column witnesses complete
in9.77 seconds. Two fresh directed evaluations certify285.081339551997
bits for r18, exactly twelve active rows per group, and all99 pairing
profiles. Each evaluation independently rechecks the entire witness
inventory. This new certified class belongs to the shear only; it must
not be unioned with double-transvection or unchanged-H bounds. Mixed
arities and other group counts remain open. The75 focused tests pass,
and no native code, production code, paper, commit, or push changed.

The next shear round completes all560 three-row sets, giving a selected
r18/Fmax<=1 FLOAT bound of-2391.10124575756 in74.77 seconds. Its exact
80-form catalog also reproduces all prior one-/two-row results. The
twelve-row directed class now clears40 bits for r16..19, with all source
profiles at each r; the fixed bound fails at r20. A separate complete
arity11..16 inventory adds the all-eleven-row r18 class at45.697980665791
bits, but does not close the mixed-arity complement.

An optimistic abstract-cover calculation proves that better independence
witnesses alone cannot make the current p2/3/4 relaxation close any of
arities4..10. This is a limit of that proof method, not the code. The
next proposed bound instead averages the relative scalar already present
between each pair of full-D1 columns. Distinct cancellation ratios affect
disjoint subsets of overlap targets. Under the positive common measure,
merging these cancellation classes yields a rigorous two-case majorant.
It must still be tested numerically; no new encoder stage is needed to
state this bound. HD_LAMBDA_RENYI.md records the exact scopes and receipts.

The first shared-scalar test now passes: two prescribed four-row sets
give-1375.05 and-1099.30 FLOAT log2 bounds under the same selected source
measure, without enumerating eight-variable relation kernels. Its cost
is0.447 seconds including a three-row control. All50 new focused tests
pass; independent review clears the ratio lemma and controller. Next is
complete four-row coverage, then larger single-column histograms. This
does not yet close the whole code or change any measured encoder timing.
