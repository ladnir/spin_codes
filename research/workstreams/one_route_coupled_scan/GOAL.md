# Goal and progress

Started 2026-10-20. Research log, not a production change. The fixed-pair
exploration ended at the redesign checkpoint below. The random-pairing
exploration is now active; the whole-code proof remains open.

## Current goal: cheap inner at near-half distance

### 2026-11-05 — Preserve the paused global-rank checkpoint

The user paused the goal after the November 1 work. This publication
checkpoint preserves sources and compact proof notes without resuming
numerical work. The fixed construction still lacks a full 45%/40-bit
certificate.

The coarse sloped-quota comparison failed. Exact first-layer checks
distinguish the sharp gaps created by restricting destinations from the
smoother, still nonconcave profiles of the full operator; see
[the scaling note](coded_state/FULL_POPULATION_SCALING.md). No general
population-compression theorem was obtained.

The [marker-specific gates](coded_state/QUARTET_GENERAL_MARKER_GATE.md)
identify two promising restricted comparisons: capacities `(1,1,2)` at
`y=3/2`, and `(2,1,4)` at `y=5/4`. These change only the proof weights.
Neither restricted result bounds all physical populations or certifies
code distance.

The general full-population Horner caller passed independent static
review and exact small checks: 4410 physical coefficients, 42 quota
rows, 42 endpoint rows, 84 cached-chain rows, 84 unit-endpoint rows,
42 scale comparisons, and four composed comparisons. Its original-network
production control was launched, but the turn was interrupted at the
pause and no completed result is retained here. Do not infer a pass or
failure from that launch.

On resumption, first obtain one complete bounded full-population result
at the smaller marker. Only a successful result would justify rigorous
outward replay and extension across the remaining row-count classes.
Local refinements remain parked; prior proof checkpoints are preserved.

### 2026-11-01 — A quota-aware comparison passes the restricted diagnostic

The [remaining-quota interface](coded_state/QUARTET_REMAINING_QUOTA.md)
now includes an exact interval-mass correction. Each physical transition
is a positive mixture of products of 128 local label-count polynomials.
An adjacent-mass Fourier argument bounds every coefficient atom by an
exact rational eta. An interval of d label counts therefore contributes
at most min(1,d*eta) of the evaluated transition. This corrects the earlier
extension from a short permitted label interval to all label counts.
The proof and its implementation pass independent review.

The exact helper passes 90 local-law checks, six central-coefficient
checks, 180 homogeneous products, 90 each for two mixed-product families
and mixtures, and 36270 interval checks. The optional caller preserves
its baseline and passes separate exact coefficient/band comparisons.
These checks establish the helper's scope; they are not a code certificate.

With singleton quota bands and this correction, the three-state control
passes for the first time. Four bounded witnesses take 41.848 seconds.
The final maximum regional growth log is 161.03135; the restricted
256-region comparison has logarithm -10463.8755, including its start,
row-set union, and restricted terminal. This remains a FLOAT upper bound
for paths confined to the three-state restriction, not for the full
physical operator. No codeword class or full distance certificate closes.

The next obstacle is explicit: extend the passing comparison to all
33153 physical populations without performing thousands of separate
Horner actions at every update. The literal singleton-quota recurrence
needs 229488 actions per witness. An independently audited sloped-band
and capped-mass interface is recorded in the remaining-quota note. Its
next bounded test uses the restricted model to check whether compression
preserves the passing bound; even width-256 bands with three slopes and
seven evaluations would still need 3024 full physical actions. Do not
authorize that full job without a further cost reduction. Preserve the
passing restricted checkpoint and exact eta helper. The construction
remains fixed; local refinements remain parked; the global goal is active.

### 2026-11-01 — Destination splitting is insufficient; retain the quota path

The [destination-split control](coded_state/QUARTET_DESTINATION_QUOTA.md)
passes independent source review and exact small checks. Its full
33153-state run completes in 197.291 seconds, with logarithmic regional
growth 642.6367 and tail upper bound +112868.759. The minimum physical
ratio for this witness is also 518.937 in floating arithmetic. Thus the
failure is not confined to a few rows that a small potential adjustment
would fix. The test closes no probability class. Further runs of this
fixed comparison are parked.

An [analytic return-family bound](coded_state/QUARTET_QUOTA_RETURN_FAMILY.md)
then stops a proposed exact-quota obstruction test before execution:
even all 277 layouts of its loops and three-step excursions contribute
less than 2^-21971 to the restricted full comparison. The discarded
helper was new and untested; existing proof and numerical sources remain
intact. This controls only that auxiliary path family.

An [exact first-moment diagnostic](coded_state/QUARTET_QUOTA_POPULATION_DRIFT.md)
shows that, when one update is forced to have a balanced quota, high
persistent overlap contracts strongly under the normalized moment
measure. All ten exact normalizations and an independent audit pass.
This explains why a large row mass alone does not establish a persistent
bad trajectory. It does not bound the actual regional operator, whose
quota is shared across sixteen updates. The next discriminating test is
a coefficientwise envelope that retains the remaining quota through
those updates, first on the restricted three-state model. Construction,
setup, and free terminal convention remain fixed. Full coverage is open.

### 2026-11-01 — Isolate the quota-comparison obstruction

The exact arbitrary-vector kernel now passes independent small-instance
checks and all 33153 production-row comparisons with the prior exponential
action. Its first full regional control gives a FLOAT logarithmic tail
bound of +137000.328. No failure class closes. An optional cache of the
positive powers T^k also passes exact small checks and cached/uncached
comparisons; one full action agrees within 1.12e-15 in log2 arithmetic.

More significantly, an
[exact three-state certificate](coded_state/QUARTET_REGIONAL_CAUCHY_OBSTRUCTION.md)
proves that the unsplit rowwise Cauchy comparison cannot close, even
with every positive quota evaluation and every positive potential.
Its regional growth exceeds 2^400, forcing the resulting full comparison
above 2^49660. The parent replay and independent source audit pass.
This is a limitation of that comparison, not an encoder counterexample
or an obstruction to exact quota conditioning. Further tuning of its
potential or quota grid is parked.

The next bounded test separates destination populations before applying
the quota bound. A three-state FLOAT diagnostic distinguishes the useful
split: separating persistent overlap gives growth log 177.1042, whereas
separating empty from nonempty hole populations still gives 258.696.
Only the former avoids the known obstruction. The full test uses the
fixed partition w<256 versus w>=256, exact physical-row actions, and
separate quota evaluations for the two destinations. Its bulk and free
endpoint comparisons pass independent mathematical review. The test
still has to control all physical populations; the three-state result
is not an upper bound on them. The full certificate remains open.

### 2026-11-01 — Exact arbitrary-vector population transfer

The main new result is an
[exact positive Horner recurrence](coded_state/QUARTET_EXACT_POPULATION_HORNER.md)
for the adaptive population operator on an arbitrary vector. For the
hole and selected quartics T and S, maintain P_k=T^k and
H_k=S H_(k-1)+binom(m,k) D_k(P_k), where D_k weights a hole monomial
by f(4k, its old-A degree). Induction gives the exact joint population
action after the one final old-multinomial division. Two independent
algebra reviews pass. Each stage uses two sparse 15-term products;
there is no dense physical matrix or general convolution. The powers
of T can also be reused between successive actions at the same x.

This removes both population-bin maxima and the minimum-envelope
projection, instead of further tuning either relaxation. The regional
quota still comes after sixteen fixed-x actions; that coefficient
comparison remains an upper bound. The next discriminating test is the
new kernel against independent small physical enumeration and the old
exponential action, followed by a bounded full regional control.

The partial-capacity envelope control is now finished: its best FLOAT
bound is +258810 bits, with about 80 bits of projection loss per update.
It closes no class. Broad capacity and slope searches in that comparison
are parked. The
[global algebra audit](coded_state/GLOBAL_RANK_ANNEALED_LIMIT.md)
also records exact limits of simpler rank reformulations, and a valid
full-availability normal form for cut witnesses. These statements pass
independent review. They do not supply the missing geometry tail.
The full certificate remains open and the construction is unchanged.

### 2026-11-01 — Row-mass control and a comparison without population bins

The nonlinear row-budget helper passes independent mathematical and
source review. Serial exact checks cover 135 capped-simplex vertices,
42 physical rows and composition comparisons, 21 regional quotas and
endpoints, two Collatz bounds, floating parity, and the workspace guard.
The full-size control preserves the earlier nine modes and three quota
evaluations. It changes the logarithmic tail bound from +275797.785 to
+275790.249, only 7.54 bits. Its best regional improvement is about
0.0264 bits. This is a noncertifying FLOAT comparison, not a failure
probability estimate or an encoder benchmark. No probability class closes.

Preserving total row mass is therefore not the main missing ingredient
in this comparison. Park further row-budget tuning. The next bounded
test removes population bins altogether: represent a positive function
as a minimum of exponential potentials, apply the exact physical-row
actions, and take the least upper envelope in that same finite family.
Compose sixteen steps at a fixed quota evaluation before minimizing
across evaluations. The interface passes independent mathematical review;
its implementation and bounded control now pass independent review and
exact small checks. At the same nine modes and three evaluations, the
two-witness envelope control gives +218581 bits; its single-mode fallback
gives +218576. It is still vacuous and does not improve the earlier
154-bin result. The fixed-slope projection loses as much as 82.74 bits
per update at the central evaluation. Thus removing bins has not removed
the need to represent the joint population transitions sharply. Broad
mode searches are parked; this result does not refute the exact operator.

The [partial temporal-capacity note](coded_state/QUARTET_PARTIAL_TEMPORAL_FLOW.md)
gives a different global proof tradeoff. Retaining the death penalty
allows an exact two-state support census over every death count. The
networks (2,1,4) and (4,1,8) would yield conditional field failure below
2^-557 and 2^-74 if their required simultaneous geometry tails closed.
Their half-row integer deficit targets are 184116 and 406616. Exact
rational and integer constants pass serial checks and independent review.
The smaller network is the next bounded control: its temporal fraction
is precisely enough to remove an identified one-target Hall overload.
The geometry tail, other row counts, and the full certificate remain
open. Neither the encoder nor its setup distribution changes.

### 2026-11-01 — Preserve transition mass and audit the proof-network tradeoff

The previous turn made concrete progress; the full certificate remains
open. The [population-bin helper](coded_state/quartet_adaptive_bins.py)
now implements the audited disjoint-lift comparison. Independent source
review and exact small checks pass, including regional coefficients,
endpoints, composition, and the workspace guard. All production controls
below are noncertifying FLOAT calculations, not encoder benchmarks.

The first 46-bin comparison gives log2 probability upper bound +275798
at y=2 and deficit threshold 60518. Refining to 154 bins, with the same
lift and auxiliary actions, improves it to +205607. The observed source-
bin variation falls roughly by half; destination-bin overcount increases.
A separate endpoint-matched lift gives +355171. Changing that lift also
changes the absolute auxiliary-action family, so it does not isolate a
pure gauge effect. None of these controls closes a probability class.

Four exact-algebra diagonal-path diagnostics retain actual adaptive
probabilities and regional quotas. At k=4,8,10,12 their conditioned
diagonal lower controls are approximately -178,-690,-1131,-1769 bits.
These lower bounds detect no obstruction in those particular subfamilies;
they are not upper bounds on the complete operator. The physical and
rational checks for [the path helper](coded_state/quartet_adaptive_paths.py)
pass independently.

The [general flow-cover interface](coded_state/QUARTET_GENERAL_FLOW_COVER.md)
now covers arbitrary eligible-row counts and integer proof-network
capacities (a,b,c), retaining the field penalty for deaths. This is a
change of proof comparison, not a change of encoder. It also exposes an
important failed shortcut: removing temporal transport entirely makes
the population calculation smaller and would permit an all-death field
bound below 2^-419, but the required geometry tail is false. Disjoint
Hall stars prove, for even one fixed half-row set, deficit-tail
probability greater than 1/256 at both relevant thresholds. The exact
quota and integer arithmetic are checked. The
[no-temporal note](coded_state/QUARTET_NO_TEMPORAL_FLOW.md) preserves both
the valid conditional field theorem and this obstruction. No such
obstruction has been proved for the encoder or the original nine-neighbor
comparison.

The next comparison retains each physical row's exact total transition
mass, rather than summing independent destination upper bounds. Its
monotone homogeneous composition and regional quota order are proved in
the population-comparison note. Keep the original linear bound as a
fallback. A partial temporal capacity remains a possible later tradeoff,
but has no certified tail and is not a replacement for the current goal.

### 2026-11-01 — Exact matching, an amplitude bridge, and adaptive global counting

The construction is unchanged and the full 45%/40-bit certificate remains
open. Local seven-type refinements remain parked. The half-row flow is
now proved **equal** to a matching with at most nine input neighbors per
source; it is not merely bounded by a restricted matching. A minimal-
transport integral-flow argument removes paths longer than one temporal
edge. Independent proof review and exhaustive small checks pass. On the
same 64-update diagnostic, both formulations give integer deficit 467.

The [flow-to-field theorem](coded_state/QUARTET_FLOW_FIELD_BRIDGE.md)
also resolves an important scope issue. Excluding the sufficient
dimension obstruction does not itself prove distance. A stronger,
simultaneous bound on the flow deficit does bound **all output amplitudes**
in a low-death class, including cancellations and the free endpoint.
For h=4096, integer deficit caps 60518, 45914, and 33518 would give field
failure below 2^-86, 2^-86, and 2^-60 for death counts at most 65536,
98304, and 131072, respectively. These constants pass exact integer
checks and independent review. The required flow tails are not proved;
other row counts and larger death counts still require coverage.

The [dual Hall count](coded_state/QUARTET_SHORT_MATCHING_HALL_COUNT.md)
normalizes away optional negative witnesses. Its smaller 129-state
comparison drops own-input coverage. Exact regional coefficient actions
are now implemented with positive polynomial convolutions, shared
scales, and no FFT. A bounded full-size action took 5.32 seconds and
improved the tested comparison, but remained positive. More decisively,
an exact rational constant-population witness proves that this pair-only
comparison at y=2 exceeds 2^42132 even with exact quota extraction.
This rules out further potential tuning of that comparison at this
marker. It is not a counterexample to the code or a lower bound on its
failure probability.

The replacement is the [adaptive Hall normalization](coded_state/QUARTET_ADAPTIVE_HALL_NORMALIZATION.md).
It processes a global maximizing witness backward, retaining the future
coverage when bounding each branch difference. A decision-tree argument
gives a positive moment bound that keeps own-input coverage and exact
regional quotas. Both independent mathematical reviews pass. Exact checks
cover 259201 physical input subsets, 1007 canonical maxima, 2014
decision-tree majorants, and 450 local and 450 endpoint coefficients.
Its exact
population transfer retains two nested hole sets; a coefficientwise
scalar maximum loses their transition constraints and already has a
separate exact lower bound of 2^135200 at y=2, valid for every positive
pair of scalar-potential parameters. Do not launch another scalar grid.

The [disjoint-population comparison](coded_state/QUARTET_ADAPTIVE_POPULATION_COMPARISON.md)
passes independent mathematical review. It uses exact positive
exponential-potential actions and a common disjoint lift to smaller
population bins. Destination tilts are minimized per physical row before
the source-bin maximum; quota minima occur only after a full regional
composition. This specifies a computable comparison without a dense
33153-state transition matrix. It does not yet supply a numerical tail.

The basic [population action](coded_state/quartet_adaptive_population.py)
is implemented and independently audited. Its positive quartic expansion
agrees with exact rational and physical-partition checks at widths four
and eight. A single width-512 action retained all 33153 physical rows and
matched the unit-potential closed form with maximum log2 error below
3.3e-15. It took 1.59 seconds in a bounded serial proof calculation;
this is not an encoder timing. The bin composition and global tail
calculation are not yet implemented. Production arithmetic in this
action remains FLOAT, so any successful tail will need directed
verification before becoming a certificate.

Next evaluate a positive comparison that preserves these population
transitions, then target the stronger deficit thresholds from the
amplitude bridge. Success requires actual flow-tail bounds and explicit
coverage of the remaining row/death classes, not another conditional
interface or an obstruction-only test. No encoder, setup, paper, or
benchmark change was made.

### 2026-11-01 — Replace the support union by a fixed-row flow deficit

The [fixed-row flow identity](coded_state/QUARTET_FIXED_ROW_FLOW.md) has
passed mathematical and source review. For a fixed eligible-row set H
with p=|H|/8192, its exact maximum-flow deficit is
max_S ((1-p)L(S)-pD(S)-|V(S) outside H|). Every dimension obstruction
requires deficit at least p(K-W)+1. This controls every support mask
simultaneously, with no additional union over ranks, volumes, births,
deaths, or footprints. The remaining union is over H and its size.

Independent whole-region windows supply feasible subflows. Their
deficits add as an upper bound for the global deficit, with boundary
losses already included. At h=4096 the full threshold is 52430; average
window budgets are 409.609375 for 32 updates and 819.21875 for 64 updates.
Two bounded, distinct seed-zero instances give exact deficits 229 and
233.5 respectively. Integer max-flow equals independently reconstructed
cut capacity in each case. These observations justify investigating this
global criterion; they do not establish a probability, typical value,
or a certified maximum deficit.

The standalone [flow checker](coded_state/quartet_fixed_row_flow.py)
passes 18432 exact tiny support comparisons, 72 exhaustive minimum cuts,
six terminal checks, five quota checks, and 76 greedy-flow comparisons.
Both diagnostics use exact regional H quotas and the actual random
quartet partitions. No encoder benchmark, setup change, or field-map
surrogate is involved.

A backward greedy feasible flow is too weak on these instances, giving
deficit upper bounds 828.5 and 1510.5. Moreover, the 15-case pointwise
exponential-potential analysis of that greedy process has an analytic
stationary-mixture obstruction; no parameter search was run on it.
The full flow and a restricted two-step matching certificate retain
additional transport. Next prove a deficit tail using that transport
and the actual partition distribution, then pay the single H union.
The full 45%/40-bit certificate remains open; the goal stays active.

### 2026-10-31 — Close a global phase family and retain whole-run geometry

The goal remains active: the full unchanged K20, rate-half, 45%/40-bit
certificate is open. Local seven-type refinements remain parked.

The new [routing-rank theorem](coded_state/QUARTET_SYNCHRONOUS_OBSTRUCTION.md)
excludes every support whose state is wholly empty or wholly full at each
update from satisfying L>R and R+D+1<=W, except with probability below
2^-119. It counts all 4096-update phase histories and uses negative
association of the regional routing samples. The number of partially
covered regions is controlled by the number of closed runs. This keeps
the relationship that a scalar one-update maximum loses.

A deterministic Lipschitz argument extends the same count to all masks
within 16r coordinate-time changes of a synchronous base with r closed
runs, with total exception below 2^-73. No perturbation-set union is
needed. The r=0 case has radius zero and is impossible deterministically.
Both statements are uniform over the quartet partitions and field maps.
They exclude these sufficient dimension obstructions; they do not exclude
every bad codeword or all partial-state obstructions.

The derivation passes independent review. Parent-run exact checks pass:
26 constant comparisons, 75 small-range inequalities, 4432 synchronous
histories, 12960 partial-region comparisons, 19 exact history-count
identities, 81 suffix rank checks, 57384 one-bit Lipschitz comparisons,
4218 obstruction-threshold comparisons, and dependent-route/moment toy
checks covering 1929 route outcomes. The source is
[quartet_synchronous_obstruction_checks.py](coded_state/quartet_synchronous_obstruction_checks.py).

Separately, [run saturation](coded_state/QUARTET_RUN_SATURATION.md) fills
every occupied allowed run and every open run. It preserves the exact
selectors, births, footprint, and outer rank; volume increases and deaths
decrease. Thus dimension obstructions have a whole-run normal form.
Ineligible-birth pruning is valid for a union over profiles, but can change
the profile. Exact checks cover 8257 independent saturation histories.
The [saturated-cover count](coded_state/QUARTET_SATURATED_COVER_COUNT.md)
retains every actual live-hit condition. Its pointwise diagnostic is
+247479; the averaged-cone diagnostic is +299226, including the required
55204-bit conditioning cost. Tight FLOAT finite-box comparisons do not
make these certificates or construction obstructions.

The new [513-population regional control](coded_state/QUARTET_ACTUAL_REGIONAL_RANK.md)
retains actual state populations, hits, births, and deaths through each
16-update region. Its exact tests and static audit pass. At the same four
markers, adding the regional footprint-capacity penalty improves the
hard-class FLOAT bound from +843844 to +768162. Both are vacuous and this
control reintroduces within-run support multiplicity. No larger grid or
claim of optimizer optimality follows from these two tests.

The next checkpoint above pursues a genuinely global replacement for
the remaining support union. No deficit-tail theorem or full-code
certificate is claimed. The synchronous proof is a usable checkpoint,
not permission to subtract its histories from an unmodified positive
moment.

### 2026-10-31 — Separate mandatory deaths from optional run capacity

The global-rank goal remains active and the construction is fixed. The
full certificate is still open. A fresh review identifies a stronger
global support filter, rather than another local trajectory refinement.

The [forced-run packing theorem](coded_state/QUARTET_FORCED_RUN_PACKING.md)
characterizes exactly which state volumes fit inside an allowed selector
history with at most a given number of deaths. Every closed run started
by a forced birth spends a death first. Only the remaining deaths can
activate other closed runs. The resulting capped inequalities strictly
strengthen the earlier unlabelled run bounds, including on small examples
that pass all of those earlier caps. The proof includes open terminal
runs and an exact support/amplitude polynomial inside the global quotient.
Two independent mathematical audits pass. Exact checks cover 2376 support
families, 19720 volume intervals, 183354 capped inequalities, 167 fixed-birth
amplitude sums, and 968 independent multi-run families with 7808 further
volume intervals. These establish the interface, not a production margin.

The full regional footprint capacity is also now implemented before
16-step composition. At fixed eligible-row count h, it gives
sum_regions(v_region-h)_+ <= R-128h, with no threshold union. Its exact
small checks and source audit pass. Four bounded marker controls remain
positive. The best is +148157, a small improvement from using the exact
free endpoint; the regional penalty itself does not improve that best
control. At one other marker choice it reduces +199198 to +162661,
showing that it removes some histories but does not close the class.

The new scalar implementation retains both old and new death caps, and
can also promote compulsory singleton hits without changing f or the
amplitude normalization. Its exact checks pass. A separate epigraph
barrier control substantially improves the optimizer diagnosis: at the
hard profile and D<=104857 its primal is +177152.55 and its FLOAT finite-box
dual lower estimate is +177140.02, a gap of 12.53 bits. This is evidence
against further tuning of that scalar relaxation in the tested box,
not a rational certificate or an unbounded-family impossibility result.
The barrier derivatives, feasibility guards, and dual comparisons pass
their checks; independent static review confirms the finite-box scope.
The earlier all-death +148157 control remains stronger and still fails.

Two exact small counterexamples isolate an additional loss: packing
state volume into allowed runs need not hit every selected quartet.
One even passes the aggregate live-volume condition. Sole-hit positions
are compulsory state positions, but they do not become extra prescribed
birth amplitudes. A mandatory-run count, fractional packing of run-hit
clauses, and a whole-run Lagrangian give valid global feasibility filters.
Their interfaces and counterexamples are preserved in the forced-run
note. Singleton promotion alone is not asserted to exclude every example.

Next retain actual live-hit constraints jointly with the regional count.
Further scalar parameter searches are parked. The scalar row-classification
filter is not a promising missing factor: conditional on covered H-slot
geometry at the hard profile,
almost all within-group row labellings already meet the degree thresholds.
Nor can a different elimination raise the state-equation rank above R(V).
The low-rank remainder still needs a successful geometric exclusion.

### 2026-10-31 — Group amplitudes globally; retain regional capacity and births

The global-rank goal remains active. The construction is unchanged, local
refinements are parked, and the full 45%/40-bit certificate remains open.

The [live-selector quotient](coded_state/QUARTET_LIVE_SELECTOR_QUOTIENT.md)
now retains actual state volume before applying a fractional power. An
exact scalar certificate bounds the profile
(h,R,L,f,F)=(1536,753664,730726,4,851968) by 2^-305581.
Its source and derivation pass independent review. The corresponding
geometry-only control is positive, so the retained field rank matters.
The [adaptive-zero budgets](coded_state/QUARTET_ADAPTIVE_ZERO_BUDGET.md)
give separate all-message exclusions and a sharper death upper bound.
The temporal selector also has a third rational control, giving -1060
at (2560,786432,779879,4,950272).

The [one-active-input-row theorem](coded_state/QUARTET_ONE_ACTIVE_ROW.md)
closes that entire message class below 2^-82, regardless of how many rows
its footprint makes eligible. A coarse state-growth level, the output
budget, and the number of input-occupied regions give a contradiction
outside three rare-event families. All thirteen integer growth bounds,
128 weight budgets, and independent structural checks pass. The proof
has been independently reviewed; it is not merely a selected tuple.

Two additional global restrictions repair the selector overcount. First,
each eligible row supplies only one position in each routing region, so
C_H-f <= sum_regions min(h,G_region). Second, a selected quartet disjoint
from the preceding selected set requires a distinct forced birth. There
are at most f such activations. Both restrictions pass independent
mathematical review. Their positive 129-state transfer passes physical
small-instance checks, including the free endpoint and regional blocks.
The combined bound now has an integer certificate: with rational markers
beta=35/16, z=459/1024, c=33/32, u=1/256, its directed exponents are -4453
at (2048,819200,819200,4,950272) and -18217 at
(2560,819200,829030,4,983040). Each includes all death counts and qualifies
for the established 96-bit profile union. Exact positive matrix products
use 512-bit upward rounding; the valid 60-bit version was too coarse.
The small-instance checks and two independent source audits pass. These
are additional whole-profile classes, not the full remainder.

The wider diagnostic grid still has 14730 open points after these new
screens (down from 20435 before the birth/region transfer and one-row
class). These counts are diagnostics, not a coverage fraction. Its worst
remaining examples now have many forced births; paying an independent
birth-mask binomial and then a birth-activation marker is too loose.
The next count integrates the birth choices inside each selector
transition. Its first bounded controls reduce several large upper bounds
substantially but remain positive, including about +148464 at
(4096,819200,838861,65536,1081344). Separately, the output-aware quotient
must retain closed allowed runs: a proved negative control rules out a
substantial gain from output markers after those runs are discarded.
The exact relation |D-C_A|<=A0+f-L now connects actual deaths to these
runs. A finite local supersolution retains that global statistic without
a large local-history dynamic program. Its exact toy checks pass, but two
bounded searches remain positive (+177477 and +191952); neither search
certifies an optimum or an obstruction to the whole scalar family.

The remaining proof has a necessary geometric component. If a support
mask satisfies L>R and R+D+1<=W, dimension and Singleton give a bad
codeword for every field-map choice. Such masks must be excluded by the
routing/partition geometry. This implication is not a claim that a mask
exists. For the current h=4096 profile, it includes deaths through 124517.
Larger-death slices can use the grouped output-weight field penalty.
The audited [death-cover interface](coded_state/QUARTET_GLOBAL_DEATH_COVER.md)
keeps both arguments in one finite, explicitly summed cover. The next
global bound now averages the actual random partition instead of
maximizing over each local quartet separately. Each selected quartet
must hit the preceding allowed state set, reducing the comparisons to
24897 reachable population pairs. This requires positive coefficient
row comparisons, not a transition matrix of that size.

The bounded averaged search did not improve +177477. A sum of two
positive potentials reduced that bound to +170881, still weaker than
the existing all-death +148464 control. Retaining regional H quotas and
closed runs together gave +177498. Exact small-instance checks and
independent derivation audits pass, but none of these positive controls
closes another class. A separate rational return-floor test is negative,
so it does not establish a limitation of the potential method either.
These searches are parked; they establish neither optimality nor a bad
codeword. The next structural review concerns regional concentration:
R(V) is at least the sum of the 128 largest regional footprint sizes.
That constraint must be retained jointly with the births needed to
start and stop dense stretches.

A joint regional-quota transfer was also derived and independently
audited. Three bounded controls remain positive and do not improve the
existing external-route controls; a fourth hits a zero-parameter proposal
boundary and gives no result. Further quota-marker tuning is parked until
the closed-run information is retained. An exact
[cut-interface obstruction](coded_state/QUARTET_DIMENSION_CUT_INTERFACE.md)
also explains why unconstrained min-cut scalarization cannot settle the
geometric remainder: empty/full masks defeat that relaxation. No flow
solver or large computation was added.

The five-state edge-potential relaxation is parked. An exact conditional
mean obstruction shows that no tuning of its continuous penalty
parameters can close its old h=1536 profile; the stronger quotient rank
already closes that profile. This is an obstruction to that relaxation,
not a bad codeword. Prior proofs and helpers are preserved. All jobs are
serial, root-owned, single-threaded and file-free on Peach. No construction,
library, paper, benchmark, commit, or push changes accompany this work.

### 2026-10-31 — Global selector capacity closes the former hard profile

The construction is fixed and local refinements remain parked. The
global-rank goal is active; the full 45%/40-bit certificate is still open.
The new [selector-capacity theorem](coded_state/QUARTET_SELECTOR_CAPACITY.md)
avoids a state-history union. An actual history must select enough
H-rich quartets to cover its eligible outer rows. A simultaneous tail
bound on the available quartet labels limits this selection and gives
lower bounds on the retained RS rank.

The entire former hard profile, including every death count, is now
excluded below 2^-4094. Exact directed arithmetic also verifies all
24573 (h,r) rank-plane thresholds, with total exception below 2^-77.
The bounded diagnostic grid loses 34422 of its 58778 previously open
points; 24356 remain open in that diagnostic. These grid counts do not
constitute a full certificate or show that the surviving tuples occur.

The remaining question is temporal consistency: selected H-rich
quartets must overlap enough at consecutive updates to carry the
state. The [temporal-selector theorem](coded_state/QUARTET_TEMPORAL_SELECTOR_CAPACITY.md)
now closes two further whole score classes. A scalar rational
supersolution replaces the exploratory matrix calculation. Directed
arithmetic gives exponents -11774 and -6002 at the two hardest
remaining sampled tuples, with all death counts included. Its fixed
markers exclude another 377 previously open grid points; 23979 remain
FLOAT-open after both selector screens. Separately, an independently audited
[quotient enumeration](coded_state/QUARTET_LIVE_SELECTOR_QUOTIENT.md)
sums state amplitudes over allowed runs after choosing live quartets
and forced births. It preserves the RS residual rank and free endpoint;
its full geometry sum is not yet bounded. Exact GF(4) checks verify
the run formulas, residual equations, and support containment.

Next test a five-state global edge-potential bound that prices both
the number of selected quartets and their temporal overlap. The
inequality and endpoint accounting pass independent mathematical
review; its bounded diagnostic implementation is not yet a certificate.
The small-h near-burst corner also needs a separate complete bound.

The earlier exact death slices, exact-birth controls, and regional-bin
controls remain preserved. The latter two do not improve the hard
profile enough and are parked. Continue with the global selector and
rank bounds, not more local history refinement. No production code,
construction, paper, performance measurements, commits, or pushes change.

### 2026-10-31 — Diagnose the global low-rank overcount before extending it

The goal remains active and the full 45%/40-bit certificate is open.
The construction is fixed; local refinements remain parked. The new
[actual-state/MDS joint count](coded_state/QUARTET_ACTUAL_MDS_JOINT.md)
retains zero outer inputs inside the support count, restores each
regional quota, and keeps nonzero zero-boundary excursions. Exact toy
checks and independent mathematical/source audits pass. Its tested
primary-profile upper bound is still positive, about +675044 bits.

The useful new result is a precise diagnosis of that failure. An
[explicit persistent family](coded_state/QUARTET_LOW_BIRTH_MOMENT_OBSTRUCTION.md)
has exactly four forced births, satisfies the relaxed local constraints,
and makes the present envelope exceed one even with exact quotas and
exact birth counting. Its footprint is nearly full: it is outside the
target remainder and its proposed eligible-row label is inconsistent
with its true row degrees. This is an obstruction to the relaxation,
not a bad codeword or a failure of the construction. Exact birth
truncation alone is therefore not the next large implementation task.

Two other proposed refinements are now ruled out quantitatively:
the combined deterministic output budget excludes no history in the
primary profile, and selected-zero MDS filters gain less than three
bits at the current weights. The
[audited derivation](coded_state/QUARTET_MDS_REFINEMENT_LIMIT.md)
also identifies the actual death threshold where the field penalty
changes. Exact coefficient extraction is not a substitute for retaining
that death count before choosing the output marker.

Bounded stronger-field and profile-matched controls remain positive.
Two initially hit arithmetic range guards and yielded no result; an
explicit positive-floor FLOAT majorant permits later diagnostics without
discarding tiny terms. Its exact-reference checks pass, but it is not a
directed certificate. Comparison-operator diagnostics move from dense
populations to short excursions as the weights change. They describe
the upper matrix, not an empirical distribution of code histories.

Next retain the already proved coverage and actual census complements
INSIDE the global count, using the positive multipliers recorded in the
joint-count note. Combine this with death-conditioned field weights,
test the same hard profile, then seek a covering family. Do not merely
subtract prior failure allowances from an unrestricted moment, restart
the local DP, or launch an unconstrained marker grid.

All numerical jobs are root-owned, serial, single-threaded and file-free
on Peach; every job has exited at this checkpoint. No construction,
encoder, paper, benchmark, commit, or push changes accompany it.

### 2026-10-31 — Actual-state classes and an exact geometry-only class

The global-rank goal remains active. The construction, once-sampled accepted
setup, free terminal, and 45%/40-bit target are unchanged. Local seven-type
refinements remain parked. This continuation proves additional classes;
it does not close the whole low-rank remainder or exhibit a bad codeword.

- Count actual nonzero state supports directly. Their field weight factors
  as lambda^Ls t^D times a free-terminal factor at most one, with
  lambda=1+qx and t=qx^2/lambda. This avoids the availability fixed-point
  relaxation and gives stronger closure factors. Six rational branches
  reuse the integer-verified 513-state moment. Whole-profile exponents
  are -34041 at the earlier test point and -134165 at a dense test point.
  Each qualifying profile class includes a global 96-bit union allowance
  and has total failure below 2^-76. Statistics are explicitly distinct
  from canonical availability statistics.
- Proved separate actual-state rank and weighted-closure census classes.
  Their complements refine the actual death interval; their probability
  budgets remain separate. Exact physical-support, rational-field,
  signed-rounding and complement checks pass, as do fresh moment and
  profile certificate runs. Independent audits confirm the interfaces.
- A second exact class uses the geometry-only polynomial, with no matrix
  witness or field penalty. Integer binomials and powers certify the
  (h,R,Ls,f,F)=(3584,524288,526950,4,954368) profile at exponent -72111.
  The winning marker is s=1, so the same bound is independent of Ls.
  Monotonicity in f extends it to every admissible f>=4 congruent to F
  modulo four at these fixed h,R,F. All death counts are included.
- The broad actual-state scout now removes the independently proved
  census classes and compares the moment with an optimized geometry-only
  control. It still finds open low-birth ranges with rank close to Ls.
  For example h=4096,R=Ls=786432,f=4,F=1048576 has a FLOAT upper bound
  around +977144. This is a failed upper bound, not a realizable-history
  assertion or a counterexample. Finite-grid counts are not coverage
  percentages or certificates.
- Implemented a joint forced-birth/run moment and the deterministic
  occupied-region threshold implied by the eligible-row quotas. Positive
  region powers pass independent exact-reference checks. Two bounded
  searches do not close the selected low-birth profiles. A separate
  two-update availability-separation count passes exact checks and gains
  about 208230 bits in its selected diagnostic, but still fails after
  routing; it is preserved as a secondary control.

Details: [actual-state classes](coded_state/QUARTET_ACTUAL_STATE_PROFILES.md),
[geometry class](coded_state/QUARTET_GEOMETRY_PROFILE_CLASS.md), and
[availability separation control](coded_state/QUARTET_TWO_UPDATE_SEPARATION.md).
Next keep the actual-state global count and retain the zero-input
constraints outside eligible outer rows jointly with state transitions
and regional quotas. A separate routing probability and a support count
still admit incompatible choices in the small-birth, near-equal-rank
range. Test that range before extending any marker search or finite cover.

All numerical jobs were root-owned, serial, single-threaded, guarded and
file-free on Peach. Every job has exited. No encoder, paper, construction,
benchmark, commit or push changes; the frozen pair helper hash is unchanged.
The full goal remains active and no genuine structural obstruction is known.

### 2026-10-31 — Exact run, quartet-occupancy and routing class certificate

Previous work was progress. An exact 513-state moment now sums all
availability closures jointly with random-quartet concentration. Integer
verification proves K v<=32v in every row and includes the free endpoint,
giving a 22323-bit moment bound. Combining this with forced-birth routing
and the grouped field bound proves a new whole-history profile class
below 2^-76, including its global 96-bit profile union.

The h=2048 profile (R,L,f,F)=(900000,900000,20000,1150000) has exact
exponent -34041; its entire half-width-512 box has uniform exponent -5285.
All D are already summed. Exact checks, a fresh certificate run and an
independent audit pass. The three older boxes centered at F=1200000 are
empty because F-R>128h throughout; their numerical bounds did not
identify a remaining obstruction. Passing the replacement box's aggregate
screen does not assert that it is realizable.

Details: [exact run/occupancy certificate](coded_state/QUARTET_RUN_OCCUPANCY_CERTIFICATE.md).
The construction is unchanged and the full 45%/40-bit goal remains open.
Next combine this class with the existing exclusions and screen the
remaining common (h,R,L,f,F) ranges before further moment calculations.
All numerical jobs were serial, guarded and file-free on Peach; no
production changes or new jobs accompany this documentation checkpoint.

### 2026-10-31 — Joint availability; accepted-map constraints in the global count

Previous work was progress. This continuation adds checked global
interfaces and necessary incidence constraints; full closure remains open.
The construction, accepted setup, free endpoint and target are unchanged.

- Implemented an exact joint availability/unavailable-coverage interface
  with safe retained subsets and positive count bins. All local, composed,
  zero-operator, terminal and quota checks pass. More retained information
  improves the selected bound, but does not close it: the broadest tested
  subset comparison gives +379257 FLOAT log2 units at h=2048. The older
  tagged control remains numerically better; these are different relaxations.
- Audited a 15-state tensor-norm comparison. Its exact checks pass, but
  the tested norm bound remains above +1.1 million. An invariant diagonal
  exposes ineligible all-full histories. No broad rescaling search is
  justified without carrying the missing row-class restrictions.
- Proved that an old-hit quartet with r available positions permits at
  most r-1 unavailable positions outside the eligible outer rows. Those
  rows have zero input, and the accepted A map is superregular. Also
  proved an output minimum max(0,j-r+1) when j quartet positions lie
  outside those rows, using the accepted A+delta I map. Forced births
  add their own disjoint output cost. Independent audits agree.
- The aggregate version of the first condition is too loose in a tested
  marker probe. Implemented the exact positive local filter instead,
  including an optional joint output-budget marker. Exact label-counter,
  finite-field, endpoint and monomial checks pass. Added a full-count
  Chernoff/bin comparison retaining every E count; its first coarse
  selected diagnostics remain positive: +716079 with a broader tilt
  family, +918902 with the default grid and exact axis slices. A tested
  piecewise potential is worse. No new whole-row-class certificate results.
- Proved a genuine obstruction criterion: arbitrary support masks U with
  L_U>R(V(U)) and L_U+D_U<=W force a low-weight codeword for every local
  map. No such masks have been found. The obvious unconstrained min-cut
  relaxation is analytically trivial (empty/full extremes), so no large
  graph search was launched or mistaken for a certificate.

Details and numerical status:
[joint global checkpoint](coded_state/QUARTET_JOINT_AVAILABILITY.md).
All numerical jobs are serial, guarded, single-threaded and file-free on
Peach. No encoder or paper edits, benchmarks, commits, pushes or raw dumps.
The frozen diagnostic hash is unchanged. Local refinements remain parked.
Every job from this round has exited. Next compare the two global
counting methods on the SAME hard (R,L,f,F) boxes, using their conditional
field bounds, before further unconditioned tuning. The output-budget
condition excludes no geometrically consistent history in the recorded
middle profile, so it is not the next closure plan.
The full goal remains active; no genuine obstruction has been found.

### 2026-10-31 — Availability recovers the footprint; weighted census proved

Previous work was progress. This continuation adds global class theorems
and a checked alternative transfer; it does not close the whole target.

- Proved V_t=Gamma_t(A_(t-1)) union A_t for every realizable footprint.
  Thus a fixed availability history determines the footprint once
  partitions are fixed. Both independent audits agree; 200704 exact
  support histories verify recovery and consistency on toy geometries.
- Derived exact and padded-run censuses. A positive 2x2 matrix gives
  sum_A 2^(-2jD)<=2^c_j, with exact rational constants. It proves two
  families of global rank/closure classes below2^-76 each, preserving
  every free terminal. Combining sixteen j values in both families
  costs at most2^-71. The full remainder stays open.
- Exact checks pass211 profile counts and20 weighted censuses. Indexed
  constants and marker substitutions in the availability-census note.
  Selected tagged-marker probes remain positive; do not confuse the
  class theorem with numerical closure of its complement.
- Implemented a direct513-state availability transfer, fixing candidate
  masks before averaging partitions. It charges closures exactly and
  has no tag approximation. Its530 exact matrix comparisons, terminal,
  partition, zero-history and regional checks pass. It nevertheless has
  a worse selected FLOAT bound (+748551) because dropped consistency
  permits artificial dense availability masks.
- Proved the canonical exclusion2F-L-f<=N-512. Its scalar marker improves
  that direct bound to+641425 at one tested weight, still much worse than
  the best tagged control (~+100811). More scalar tuning is not justified
  without retaining the missing joint consistency information.

Details: [availability census](coded_state/QUARTET_AVAILABILITY_CENSUS.md).
Next retain unavailable covered positions jointly with availability in
the global count, using their required disjointness from old-hit quartets.
An audited two-count positive polynomial now specifies this interface
exactly at the availability fixed point. Its 131841-state dense basis is
not implemented; the next step requires positive compression or a
matrix-free evaluation, not another independent scalar penalty.
Keep exact quota conditioning, full row/message coverage and the free
terminal. Local seven-type refinements stay parked; construction fixed.
All jobs were serial, bounded, single-threaded and file-free on Peach.
No encoder/paper changes, benchmarks, commits, pushes or raw dumps.

### 2026-10-31 — All-run closures give another global class certificate

Previous goal work was progress: exact run bounds, implementation, audits,
and completed serial diagnostics. This continuation preserves the entire
fixed-construction 45%/40-bit goal. The full low-rank remainder is open.

- Implemented joint (R,L,f,F) boxes with optimized conditional field
  bounds. Forty-five field-optimizer checks and an independent moment
  recomposition pass. One previously troublesome h=2048 profile has
  FLOAT log2 product below -641838; a nearby low-birth box still fails.
  These examples are not a class cover or a new certificate by themselves.
- Proved a joint all-run bound c^f lambda^(L-f)d^D_A with
  d=max(min(1,lambda/(q+1)),(1+qx^2)/lambda). Two independent audits agree.
  It includes every run and free endpoint and does not duplicate the
  earlier forced-run closure penalty.
- At x=2^-16 this proves the complete-history class
  32R-16L+15D_A>=17196732 below 2^-76, using only one 2^N footprint union.
  Exact rational constants pass. The remaining integer inequality is
  implemented as a canonical marker before relaxation and tag rounding.
- Derived a positive tagged transfer with a common quota-conditioned
  basis. Tags are auxiliary proof variables, not encoder randomness.
  Coarse thinning and sharper barycentric rounding preserve upper bounds.
  Checks pass 616 matrix entries, 48 paths, 68 terminals, 26 old-transfer
  identities, 1794 rounded probabilities, 2184 gate inequalities and
  384 run cases. The 129-state relaxation and all old checkpoints remain.
- Full-width FLOAT bounds improve but remain positive. A broad eight-level
  tag grid gives +100811 at the selected h=2048 point. A grid capped at16
  shifts to larger artificial histories and regresses to +116935. Stronger
  output penalties also worsen selected bounds. None proves full coverage.
- An aggregate row-classification filter can save fewer than193 bits for
  the old illustrative middle profile, so it cannot repair the roughly
  100000-bit gap there. Do not restart that refinement as the main plan.

Details: [all-run closure checkpoint](coded_state/QUARTET_ALL_RUN_CLOSURES.md).
Next evaluate the closure-restricted remainder with broad tag counts,
then retain joint rank/availability/closure classes rather than interpreting
upper-envelope sensitivities as physical trajectories. Construction fixed;
local seven-type refinements remain parked. No benchmark, production,
paper, commit or push changes; jobs stay serial, guarded and file-free.

### 2026-10-31 — Joint diagnostics corrected; global run-end cost retained

Previous turn was progress. This continuation also makes proof progress;
the full fixed-construction goal remains active and incomplete.

- Instrumented the global transfer. Coarse marker perturbations crossed
  different dominant contributions, so their inferred combined trajectory
  was misleading. Retracted the supposed h=512 forced-birth quota excess.
  Tiny perturbations identify an eligible-degree lower-bound violation
  there; the h=2048 profile does not violate aggregate capacities.
- Wider quota tilts, lower/upper eligible-degree markers, and joint rank
  caps were implemented and tested. They do not close the hard classes.
  An actual-active-row rewrite was ruled out as negligible: summing that
  decomposition reconstructs the existing shortened-row rank moment.
- Derived an accepted-map singleton propagation filter. Its positive
  258-state factorization and free terminal passed two audits and 34 exact
  comparisons. The coarse version remains a separate failed control.
- Proved that canonical forced births start distinct availability runs.
  At most 512 such runs remain open at the free terminal. A joint run
  argument gives state sum <=e^-512(ce)^f lambda^(L-f), where
  e=min(1,lambda/(q+1)). This fits the existing 129-state transfer, using
  m=z*b when e<1. Two independent audits and 960 exact rational run checks
  pass. It must not be combined with a duplicate reset penalty.
- The new global control improves selected h=2048 and h=4096 bounds to
  +117834 and +289604 bits. A bounded coarse search also finds +283412
  for h=4096. These remain positive FLOAT upper bounds; no whole-code or
  row-class certificate is claimed.
- All previous helper checks pass. Jobs remained serial and file-free
  on Peach; no construction, production, paper, benchmark, commit or push
  changes were made. Local refinements remain parked.

Details: [global run-cost checkpoint](coded_state/QUARTET_GLOBAL_RUN_COST.md).
Next count joint rank/availability boxes with their own conditional field
bounds, or retain additional availability closures. Keep all earlier
certificates and the invariant-plateau control; preserve the full target.

### 2026-10-31 — Rank-gap class proved; forced activation retained globally

This continuation made proof progress. The full fixed-construction goal
remains active; it is not complete or blocked.

- Inserted the previous minority exclusions into the global transfer.
  Exact small checks and independent audits pass, but these markers do
  not improve the tested middle classes.
- Found the old relaxation's impossible 1->512->1 coverage cycle. Newly
  full quartets without incoming coverage contain forced births, which
  cannot disappear in the next update. Derived audited 307-state and
  129-state controls retaining this propagation. The latter uses only
  degree-512 univariate polynomials and complete regional composition.
- Proved a new global class: R-L>=65537 has failure below 2^-63, with L
  the number of available state positions. The exact rational upper is
  (1048576/1011747)*2^-64. Only the 2^N footprint union is needed; all
  state/output supports were already grouped. Two independent audits
  confirmed the probability space, projective factor, and class marker.
- Applying the new remaining-event constraint improves the selected
  4096-row bound from +454173 to +411663 bits. It remains positive, as
  do the selected 512- and 2048-row bounds. These are floating upper
  bounds, not counterexamples or certified margins.
- All helper checks pass, including the width-8 forbidden-cycle control,
  exact terminal/occupied-region identities, and rank-gap arithmetic.
  Jobs remained serial, bounded, single-threaded and file-free on Peach.
  No construction, production, paper, benchmark, commit, or push changes.
  The frozen pair-helper hash is unchanged.

Details: [rank-gap checkpoint](coded_state/QUARTET_GLOBAL_RANK_GAP.md).
Next isolate the joint rank/available-position classes behind the global
positive bound, separating fractional-moment slack from geometric
overcount. Keep the invariant-plateau control and all earlier certificates.
Do not restart local seven-type refinements or change the construction.

### 2026-10-31 — Global composition reduced; covered-footprint classes closed

Previous goal turn: progress, with the expansion theorem and audited
grouped-footprint reduction. This continuation preserves the full target.

- Factored the forced-footprint kernel through130 nonnegative states.
  The implementation composes all4096 updates and full16-update regions,
  preserves the exact empty-input row and free terminal, and subtracts
  the zero-coverage operator through a positive recurrence.
- Combined it with a fixed candidate eligible-row set and quota-correct
  regional conditioning. Three- and17-tilt diagnostics sharply reduce
  some earlier gaps; a64-row whole-history proposal is negative. The
  middle classes still have positive bounds. No floating proposal is
  recorded as a certificate or full row-count cover.
- Passed261 exact factorization-entry comparisons with explicit small
  quartet enumeration,8 composed paths,25 eligible-rank substitutions,
  and6 rank-one renewal checks. Numerical jobs remained serial and bounded.
- Independently identified a method obstruction: an unsplit theta=1/4
  moment containing the preserved invariant plateau exceeds2^33152 for
  every output marker. Smaller lower-footprint tilts remain possible;
  this does not refute the code or the high-coverage theorem.
- Extended the burst theorem from binary STATE supports to binary COVERED
  footprints, permitting arbitrary internal states. A run generating
  function reproduces the previous scalar costs; projective normalization
  absorbs the extra factor. All4095 exact branches remain valid and give
  2^-148. Two independent audits and256 rational run checks agree.
- Proved a reverse-hole identity converting footprint loss F-L to quartet
  suppression. The resulting output-assisted all-footprint class
  m_V+11max(0,L-W)>=32d+7217 has failure below2^-76. Independent audit,
  513 quartet checks,505 path counts, and4096 exact reverse geometries pass.
- No setup, construction, production, paper, benchmark, commit, or push
  changes. Only small source/proof notes were added locally; numerical
  runs were file-free on Peach, single-threaded, below180s/2GiB limits.

Details: [covered-footprint theorems](coded_state/QUARTET_GLOBAL_COVERAGE_RUNS.md)
and [global transfer](coded_state/QUARTET_LOW_RANK_TRANSFER.md).
Next impose the certified-class exclusions in the130-state calculation.
The remaining event has partial coverage, F<53N/64, R<1009256 (h<=7884),
and m_V+11max(0,L-W)<32d+7217. Retain these jointly with actual row-degree
classification rather than recharging the closed classes. The full
45%/40-bit certificate remains open; local seven-type work stays parked.

### 2026-10-31 — All-history expansion class; low-rank remainder still open

The full fixed-construction target remains active. This continuation made
proof progress, not a whole-code closure or a construction change.

- Proved the exact uniform-quartet moment E[4^-G]<=2^-minority at all513
  incoming weights. Combined expansion with state deaths and summed every
  support history analytically, including the free terminal state.
- The new theorem closes the event
  m+32z-ceil(log2 binom(I,z))>=32d+7216 below2^-76 over the actual accepted
  setup. The setup/order-of-averaging audit passed. The bound concerns a
  J=GF(2^32) support class; its labels do not automatically survive a
  larger-alphabet basis projection.
- Derived D<=E<=D+3m and a sharper activity/death history count. Selected
  dense partial-state profiles pass a joint floating bound, but a broader
  3392-point scalar screen still fails, often at E=D. No point screen is
  represented as a cover or certificate.
- Implemented both accepted-map support tests inside the513-state global
  kernel. Regional conditioning leaves the empty path exact and charges
  its quota correction only on positive-footprint regions. Entrywise
  minima across four tilts improve larger-row probes substantially, but
  do not close: selected h=512/2048/4096 bounds remain roughly
  +164k/+586k/+874k bits. The h=1 tested bound is still+224 bits.
- Identified the remaining relaxation explicitly: the kernel retains a
  common eligible row set and its regional quotas, but drops the actual
  individual covered-degree requirements defining that set. Aggregate
  coverage markers did not improve the tested one-step potential. A
  valid rowwise Poisson-binomial filter is derived, not yet composed.
- Made a second global reduction: fix the footprint and sum ALL compatible
  state histories before the fractional moment. Their common RS rank
  permits this grouping. A two-state positive-vector bound gives the
  grouped factor (1+(q+1)x)^L. Two independent audits agree, including
  projective normalization and the free terminal state.
- Retaining covered positions in partial quartets as forced births gives
  a global kernel with forced/optional/uncovered labels. A further safe
  strengthening marks full quartets with no incoming coverage as forced
  too. The kernel is derived and audited; the full regional composition
  is not yet evaluated. Its coarse product-potential probes still fail.
- Proved a routing-variance interface using actual covered sizes:
  variance<=D_V/1024+1/2+Q/8192. Dense coverage with state holes keeps
  Q=0 while still paying field suppression. The new helper passes504
  exact chain checks,106 footprint-kernel comparisons, and1024 forced-
  birth grouping checks. These validate the reduction, not the distance.
- Numerical work stayed serial, bounded, single-threaded and file-free
  on Peach. No benchmark, production, paper, commit, push, or setup change.

Details: [global expansion checkpoint](coded_state/QUARTET_GLOBAL_EXPANSION.md)
and [grouped-footprint reduction](coded_state/QUARTET_GROUPED_FOOTPRINT.md).
Next compose the stronger forced-birth footprint kernel with temporal
coverage and row classification, using a dense burst plus sparse fringe
as the structural control. Park local refinements and avoid repeating
common-potential tuning. The whole45%/40-bit certificate remains open;
all earlier checkpoints are preserved.

### 2026-10-31 — All-row burst class closed; regional operator diagnosed

Previous goal turn: progress, with the rational53N/64 coverage theorem.
This continuation preserves the full target and parks the local refinements.

- Implemented the exact513-state quartet kernel and16-update regional
  composition. Twenty renewal comparisons and106 independent kernel
  comparisons pass. A renewal identity enforces at least129 occupied
  regions without an extra large matrix counter.
- Fixed an underflow in the first full-width diagnostic: ordinary scalar
  coefficients vanished, so their negative-infinity results were rejected.
  Normalized vectors and log-domain convolutions give positive bounds.
  Four fixed-marker and four class-specific-marker probes do not close.
  This is a diagnosed bound loss, not a construction counterexample.
- Derived an all-row analytic proof for every full/empty-state history,
  with no2100-inactive-row restriction. It combines projective selected-zero
  counting, exact burst multiplicities, regional quota coupling, and
  a shifted positive-part MGF using the actual small coverage variance.
- Exact rational checks cover4095 branches and give a total failure below
  2^-148 at45% distance; an independent audit agrees. Worst low-coverage
  branch is r=7; worst high branch
  is r=43. The earlier high-coverage and rank-class theorems remain intact.
- Root-owned remote runs were serial, bounded, file-free, and guarded.
  No benchmark, production, paper, construction, commit, or push changes.

Details: [global burst theorem](coded_state/QUARTET_GLOBAL_BURSTS.md).
Remaining: partially occupied state histories outside the previous rank
and coverage classes. Next combine their suppression/expansion cost with
the burst proof and actual regional variance, rather than refine a loose
common scalar potential. The full45%/40-bit goal remains active.

### 2026-10-31 — Global low-rank goal active; high coverage certified

The user explicitly requested a goal to pursue the global rank bound,
park local refinements, and keep the construction fixed. The active goal
retains the full 45%/40-bit target, not a restricted-class substitute.

- Sum optional output supports within a state history before applying
  the fractional moment. Projective normalization saves one field
  dimension; the survivor factor is now (1+m^(1/theta))^theta rather
  than 1+m. Three reviewers contributed independent derivations and
  structural checks; only the root ran numerical work.
- Proved an all-history class theorem for |V|>=53N/64=1736704, with
  failure strictly below 2^-6359. Its rational route and history bounds
  and final integer power comparison pass exactly. Both it and the
  preserved 7N/8 bound below 2^-56023 were independently reviewed.
- Derived the exact quartet transfer with incoming weight and live/birth
  footprint counters. Passed 212 exhaustive small all-destination checks
  and 1536 exact core/fringe inequalities.
  Floating scalar probes improve substantially but do not close the
  smaller-footprint range; no floating result is labeled a certificate.
- Identified a real rare invariant-subcode setup as a mandatory control
  against uncapped first moments. Derived a self-contained core/fringe
  bound retaining distinct forced-birth rows and exact sampling without
  replacement. Accepted-map filters remove raw witnesses that require
  a rejected minor.
- All remote runs were sequential and file-free, with a 180-second
  watchdog, 2 GiB address-space limit, and output-volume guard. The first
  probe found scipy unavailable; the helper now uses numpy only. No
  installation, benchmark, construction, library, paper, commit, or push.

Details: [global low-rank checkpoint](coded_state/QUARTET_GLOBAL_LOW_RANK.md).
Remaining: |V|<1736704 together with failure of the earlier global rank
threshold. Next compose the grouped count while retaining regional
concentration, using the exact quartet kernel and the necessary fringe.
The goal remains active. Earlier local checkpoints are preserved.

### 2026-10-31 — Global inverse-rank reset, broad class certified

The user requested a larger proof rethink. Three independent mathematical
reviews examined reductions, global alternatives, and structural obstructions.
They made no edits or numerical runs. The construction remains fixed.

- Derived the direct inverse-mask bound2 q^(j-d-R(V)), including residual
  outer equations on dead coordinates. Their rank is independent of death
  equations because the latter have distinct live pivot columns. All
  reviewers agree; forced nonzero births in short rows are excluded.
- Output support plus death positions determines the state support. This
  gives binom(N,j)binom(j,d) masks and a global profile-sensitive rank
  threshold with failure below2^-63. No active-row union, pending-run
  truncation, or regional stationary assumption is used.
- Exact corollaries certify rank>=1009256 below2^-62 and
  rank+deaths/2>=1009257 below2^-71, for every fixed geometry. These are
  broad classes, not the whole-code claim. Corrected an overly strong
  intermediate2^-72 rounding before recording the theorem.
- Retained capped RS rank through routing. A conditioned row-census bound
  preserves every regional count and the total footprint. Its low-rank
  all-history sum is not yet evaluated or certified.
- Passed65535 exhaustive small-field outputs,624 joint-rank classes,
  255 short-row rejection checks,1397 support-census checks,36 conditioned
  route-census checks, and exact integer budget comparisons. One bounded
  remote check at a time; no native benchmarks or raw files. Frozen
  helpers, production, paper, and prior checkpoints are unchanged.

Details: [global rank reset](coded_state/QUARTET_GLOBAL_RANK_RESET.md).
Next evaluate the remaining low-rank event globally with actual footprint
and forced-birth constraints. Do not claim the high-rank theorem closes
the45%/40-bit minimum-distance target. The prior automated goal stays cleared.

### 2026-10-31 — Bounded regional composition diagnostic; no closure

The user requested a bounded one-region test after clearing the automated
goal. Construction unchanged; no new goal or benchmark was started.

- Added a full-width sequential importance probe retaining all sampled
  seven-type populations, the exact2100-input quota over16 updates, and
  pending types across region boundaries. No population truncation or
  internal terminal test. Estimates vary by hundreds of bits and fail
  exact lower controls; they are rejected, not recorded as certificates.
- Derived a full-width closed-suffix matrix-entry lower bound. Eight fully
  allowed suffix updates give all-Z to all-L2 weighted entry at least2^75,
  proved by rational coefficient evaluation with the actual regional quota.
  This rules out uniform contraction below1 at the retained weights, not
  the global45%/40-bit target. The transition cannot repeat from its source.
- Independently, the empty-paid-set subset gives all-Z row lower2^-442;
  with flat boundary weights it gives1. Both expose missed weighted mass.
- Exact width-four,16-update regions retain75 endpoint populations and
  peak1115 states. One- and two-region log2 sums are12.6466 and18.1430,
  with integer brackets proved exactly. The second region carries the
  actual outgoing distribution rather than repeating the first scalar.
- Passed718 normalization checks,1232 exact suffix checks, and negative
  controls for global-quota substitution and an extra terminal test.
- Source-only, file-free bounded remote work; no construction, library,
  paper, frozen-helper, commit, or native benchmark changes. Two early
  diagnostic sessions briefly overlapped; subsequent checks were serial.

Details: [one-region checkpoint](coded_state/QUARTET_ONE_REGION_TEST.md).
The full-width composable upper bound remains unproved. Next use
deterministic boundary-cell bounds and count every exit; retain the exact
quota and use the suffix entry as a mandatory lower-bound control.

### 2026-10-30 — Exact pair rows and the closed-set regional-quota obstruction

Previous turn: progress. This turn adds exact two-update identities and a
directed obstruction that changes the next calculation. It does not close
the45%/40-bit goal or change the construction.

- Derived exact all-destination pair rows for all seven homogeneous source
  types. The first partition's non-full/target constraints and the second
  partition's blocked/selected incompatibility are retained.22 exact small
  comparisons pass;182 rational evaluations verify all original-count
  coefficients for14 small homogeneous profiles.49 derivative checks pass.
- At the old weights, the L2 row drops from1289.76 to1080.96 log2 units.
  A bounded71-second homogeneous-source optimization still gives a
  positive full-stream diagnostic near638052. Mixed B/L2 bounds are
  improved but do not constitute all-source or all-history coverage.
- Rechecked the old global-quota problem after paid-set closure.3388
  suffix witnesses are exactly closed. Directed lower endpoints exceed
  22406 for log2(b)>=1.5 and136513 at the current log2(b)=2.5, even for
  exact global coefficient extraction and with the payment cap retained.
  These suffixes are impossible under the real regional quotas. This is
  a relaxation obstruction, not a codeword or failure-probability claim.
- Added pair composition with the original-inactive counter reset only
  at real region boundaries and one terminal test. It reproduces all11
  exhaustive small closed-set coefficients. Global-quota and extra-terminal
  negative controls are detected. Peak69 states; no large table launched.
- Serial bounded remote runs wrote no files. No construction, production,
  paper, frozen-helper, commit, or benchmark changes. Prior checkpoints
  remain indexed and the full goal remains active.

Details: [pair/quota checkpoint](coded_state/QUARTET_CLOSED_PAIR_QUOTAS.md).
Next derive quota-conditioned pair bounds for mixed populations and their
exits, retaining the first-step deaths as second-step blocked mass. Do not
run a larger temporal transfer with only a common global marker: at the
current moment parameter that relaxation is now rigorously ruled out.
Other row counts and the final failure-budget audit remain outstanding.

### 2026-10-30 — All-destination factorization and a certified envelope obstruction

Previous turn: progress. This turn makes proof-interface progress and
rejects a global shortcut; it does not increase certified coverage or close
the45%/40-bit goal. The construction is unchanged.

- Factored the exact seven-type quartet kernel into non-full, full without
  a selected predecessor, and full with a selected predecessor cases.
  Every destination is included.420 rational local checks and eight
  width-eight all-destination coefficient comparisons pass.
- Derived a source-population-uniform bound using multinomial mean atoms
  and positive Bernstein degree elevation.8550 small mean-atom checks
  and691 exact normalization checks pass. Five-/seven-category FLOAT
  proposals fail badly; increasing the tested degree does not help.
- Replaced numerical suspicion by a directed AM-GM/flow obstruction.
  Exact stationary flows over132 monomials cancel all product weights;
  a rational mixture exactly matches the regional inactivity fraction.
  The mean paid count is below the allowed cap. Every uniform one-step
  product envelope with common markers fails for1<=b<=8 at the current
  field weights, even with a payment-cap marker. Trusted log2 lower
  endpoint3856048867539029/4294967296 exceeds897806. This is NOT a
  construction obstruction or a lower bound on its failure probability.
- Tested a population-dependent blocked/selected interaction potential.
  Eight rational polynomial checks pass. Homogeneous-row FLOAT proposals
  improve, but mixed blocked/selected sources defeat the sampled choices.
  Do not spend another sweep tuning this particular failed proposal.
- All computations were serial, bounded, remote, and file-free. No
  construction, production, or paper changes; no native benchmarks,
  commits, or raw dumps. Prior certificates and checkpoints are preserved.

Details: [global envelope audit](coded_state/QUARTET_CLOSED_GLOBAL_ENVELOPE.md).
Next use the exact factorization in a coupled two-update population bound.
Retain deaths as blocked source mass in the next step and cover exits,
starting with the pure-selected and mixed blocked/selected hard sources.
Do not replace this with another uniform one-step product norm. Other row
counts, initialization, regional quotas, and terminal acceptance remain
explicit whole-code obligations.

### 2026-10-30 — Directed regional transfers and an exact global score

Previous turn: progress. This turn replaces FLOAT-only local evidence by
directed coefficients, composes full regional quotas, and allows moving
populations. The full45%/40-bit goal remains open; construction unchanged.

- Joint positive extraction passes39 exact recurrences and61 full quartet
  checks. Two conditional proposals improve from331.49/322.29 to163.27/
  177.23 in FLOAT. Those comparisons still fix the131/132 label schedule.
- Added integer/rational upward powers, arbitrary rounded log markers,
  and outward uint64 coefficient recurrences with overflow checks. Eight
  power tests and54 rational comparisons pass. A cheaper directed version
  gives224.982322/230.837126 for the same conditional comparison.
- Summed every quota allocation for the fixed populations: directed
  regional bounds268.207363/274.473345. Original inactivity remains distinct
  from paid count; there is no internal terminal reset.
- Built directed polynomial matrices for9 and25 population states. All
  internal paths are included, but exits are not yet covered. Exact tests:
  28 matrix comparisons and53 modular coefficient comparisons.
- Retaining incoming L1 modulo4 gives9-state bounds261.880049/278.324275.
  Widening the second neighborhood to25 states and modulo8 gives291.230705.
  This is not a passing global certificate; the larger range is harder.
- A different25-state neighborhood near a score-compatible profile has
  directed weighted operator log2 bound-28.661259882. Boundary potentials
  are positive rationals and must be carried between regions.
- Derived and tested2245 local instances of the exact global condition
  sum(707*S-1905*L0+3072*Z-403200)+1905*S_terminal>81074176.
  This retains the one free terminal correction. The earlier low-score
  25-state neighborhood has maximum per-step score-19886, so trajectories
  confined there cannot fail. They can occur within bad trajectories that
  compensate elsewhere; do not discard those mixed paths.
- Six higher-score stationary probes have negative FLOAT conditional
  bounds. This guides coverage; it is not a proof over that entire range.
- All runs were serial, bounded, remote, file-free. Longest selected run
  63s, wider modular run28s. The first launch exceeded Windows command
  length before starting remotely; a smaller in-memory loader succeeded.
  No construction/library/paper changes, commits, or generated dumps.

Details and exact endpoints:
[directed regional checkpoint](coded_state/QUARTET_CLOSED_DIRECTED_REGIONS.md).
Next: use the telescoping score to cover transitions between population
sets, especially exits and compensation between low- and high-score
regions. Preserve all other row-count and final-budget obligations. Do not
substitute the current restricted regional bounds for full closure.

### 2026-10-30 — Population bounds expose the necessary run history

Previous turn: progress. This turn derives and tests population-conditioned
bounds, rejects a coarse regional transfer, and obtains an exact seven-type
interface. Full45%/40-bit closure is still open. No construction change.

- Positive population coefficients:60 exact checks and60 tilt replays.
  Original-inactive coefficients:480 exact comparisons and480 tilt replays.
  Uniform source cells:585 small comparisons. Local improvements are not
  global margins; full-width numerical values are FLOAT proposals.
- Built45 coarse cells with all513 original-inactive labels and nine
  outgoing marker pairs. An exactly quota-compatible path gives+622862.10
  in this upper relaxation. This rejects that fixed proposal, not the code
  or exact moment. Runtime112.73s. No further sweep of this coarse matrix.
- Joint endpoint populations remove false diagonal transitions. If U=S at
  both endpoints, selected support cannot shrink. A constant population
  requires its zero count divisible by4 and has an exact routing factor
  binom(128,z/4)/binom(512,z). All275 small occupancy checks pass.
- A blocked pending coordinate and a selected pending coordinate require
  opposite fullness conditions and cannot share a quartet:446 exact local
  checks. The old type16 count equals current deaths in6355 local rows.
- Derived an exact quotient with seven types and210 quartet rows. Every
  one of1365 original local rows, all twelve terminal types, and all eleven
  small regional coefficients agree. No histories were maximized away.
- Exact age-flow identities hold on680 quotient rows. They show that
  constant-population trajectories require l0=l1=deaths, selected>=3*l0,
  and blocked>=l0. Two very negative numerical fits violate these exact
  conditions and are excluded, not counted as positive proof evidence.
- Seven-type fixed-population regional proposals now include-83.47 and
  -414.85, but neighboring profiles consistent with age flow remain
  +331.49 and+322.29.
  These omit other histories and the final global factors. None is a full
  certificate or a lower bound on failure. The exact quartet polynomial
  has2245 coefficients and is ready for sharper joint extraction.
- All runs were serial, bounded, remote, and file-free. New source only;
  no library/paper changes, commits, or generated dumps. The frozen helper
  is preserved.

Details: [population checkpoint](coded_state/QUARTET_CLOSED_POPULATIONS.md).
Next: sharpen the two positive seven-type profiles consistent with age
flow, then retain age flow in the regional transfer. Preserve quotas and the single terminal
test. Cover moving populations, other row counts, and directed final sums
before claiming closure. Do not return to the obstructed all-paid-set
moment or the rejected age-forgetting coarse transfer.

### 2026-10-30 — Exact used-payment transfer and a necessary marker change

Previous turn: progress. This turn builds the replacement counting
interface and tests two possible shortcuts. The full45%/40-bit goal is
active and open; construction and prior failure allocation unchanged.

- Derived a six-state automaton with exactly one accepted canonical support
  per allowed stream. Its pending-input form has12 live types. Tested all
  2046 nonempty bit streams through length10 against direct run selection.
- The quartet factor is x^u(1+xy)^v, where u counts disallowed inputs and
  v counts actually touched inputs. Original inactivity and payments are
  now coupled to coverage, not chosen from all allowed inputs.
- Implemented the exact width4 two-region transfer. It reproduces all11
  exhaustive paid-count/score coefficients with131 maximum states; paid
  count zero has the required binom(8,2)^2 normalization. Terminal acceptance
  is applied once and pending states cross the regional boundary.
- Derived the full-width random-partition kernel as a coefficient of a
  quartet polynomial power, divided by the incoming multinomial count.
  Width8 labelled partitions and ordered-group composition agree on126
  coefficients;45 additional coefficients check the free terminal. This
  is an exact interface, not an evaluated full-size certificate.
- Tested a uniform scalar quartet envelope over919 nonzero local rows and
  20723 marked monomials. Its FLOAT upper candidate is+912692.56, despite
  a small optimized product potential. The worst row is the initial-zero
  quartet. Do not infer a structural lower bound from this failed proposal;
  retaining finite populations and regional coefficients is the next step.
- Directed closed-burst counting removes over101214bits of redundant
  payment count from the seven-update witness. Longer bursts still rule
  out b=16. Nine-update halves give the necessary restriction
  log2(b)<2.88903170 at the selected field weights for a40-bit target.
  This is a moment-method restriction, not a code-distance counterexample.
- Runs were serial, remote, bounded, file-free:0.29s regional tests,2.39s
  partition tests,0.69s envelope proposal. The first optimizer import found
  SciPy unavailable; no packages were installed. A NumPy-only proposal
  replaced it; the reported bound uses an explicit maximum and does not
  assume optimizer convergence. No production/paper changes or commits.

Source and exact scope: [closed-set transfer](coded_state/QUARTET_CLOSED_TRANSFER.md).
Next: bound the population-conditioned kernel while retaining the original-
inactive quota across16 updates. Begin with the hard middle class and a
moment below the proved ceiling, then cover other paid/row counts and do
directed verification. Do not return to the obstructed all-paid-set moment
or spend a large sweep on its scalar replacement. All checkpoints preserved.

### 2026-10-30 — Exact moment obstruction; count closed paid sets instead

Previous turn: progress. This turn identifies a proved obstruction to the
current proof interface and derives its replacement. Construction, original
full45%/40-bit objective, certified exclusions, and failure budget unchanged.

- Completed marked arbitrary-pair regional composition. Exact tests pass
  5979 moments,879 FLOAT comparisons,6732 bucket envelopes. Reused201 coarse
  source caps over2502 fine boundary states;1,405,231 pair edges tighten.
  Runtime729.66s. FLOAT upper candidate+255263.80682 and relaxed lower
  diagnostic+212380.35033 still fail. The previous upper+255164.00196 remains
  applicable because the matrix only decreased; weight fitting is not exact.
- Found a feasible full-input burst crossing each pair of regions. Seven
  full updates on each side support thirteen selected updates plus a real
  death. Total A=851968,D=65536, without any internal terminal reset.
- Directed rational counting proves the EXACT all-paid-set moment, after
  its row-set union, exceeds2^513.62928 for every b>=1 and every current
  mu,nu,eta>=0. This does not require onset three. Thus more exact histograms,
  finer buckets, or scalar parameter tuning cannot close this interface.
  It is not a lower bound on code failure. Eight-update halves give a
  stronger obstruction for the selected weights and onset-three family.
- The obstruction check passes28 exact event counts, three support sizes,
  rational multiplier inequalities, and directed integer/rational logs.
  Runtime0.46s. A preliminary Jensen lower diagnostic is negative(-14680.29)
  because it discards cross-boundary rewards; its3130 exact mean checks
  and negative result are preserved separately.
- Proved a canonical closure: delete paid positions not touched by the
  maximizing support. Its support and H value are unchanged, its size only
  decreases, and it is fixed after one reduction. Since B(c) increases,
  every bad event still yields a closed witness. The original lower bound
  on c must not be imposed on that witness; the upper bound survives.
- Derived the exact regional factor binom(v,c)/binom(8192,r), with v the
  actual touched-input count, replacing binom(a,c)/binom(8192,r). Exhaustive
  toy checks pass12544 closures; summing65536 allowed masks reproduces all
  11 paid-count/score coefficients for4592 closed pairs. This is not yet
  a full-width probability calculation or certificate.
- Runs were serial, remote, bounded, one-thread, and file-free. A Windows
  command-length failure occurred before launch; unused in-memory source
  bodies were removed before retrying. A local diagnostic's positive
  underflow was repaired using an explicit conservative floor. Frozen
  helper hash unchanged; no library/paper edits, commits, pushes, promotion.

Sources: [exact obstruction](coded_state/QUARTET_EXACT_MOMENT_OBSTRUCTION.md),
[closed paid sets](coded_state/QUARTET_CLOSED_PAID_SETS.md),
[preserved pair calculation](coded_state/QUARTET_ARBITRARY_PAIR.md).
Next: derive the local delayed-coverage transfer for canonical closed sets.
Retain the actual quartet incidence, regional quotas, and one free terminal;
test the crossing-burst family before a large calculation. Do not keep
optimizing the superseded all-paid-set moment. The full goal stays active.

### 2026-10-30 — Couple two updates before discarding the third age

Previous turn: progress. This turn makes a bounded dependency repair, not
a construction change. The full45%/40-bit goal remains active and open.

- Derived two-update coefficients when either input count is512. That
  count refers to the allowed-input set, not actual nonzero state support.
  Exact tests pass785 coefficients,777 comparisons with FLOAT uniform
  bounds, and15 blocked-counter identities preserving regional quotas.
- Composed the tightened pairs over all paths on the same eight selected
  counts and2502 states. FLOAT upper+255164.00196 versus+263756.41563
  previously; balanced relaxed lower+216759.24666. Runtime302.90s. This
  diagnoses insufficient improvement, not actual low-distance codewords.
- Generalized the calculation to arbitrary consecutive counts using a
  bivariate positive polynomial. It retains the actual intermediate third
  age without introducing a third boundary-state index. Exact tests pass
  5979 marked moments;879 FLOAT bounds dominate exact toy moments.
- All12 selected full-width total-row bounds improve. The408->312 rows
  improve by15.63 to56.88bits, while408->408 improves only about1.5bits.
  These are local coefficient changes, not accumulated security margins.
  Tests plus local diagnostic take13.01s; no full regional run of the
  arbitrary-pair formula is claimed yet.
- Runs were serial, remote, one-thread, bounded, and file-free. Prior
  controls and frozen helper hash are unchanged. No library/paper edits,
  commits, pushes, promotion, or new failure allowance.

Sources and scope: [full-input pairs](coded_state/QUARTET_FULL_PAIR.md),
[arbitrary pairs](coded_state/QUARTET_ARBITRARY_PAIR.md). Next turn: use
outgoing-age markers to form arbitrary-pair edge bounds and compose the
matched restricted regional calculation. This directly tests the remaining
partial/partial dependency rather than retuning the failed transfer.
All-count/row/covered-count coverage and directed arithmetic remain open.

### 2026-10-30 — Regional weights do not fix the tail relaxation

Previous turn: progress. This turn tests the proposed regional weights and
separates their quality from an intrinsically large relaxed coefficient.
The full45%/40-bit objective remains active; construction and prior controls
are unchanged. No restricted result is substituted for full coverage.

- Intersected the original edge bounds with25 joint moment bounds. Fitted
  a second positive vector to the quota-weighted sixteen-update operator.
  The true initial region is evaluated directly;255 subsequent regions
  use a maximum row ratio, with the free terminal applied once.
- The165-cell all-count calculation settles at FLOAT+407961.58605, compared
  with+419478.71616 before regional fitting. It still fails. A new balanced
  coefficient lower diagnostic is+404094.84018. Its exact formula retains
  every ordering of44 q=-1 regions among255, after an initial q=0 region.
  This is strong numerical evidence against fixing THIS linear transfer
  by further boundary/payment tuning, not a lower bound on actual failure.
- Narrowing only incoming count cells64->16, with age buckets still64,
  gives525 states and+382506.65796 upper/+377771.82876 lower. This removes
  only a modest part of the gap. Runtime323.66s; all513 counts retained.
- An exact selected-count diagnostic uses0,312,384,392,407,408,424,512 and
  width16 age buckets:2502 states,+263756.41563 upper/+221544.01084 lower,
  runtime111.39s. This adds extremes to the passing-band control but does
  not cover other counts. It isolates a remaining tail relaxation.
- Tested one old-age potential that equalizes selected type3/type4 slopes.
  It changes proof weights only. The restricted upper worsens to+274462.68;
  lower diagnostic+249420.52, runtime110.24s. Preserve it as a failed test,
  not an improved construction or an optimum over all proof weights.
- Exact toy checks pass nine regional upper bounds and their numerical
  replays, nine balanced coefficient lower bounds, and an all-covered-count
  generating-function identity. The latter supports an all-c route: at
  log2(payment marker)=-24/5 the threshold's c coefficient cancels. This
  could avoid a separate union over c, but the full bound is not evaluated.
- All runs were serial, remote, one-thread, bounded, and file-free. Coarse
  runs169.72/167.46s;1200s watchdog. Counter floors are tracked separately;
  their reported output gaps exceed288bits in every retained run. The
  current source allows32 propagated floors for row-sum tolerance. A
  Windows command-length failure occurred before remote launch; unused
  source modules were removed from the in-memory loader before retrying.

Sources and scope: [quota-potential checkpoint](coded_state/QUARTET_QUOTA_POTENTIAL.md).
Next test a two-update bound that retains the actual third-age count. The
current transfer maximizes the type3/type4 split independently, although
the next n3 equals the preceding selected type2 count. Test this dependency
on a tail transition before building a much larger grid. Full count/row
coverage and directed verification remain required; no redesign, new
failure allowance, production/paper change, commit, or push.

### 2026-10-30 — Compose total caps; retain joint outgoing moments

Previous turn: progress. This turn completes the all-count capped diagnostic
at r2100,c111572. The full45%/40-bit goal remains active and open; the
construction, failure allowance, and earlier checkpoints are unchanged.

- Implemented the row-cap support function by its positive dual, not by
  scaling down individual edge bounds. Exact tests pass81 dual evaluations,
  81 two-region quota cases, and170 comparisons of the production counter
  with exact capped coefficients. The true initial row and one free terminal
  are retained in the full calculation.
- The165-cell,7725-edge all-count run completes in750.39s. Its FLOAT log2
  expression is+429554.70068170864, improved from+688051.0077459727 but still
  failing. It covers all counts0..512 and all regional payment allocations
  only at the selected r,c. It is not a whole-code or directed certificate.
- Its proposed vector residual is1.00668; the actual maximum ratio is used.
  The worst regional source cell is a512,n1 in128..191,n2 in64..127. The
  old linear counter-floor error estimate is not claimed for this operator.
  All added floors remain in the reported expression.
- Derived a joint moment bound using selection weights
  (d0,xy*d1,x*d2,x*d3,x*d4). It bounds outgoing first/second age cohorts
  together with the fixed population. Checks pass3393 exact identities,
  3555 exact destination bounds, and3393 FLOAT coefficient comparisons.
  A25-marker local screen lowers the l408 edge-sum log bound from477.17620
  to430.57368. This is local, not a global margin; full composition is open.
- An activity potential kappa^a telescopes exactly along actual paths.
  Checked966 transitions and1250 paths. It is an optional proof-weight
  refinement, not new encoder randomness or a claimed distance gain.
- A fixed-path diagnostic inside the capped relaxation gives-27937.69,
  with exact regional payment counts. This negative diagnostic neither
  closes the transfer nor establishes an intrinsic positive obstruction.
  It cautions against declaring weight optimization exhausted.
- Runs were serial, remote, one-thread, bounded, and file-free. Full run:
  1200s watchdog,2GiB array bound,512MiB counter guard; observed RSS~114MiB.
  Joint local screen10.73s; fixed-path diagnostic27.86s. A lost earlier test
  result was replayed only after confirming no remote computation remained.
  No local numerics, native benchmarks, paper/library changes, or commits.

Details and sources: [capped regional checkpoint](coded_state/QUARTET_CAPPED_REGIONAL.md).
Next test a boundary potential fitted to the full sixteen-update regional
quota. The current proposal comes from unrestricted updates, which admit
long full-input runs forbidden by that quota. Use the coupled joint-moment
edge bounds for a cheaper proposal, then validate the resulting vector with
the actual upper operator. Do this before a much larger state grid. Other
covered/row counts and directed whole-code verification remain required.

### 2026-10-30 — All input counts, and a coupled-population repair

Previous turn: progress. This turn evaluates complete input-count coverage
at the selected r2100,c111572 and diagnoses why the coarse bound fails.
The encoder and full45%/40-bit objective remain unchanged.

- Implemented165 incoming cells and7725 exact outgoing count labels covering
  every count0..512. The backward transfer includes every count path,
  all regional payment allocations, the true initial row, and free terminal.
- A reused narrow-band normalization gives FLOAT+799982.94. Retuning the
  count/payment markers gives+688051.0077459727, still far too loose. The
  upward counter-floor contribution is over750bits below the regional norm;
  it is not the dominant loss in the retuned calculation. The final replay
  uses a common tilted-matrix log scale and reproduces the same result.
- Found an exact coarse-count defect: independent maxima for inactive and
  active populations can count575 positions instead of512. For incoming
  a320..383, outgoing l408, selected k208, one unweighted relaxed term is
  at least2^40 times the entire true unweighted row. This is a proof-bound
  obstruction, not a low-weight codeword or evidence against the encoder.
- Derived a coupled coefficient cap keeping n0+...+n4=512. At the diagnostic
  incoming cell it reduces the l408 total-row log bound by48.19bits; other
  selected counts show1.01--86.01bit local reductions. These are local
  factors, not improvements that can be multiplied across4096 updates.
- Derived the exact row-cap support function: maximize continuation weight
  subject to individual edge caps and the coupled total cap. Equivalently,
  minimize theta*g+sum m_e max(f_e-theta,0). This retains the outgoing
  cohorts without rescaling edge bounds unsafely. Regional composition
  of this constrained operator is not yet implemented.
- Checks pass8692 exact coarse edge bounds,81 exact regional quota pairs,
  10107 interval polynomial checks,7725 full-width edge labels,4805 coupled
  coefficient checks,180 exact vertex objectives, and81 exact cap programs.
  The population-mismatch lower bound is exact integer arithmetic.
- Full runs are serial, remote, one-thread,300s/2GiB bounded and file-free.
  Final full replay116.60s; local cap screen5.31s; audit0.33s. A Windows
  command-length error occurred before remote launch; the loader now omits
  unused helper bodies. All numerical jobs joined; no concurrent jobs.

Details: [coupled row caps](coded_state/QUARTET_COUPLED_ROW_CAPS.md).
Next compose the row-cap operator through complete regions on the current
all-count cells. Test this specific population defect before allocating a
large finer grid. Preserve the-3907.23 fine-band checkpoint. Other covered
counts, row counts, and directed verification remain required. Goal active;
no external blocker, construction change, or new failure allowance.

### 2026-10-30 — Fine joint cohorts clear the previously failing band

Previous goal turn: progress. This turn first joined its verified live
full-band job; it also derives and checks the interface needed for tails.
The construction and complete45%/40-bit objective remain unchanged.

- The full392..424 input-count band with width16 joint buckets and11618
  states completes at FLOAT log2 expression-3907.2296673407545. This sums
  all4096-step paths inside that band and all regional payment allocations
  at r2100,c111572. It preserves inherited cohorts and the free terminal.
  It improves the coarse joint+4184.67 and fine one-cohort+7566.84 screens.
- The computation takes967.50s, with1200s watchdog and8GiB conservative
  array budget. Observed RSS was about3.1GiB. It uses the maximum row ratio
  of a positive vector, not an assumed converged eigenvalue. Remote,
  serial, single-threaded, no generated files; the job is joined.
- Derived an exact two-neighbor recurrence for the routing coefficients
  A_(a,l,u). It also proves A_(a,l,u)<=A_(max(a_lo,u),l,u) when a>=a_lo.
  The max with u is necessary to retain feasible coefficients above a_lo.
- Derived a rectangular labelled-edge transfer: incoming states can group
  input counts, but every edge retains its exact outgoing count. Regional
  coefficient extraction remains exact; no quota is replaced by a midpoint.
- Exact tests pass680 routing recurrence checks,3770 coarse-count geometry
  checks,324 labelled-edge paths, and the prior square/joint/regional tests.
  Another12023 FLOAT coefficient regressions compare the recurrence with
  the direct sum at counts including0,1,511,512. Audit runtime0.43s.
- The joint driver now uses the equivalent coefficient recurrence. The
  four-count width16 result replays within1.4e-8 in log2 at-9577.53846.
  This is a proof-evaluation optimization, not an encoder timing claim.

Sources: [joint-cohort result](coded_state/QUARTET_TWO_COHORTS.md) and
[exact edge counts with coarse states](coded_state/QUARTET_LABELLED_EDGES.md).
The full-band pass is FLOAT, not a directed restricted certificate or a
whole-code certificate. Any path leaving392..424 remains uncovered, as do
the other covered counts and row counts. Next construct the mixed-resolution
tail operator using exact edge counts; do not merely widen another selected
path family or change the encoder. Directed verification follows a complete
coverage expression. Goal active; no external blocker.

### 2026-10-30 — Compose joint age cohorts through regional quotas

Previous turn was progress; this turn extends the local identity to a
uniform joint-bucket operator. The encoder and full45%/40-bit goal remain
unchanged. No additional failure allowance is introduced.

- Incoming labels retain input count a and intervals for n1,n2. For each
  exact n1 in its interval, sum outgoing joint coefficients before taking
  a coefficientwise maximum. The old coefficient bound is uniform over
  the permitted n2 interval and both remaining age types.
- Positive convolutions reuse the selected-fresh factor across outgoing
  input counts. A cancelling coefficient tilt controls dynamic range;
  assertions reject underflow instead of discarding positive terms.
- Reused the sixteen-update counter DP with the joint kernel. Both cohorts
  cross regional boundaries; only the regional input counter restarts.
  Initial and terminal factors occur once, and all payment allocations
  are summed at covered count111572 for inactive-row count2100.
- Exact tests pass8692 outgoing bucket comparisons and all81 two-region
  input-total pairs. Another1155 exact toy polynomial coefficients pass
  the FLOAT coefficient implementation check. These are not full-width
  directed certificates.
- On the four-count domain {392,407,408,424}, every count path and payment
  allocation is included. The one-cohort width8 expression is+2724.82424;
  the joint-cohort width32 expression is-2697.30787. This is a restricted
  FLOAT improvement, not a whole-code security margin.
- Full392..424 coverage with width32 joint buckets gives+4184.66545, still
  failing. Its marked upper-expression lower diagnostic is+604.29, so
  boundary tuning alone does not repair that particular expression.
- Halving joint buckets to16 improves the four-count control to-9577.53846.
  This identifies substantial discretization loss but cannot be extrapolated
  to the full band. The next test uses width16 on the complete band.
- Old coefficients depend only on a-r. Caching them reproduces the fine
  four-count result exactly and reduces its computation20.93s to16.22s.
  This is not an encoder timing. All runs are remote, serial, file-free.

Details and reproduction sources:
[joint-cohort analysis](coded_state/QUARTET_TWO_COHORTS.md).
The next discriminating test is the fine full392..424 band. All count tails,
other row/covered counts, and directed verification remain required.
Goal active; no external blocker.

### 2026-10-30 — Joint age cohorts and selected-fresh multiplicities

Previous turn was progress. This turn adds a causal local refinement;
the construction and full45%/40-bit objective remain unchanged.

- Derived the exact joint transition n1'=l-k+v0,n2'=v1. Its coefficient
  retains the selected-fresh binomial factor and the remaining old-age
  polynomial, rather than assigning the unresolved age population anew.
- Proved the joint feasibility polygon max(0,l-n1'-a+r)<=n2'<=min(r,l-n1').
  Its vertex envelope sharpens selected two-update paths, but a range
  alone permits zero selected-fresh positions without paying their rarity.
- Factored the joint coefficient into A_(a,l,u), a selected-fresh binomial
  factor, and an old polynomial coefficient at u-v1. This avoids a large
  dense two-cohort matrix for the bounded diagnostic.
- Exact checks pass13200 joint coefficients against the independent
  five-type kernel,13200 Cauchy bounds,565 normalizations, and1349 nonempty/
  2476 empty polygons. The factored and unfactored exact formulas agree.
- At full width, summing joint intermediate buckets improves three local
  two-update comparisons by1.118,6.586,5.969bits against the same first-step
  control. Twelve range-only comparisons also complete. These are FLOAT
  local factors, not whole-stream failure bounds or security margins.
- An initial range run stopped on an unreachable bucket. The fix requires
  its incoming coefficient to be zero before dropping an empty polygon.
  Final tests/screens pass. All runs serial, remote, one-thread,90s bounded,
  below64MiB, and file-free. No frozen helper, construction, library, paper,
  failure-budget, commit, or push changes.

Details: [joint cohort analysis](coded_state/QUARTET_TWO_COHORTS.md).
Next compose this bound across a whole region with the actual quota and
inherited endpoint state. Test the failed392..424 band before expanding
other row/covered counts. Directed verification and all omitted histories
remain required. Goal active; no external blocker.

### 2026-10-30 — Complete regional count paths; diagnose the local relaxation

Previous turn: progress, not a wait or blocker. This turn also makes
progress; construction and full45%/40-bit objective remain unchanged.

- Replacing each fresh-count matrix by a one-update weighted norm fails
  badly (+148295.79 at392..424). Preserved that rejected shortcut rather
  than hiding it or using it to claim that the construction fails.
- Derived a positive-vector regional bound that retains the fresh-count
  distribution through all16 updates. A maximum row ratio is sufficient;
  numerical convergence to a Perron vector is not a proof assumption.
- Summed EVERY input-count path inside each tested band and EVERY regional
  payment allocation, bounding the global coefficient at c111572. FLOAT
  expressions pass400..416 (-12064.53) and396..420 (-3152.64), but fail
  392..424 (+9357.30). This is still a restricted row/covered-count test.
- Halving fresh-count bucket width to4 gives+7566.84, still positive.
  Payment-marker tuning saves less than0.11bits. Twelve regional potential
  iterations on the balanced-payment subproblem also fail to close it.
- Added LOWER diagnostics on the selected upper-transfer expression.
  For the balanced regional allocation family, these give+5130.99 with
  width8 buckets and+5529.10 with width4. In the latter test the regional
  weighted row ratios differ by only0.238bits. This points to local
  histogram relaxation, not boundary maximization. These FLOAT diagnostics
  are not certified obstructions or lower bounds on actual code failure.
- The exact transition also forces outgoing n2=v1, where v1 counts selected
  fresh coordinates. The present envelope discards that next age cohort.
  Retaining it is the next structural proof refinement; no code change.
- Exact tests:64 backward count paths with two terminals; two weighted
  regional sums; a non-eigenvector maximum-ratio check. The full earlier
  exact audit, including8800 outgoing fresh-count coefficients, passes.
  Cached equivalent active-degree coefficients reproduce the old rotation
  result exactly and reduce that screen from7.57s to1.08s.
- Runs were serial, remote, single-threaded,90s bounded, and file-free.
  A proposal-convergence assertion stopped an initial wider run; it was
  replaced by the valid maximum-row-ratio calculation, not an assumed
  eigenvalue. The widest successful run takes44.84s, with arrays bounded
  by768MiB. Frozen helper hash unchanged. No library/paper edits, commits,
  pushes, generated data, or construction/failure-budget changes.

Details: [regional fresh-count analysis](coded_state/QUARTET_REGIONAL_FRESH_TRANSFER.md).
Next derive and exactly test the joint n1/n2 transition, then test whether
it repairs the local wide-band loss before allocating a large transfer.
All count tails, other covered counts, other row counts, and directed
verification remain required. Goal active; no external blocker.

### 2026-10-30 — Count-path entropy and actual fresh counts

Progress, not closure. Full45%/40-bit goal and construction remain fixed.

- Summed bounded input-count bands and every balanced regional-total
  assignment. A refined auxiliary point-mass coefficient bound passes
  400..416 numerically, but392..424 gives+33746.05. Independent per-count
  potential tuning is worse. No dropped/underflowing contribution is used
  to claim success; rejected wider adaptive runs are explicitly skipped.
- Proved the exact outgoing fresh-count identity r1'=l-k+v0. The earlier
  lower bound r1>=a-k loses how many inactive coordinates were selected.
  Derived a bucket transfer retaining r1' and summing k,v0 jointly. Cauchy
  maxima are taken before selecting the positive evaluation point.
- The clustered424/392 diagnostic falls from+14421.31 with the preceding-k
  envelope to-3609.80 with width8 fresh-count buckets. Summing15 rotations
  per region and all44-low/212-high regional assignments gives-3133.02.
  All are FLOAT restricted expressions, not whole-code security margins.
- Exact tests pass8800 outgoing fresh-count coefficients for all width4/8
  histograms, plus the preceding tests.64 exact matrix paths check the
  regional DP;30 rotated words have the required totals and distinctness.
- Sources and compact results are preserved; no raw data, construction,
  library, paper, frozen helper, failure-budget, commit, or push changes.
  Final replays are serial, remote, single-threaded, bounded, and file-free.
  An earlier test launch returned a live session before another screen
  was launched; possible overlap cannot be excluded. Both were joined,
  and the final exact audit and all subsequent screens were run serially.

Details: [fresh-count analysis](coded_state/QUARTET_FRESH_COUNTS.md).
Next sum general regional count paths with the fresh-count transfer. Then
bound omitted count tails, use directed arithmetic, and cover other payment
allocations, covered counts, and row counts. Full closure remains open;
goal active, no external blocker.

### 2026-10-30 — Condition on input history; separate the type coefficients

Progress, not closure. The construction and full45%/40-bit target remain
fixed. The preceding turn supplied the potential interface; this turn
audited it and tested a replacement at full width.

- Proved stationary scalar factors cannot close for ANY positive potential:
  gamma_l>=binom(w,l) and gamma_w>=b^(lambda*w) recreate the preceding
  2^60214 free-suffix obstruction. This is not merely poor h tuning.
- Retained the preceding active-input count and the capped consecutive-full
  run length. These labels impose exact restrictions on incoming types and
  persist across region boundaries. Toy full/empty sequences become exact;
  merely remembering the last count still loses the exit charge.
- Derived a type-separated Cauchy coefficient bound, avoiding the previous
  arithmetic-mean relaxation. Its per-k positive evaluation gives a uniform
  bound on the permitted histograms while retaining actual count labels.
- Tested the fixed407/408-count path at full width. At log2(b)=4, the
  single-column Perron candidate gives previous-envelope lower+375024.36,
  sampled-exact-kernel lower-33128.59, and new-envelope path term-9333.94.
  These are FLOAT diagnostics. A negative single path is NOT a certificate
  or margin: all other count paths and classes still have to be summed.
- Tests:15 stationary floors;58 exact new coefficient bounds;65536 input
  streams/625 count paths;1024 full/empty paths including the terminal.
  Toy audit0.18s; full-width screen2.78s. Serial/remote/bounded/file-free.
  No construction, paper, library, frozen helper, or failure-budget changes.

Details: [history-conditioned bound](coded_state/QUARTET_HISTORY_POTENTIAL.md).
Next sum a bounded neighborhood of the407/408 path while enforcing each
regional total. Decide whether count-path entropy defeats the candidate
before attempting a large directed calculation. If it survives, certify
it and cover count tails, other covered counts, and other row counts. The
full goal remains active, with no external blocker.

### 2026-10-30 — Free regional ends fail; retain state with a product potential

Progress, not closure. Tested the preceding regional shortcut before a
large computation. Full45%/40-bit objective and construction are unchanged.

- A bounded full-width FLOAT read-four screen failed widely, best tested
  log2 expression+328301.86. It is not a directed certificate or a proof
  of universal failure. A separate exact argument gives the diagnosis.
- Proved the EXACT one-region free-terminal product is at least2^60214.26
  at c=111572, for every current multiplier with tau>=2lambda and every
  b>=1. Balanced regional payments and six fully active terminal updates
  already contribute that much. No scalar marker or read-degree loss is
  used. The method grants terminal freedom256 times instead of once.
- Two-region free-terminal blocks also fail at the selected weights: a
  fixed six/seven-update suffix pattern yields score165594.4>B165556 and
  a directed log2 expression lower bound26601.87696. Count assignments
  are fixed to avoid double counting nested suffix events.
- At a real regional boundary the last state bit additionally needs its
  next quartet active. A guaranteed feasible concatenation deletes that
  last update AND pays its death. The correction is lambda+tau, not tau.
- Derived an exact formula for the transfer acting on a product potential.
  Elementary symmetric averaging gives a coefficientwise spatial envelope
  that retains the active-input counter. A lower convex hull reduces
  all-histogram verification to a one-dimensional supremum for each count.
- A stationary potential telescopes across regions and incurs its terminal
  factor only once. The resulting scalar polynomial formula retains each
  regional quota and the total paid coefficient. It has not yet been
  numerically closed; stationary potentials may themselves be too loose.
- Tests:108 exact potential coefficients, all previous12544 quota cases
  and65536 boundary masks, width-eight partition checks and normalization.
  Directed audit24.22s; transfer0.34s; float screen1.17s. Serial/remote/
  bounded/file-free. No construction, paper, library, or frozen helper edits.

Details: [boundary-state potential](coded_state/QUARTET_BOUNDARY_POTENTIAL.md).
Next test stationary product potentials at c=111572, keeping state and
regional counts. Certify any promising scalar suprema and coefficients
with directed arithmetic before widening coverage. The full certificate,
other covered counts, and other row counts remain open; no external blocker.

### 2026-10-30 — Regional quotas are essential even with exact coefficients

Progress, not closure. Audited the proposed global-count transfer before
attempting a larger calculation. The unchanged45%/40-bit goal remains active.

- Proved that even the EXACT uniform fixed-count moment expression is at
  least2^38647.20358 at c=111572, for all current nonnegative field
  multipliers and b>=1. A suffix of825856 active inputs is sufficient to
  exceed every corresponding threshold. Its probability under the relaxed
  uniform-subset law gives the directed lower bound. This strengthens the
  preceding obstruction beyond Bernoulli conditioning or Cauchy evaluation.
- That event is impossible under the actual regional quotas:100 full
  regions require210000 payments, more than111572. The sharp possible
  all-active suffix length under those quotas is440540. This is a
  proof-comparison failure, not a bad codeword or a construction obstruction.
- Corrected the exact five-type transfer: count allowed inputs inside each
  region, multiply by binom(a,a-m+r)/binom(m,r) at its boundary, and carry
  the histogram and total paid count onward. This computes the actual
  expected paid-set polynomial; unweighted coefficients are binom(256r,c).
- Proved a pointwise regional block upper bound and its coefficientwise
  product formula. It removes cross-region state only as an explicit
  relaxation, with worst-case cost(lambda+tau)*512*255. At the selected
  middle profile that cost is133212, so usefulness remains to be tested.
- Checks:5632 nested suffix cases;12544 direct quota/paid-set cases;
  65536 boundary masks; width-eight two-region transfer and coefficientwise
  block checks; all normalizations and terminal weights. Max805 toy states.
  Directed audit24.22s; transfer audit0.34s. Serial/remote/bounded/file-free.
  Frozen helper unchanged; no construction/library/paper edits or commits.

Details: [regional-quota audit](coded_state/QUARTET_REGIONAL_QUOTAS.md).
Next bound the full-width sixteen-update regional polynomial, keeping its
quota and paid coefficients, and test c=111572. If the free-boundary block
relaxation fails, retain boundary state rather than tune the lost global
quota comparison. Other covered counts and row counts remain open. No new
whole-code failure-budget contribution or external blocker was introduced.

### 2026-10-30 — Exact quota audit rules out the Bernoulli shortcut

The preceding turn was progress: it gave an exact window-cover interface
and tight checked allocations. This turn diagnoses a relaxation that would
still defeat its moment calculation, then supplies a count-preserving
replacement. Construction and the entire45%/40-bit objective remain fixed.

- Proved that the middle covered-count111572 has B/lambda<=825555 for all
  current nonnegative field multipliers. Nonnegativity and the all-active
  stream give a lower bound on the EXACT Bernoulli moment expression. For
  every paid marker and moment base, its log2 is at least43307.1183 before
  the row-set union. Thus a sharper spatial/temporal moment cannot close
  this comparison. It is a proof-method obstruction, not a bad-codeword.
- At the old frozen marker1/6, even the actual unconditioned paid polynomial
  contains the wrong-count P=I endpoint, giving log2 lower>133182.4. This
  explains why the negative independent-column diagnostics were misleading.
- Proved an exact density comparison from nested regional unpaid sets to
  uniform global subsets of the SAME unpaid size. Log-concavity makes the
  density maximum occur at balanced regional counts. Directed log2 losses
  at covered counts27476,111572,255463 are495.103,247.360,110.254 or less.
  This comparison is available but has not been charged to a closed tail.
- Also ruled out evaluating the uniform fixed-count moment by one positive
  Bernoulli moment divided by its binomial point mass: its expression has
  log2 lower>38824.4121 throughout the current middle multiplier family.
  The count must constrain the computation, not merely its denominator.
- Derived an exact five-type transfer for onset-three free runs, retaining
  the active-input counter. Uniform quartet occupation and hypergeometric
  type intersections give the kernel; a separate exact terminal factor
  preserves the free final state. This removes overlapping-window and
  independent-column assumptions. The full width512 table would be too
  large and was not launched; it needs a quota-aware analytic envelope.
- Tests:2673 exact nested-quota subsets;26880 direct partition/input kernel
  cases;65536 exact temporal streams;17 explicit count coefficients; full
  binomial normalization through width8/five updates, with5101 states.
  Directed obstruction/density run20.33s, toy transfer under1s. All work
  serial/remote/bounded/file-free. No construction/paper/library changes,
  commits, pushes, raw data, or new whole-code failure allowance.

Details: [quota audit](coded_state/QUARTET_WINDOW_QUOTA.md).
Next derive a bounded quota-aware envelope or block reduction of the exact
transfer at covered-count111572. Test that it avoids both proved moment
obstructions before widening to all covered and row counts. The full goal
remains active; there is no external blocker and no proof of bad distance.

### 2026-10-30 — Eliminate histories through an exact window-cover formulation

Classify the preceding turn as progress: it proved a first-moment obstruction
and a sharper restricted certificate. This turn changes the proof interface,
not the construction. The full45%/40-bit goal stays active and unclosed.

- For a set P of inactive inputs allowed to be covered, optimize the state
  support exactly as independent allowed runs. The value is a positive
  polynomial of length3/4 windows and length1/2/3 terminal windows. Proved
  max_S F(S)=max_P[H(P)-|P|]. This removes support-history multiplicity
  before taking a moment, while retaining the actual coverage count.
- Every allocation of a window's weight among its inactive blockers gives
  a deterministic upper certificate. On the existing fixed sample, the
  empty-window value is3816.84. Equal splitting gives46354.2; least-load
  allocation gives8545.79; three reassignment passes reach4059.83, near
  the known exact3988.44. All loads and source mass independently replay.
  Sequential assignment dependence is not local merely because windows
  are short. These are fixed-geometry facts, not an ensemble certificate.
- Retained the field constraints before maximizing allowed runs. The
  remaining actual covered count lies in27476..255463. Nonnegative
  multipliers of activity, A-C, and C+D give valid count-specific run
  thresholds. Directed rounding is checked against selected_prior_weights.
  This avoids restricting a geometric maximizer by codeword-only facts.
- Derived a quota-preserving paid-set generating function by independent
  clearing of the original inactive labels. Positive-polynomial comparison
  gives a product-measure upper bound without a conditioning penalty.
  The read4 version remains positive at three discriminating covered counts;
  it is not a completed cover. Hypothetical independent columns give negative
  expressions at selected parameters, but are explicitly NOT valid bounds.
- Derived an exact one-start spatial moment using uniform quartet occupation
  coefficients. At the middle field/coverage parameters, FLOAT evaluation
  gives13.89185 bits versus11.88818 under column independence and26.96341
  under read4. The spatial correction is real. Neighboring start times
  share partitions, so multiplying these moments would be invalid.
- Checks:524288 exact support histories,128 paid-label-set optima,448 thinning
  subsets,256 product-measure masks,315 small-partition occupation cases,
  and256 exact one-start category paths. Optimized and UBSan tests pass.
  Replayed the frozen zero-coverage certificate unchanged at -60.7596261.
  All numerics serial/remote; temporary builds removed; no raw files,
  production/paper edits, commits, pushes, or new full-proof budget charges.

Details: [window-cover proof interface](coded_state/QUARTET_WINDOW_COVER.md).
Next retain short temporal overlap together with the actual partition law,
starting at covered-count111572, or derive a valid allocation-tail bound.
Do not substitute independent columns or a product of overlapping one-start
moments. After a successful bound, cover all covered counts and then the
other unresolved row classes. There is new mathematical and checked numerical
evidence; no external blocker or basis for declaring full closure.

### 2026-10-30 — A proof-method obstruction and a sharper column certificate

The preceding turn left coefficient loss versus witness multiplicity
unresolved. This turn resolves that question and improves a separate
restricted bound. Construction and the full45%/40-bit goal remain fixed.

- Fixed incoming quartet types but summed every outgoing assignment. A
  directed coefficient subcount has regional log2 expectation in
  [43.1050128316,43.1091433871], with A4096,C432,D704 and the exact2100-label
  quota. Earlier narrow integer-type proposals were too restrictive; they
  are preserved as proposals, not treated as full-count estimates.
- Added valid initial and terminal bridges. The complete4096-update
  subfamily has A1048832,C111572,D178860,score18279.6>18171, empty initial
  state, free final weight512, and no terminal deaths. It passes every
  implemented local filter. Its expected count exceeds2^10134.95. Hence
  no exact evaluation of this filtered geometric first moment can close
  the proposed cutoff. It is a counting obstruction, not a bad-code proof
  or lower bound on the probability of a bad geometry.
- Every counted geometry is full-width-burst-good for its counted row set
  (score<=2968<3702). However C+D290432 exceeds the field cutoff255463.
  Thus the obstruction does not survive as an actual bad-codeword family.
  Do not impose that field cutoff on an unconstrained geometric maximizer.
- Grouped entire vertical coordinate sequences before Holder: each input
  label affects at most4 column functions, rather than17 local windows.
  An eight-state exact hidden-bit recurrence retains temporal overlap.
  Positive-monomial comparison transfers the bound to actual quotas.
  For slope1/4 and death weight.74, the all-size2100-row-set zero-coverage
  offset15877 has directed log2 tail<=-60.7596261. The required residual
  bound for positive-coverage components is2294. It is NOT yet established.
- Derived a sound field-cutoff interface F_eta=F-eta(C+D). Its required
  offset is18171-255463eta. Four bounded checks at b=2 show the tradeoff:
  eta.004 removes the frozen subcount's score violation, but leaves only
  1471.148 for positive-coverage components after the zero-coverage bound.
  eta.01 leaves232.37; eta.02 fails this decomposition already at C=0.
  This is not a successful whole-code certificate or a parameter sweep.
- Checks:150 exact coefficient cells,401 bridge masks,35 small partitions,
  1092 exact column paths plus3276 generalized paths and12288 run identities.
  Frozen subcounts replay without proposal optimization. All numerical
  work serial/remote/file-free; no new local build or raw data. The frozen
  diagnose_pair_defects.py SHA256 remains7238FEC9...BBAED3. No changes to
  construction, production code, paper, commits, pushes, or full failure budget.

Details: [counting obstruction](coded_state/QUARTET_FILTERED_COUNT_OBSTRUCTION.md)
and [whole-column certificate](coded_state/QUARTET_COLUMN_RUN_BOUND.md).
Next bound collections of positive-coverage components under the original
score or small field penalty, retaining their overlap and actual quotas.
Do not run the full131841-state filtered transfer: its first moment is
provably too large. Full r=2100 closure and all other unresolved row classes
remain required. There is progress, no external blocker, and no justification
to mark the goal complete or redesign the construction yet.

### 2026-10-30 — Joint marked-region test remains positive

The preceding turn made progress by deriving the hole-birth gate. This turn
implements its joint regional test and establishes a stronger local deletion
rule. These are bounded diagnostic/proof advances, not distance closure.

- Implemented directed coefficient bounds retaining state weight256,
  inherited births44,outgoing deaths44,and the exact2100-label regional
  quota. The birth gate is inside the same polynomial; no marginal
  probability is multiplied into the old kernel. Its tested upper bounds
  remain positive. The complete131841-state transfer was not launched.
- Proved k>=lambda*e+tau*f for deletion of the new holes in a quartet.
  Here k is the following uncovered-inactive count and f counts new holes
  that disappear at the next update. Together with the prior canonical
  filters this leaves153 local terms. Directed regional log2 bounds remain
  1064.228 and2205.227 at log2(b)=4 and8.
- Coverage-prefix bounds also remain positive:1066.425,819.771,640.374
  for caps530,472,432. The last probes a genuinely high-score repeated bulk
  profile; its initial/terminal completion is not asserted. These are upper
  bounds on expected counts, not evidence of bad codes or proven large
  first moments. Need a tighter coefficient or exact-subcount audit next.
- Added a bounded exact row-set best response to the min-cut diagnostic.
  On the existing slope0.25 sample, the original row set is already the
  best response to its support. No rows change, so iteration stops after
  one solve. This is a local check, not worst-row-set or ensemble coverage.
- Verified82944 colored local cases and four independent two-quartet
  normalized calculations over35 partitions. Replayed the exact optimal
  triple checks including the stronger birth budget and existing flow
  verification. All runs serial/remote, no raw files or persistent builds.
  No construction/paper/library/commit/push changes or new failure allowance.

Details: [joint marked-region audit](coded_state/QUARTET_MARKED_TRANSITION.md).
Next: inspect a fixed-type subcount at the432 coverage cap, or sharpen that
coefficient, before enlarging the transfer. Determine whether the remaining
gap is actual history overcount or coefficient-envelope loss. Keep the full
45%/40-bit goal active; all unresolved row classes remain required.

### 2026-10-30 — Component audit and a necessary three-support gate

The preceding turn made progress by proving a harmless subcount and a
component-erasure reduction. This turn tests that reduction and proves a
new constraint on maximizing geometric witnesses. Construction unchanged.

- Full-size fixed-geometry audit: at slope0.25 the3988.44 optimum has30324
  positive components. Only827 cover inactive inputs, with total score447.30;
  all but one component have at most64 vertices. The largest spans10 updates.
  These remain ONE-row-set/ONE-geometry facts, not ensemble evidence.
- At slope0.26 the7832.88 optimum is near the all-selected baseline7659.52.
  The inclusion-minimal complementary hole set forms one component with
  A_H16266,U_H4589,B_H252, span40, and gain173.36. Both extremal tied cuts
  have independently checked objective and feasible-flow certificates.
- Proved a selective zero-coverage run bound using a read17 Holder inequality
  and a positive-monomial comparison to the exact regional quotas. At
  slope0.26,B25604,log2(b)=5/8, the ALL-size2100-row-set tail is at most
  2^(-1515.05). Slope0.25/B18171 does not close by this bound. C=0 was
  already excluded by the field/activity argument; this validates the
  counting method rather than adding new whole-distance coverage.
- Proved a new gate for EVERY geometric maximizer: each new noninitial
  hole must parent an uncovered inactive input at the next update. Hence
  its next quartet is entirely holes and has an inactive hole child.
  No terminal hole births are possible. This is a three-support constraint,
  not an assumption about arbitrary codeword trajectories.
- Defined the exact marked local kernel retaining preceding hole births.
  Independent enumeration verifies265 terms over20736 colored-mask cases.
  The state is(u,b), initially(0,0), with final b=0 and free final u.
  The complete131841-state transfer has NOT been run. Actual quota and
  outgoing deaths must be retained jointly; marginal gate probabilities
  cannot be multiplied into the old kernel.
- Exact gate probabilities at(u,b)=(256,43),(256,44),(256,64) have log2
  upper bounds -98.996,-100.914,-137.263, before requiring the inactive
  child. These are diagnostic marginal gates, not full trajectory bounds.
- Tests:8192 free-run identities,45 colored-partition checks,65536 hole
  score identities,16 exact optimal toy triples including2 actual births,
  30631 full-size birth-bearing quartet gates, and the existing flow tests.
  All numerical jobs were serial/remote. Temporary builds removed; no raw
  files, paper/library changes, commits, or pushes. No new failure-budget
  allowance added. Frozen helper preserved.

Details: [components and the marked kernel](coded_state/QUARTET_COMPONENTS.md).
Next: test that kernel on the failing partial-state family with its actual
quota and inherited births, before any large transfer. This is the next
discriminating test, not another scalar-marker sweep. Full r=2100 and the
other unresolved row counts remain open; the full goal remains active.

### 2026-10-27 — Exact partial-state audit and harmless-history overcount

The previous turn proved the full-width burst certificate and a deterministic
partial-state extension. This turn makes diagnostic and structural progress,
not full closure. The construction and active45%/40-bit goal stay fixed.

- Directed isolated kernel intervals show that the chord loses up to131.46
  log2 bits on tested cells. At256/256 there is essentially no chord loss;
  its large entry persists. Quota-aware16-update bounds remain positive,
  even after reducing log2(b) from8 to4. No broad scalar sweep was run.
- Proved an exact harmless subcount with A4096,C532,D700 and2100 inactive
  labels per region. Its averaged-reference excess is0 and cut score-26.
  Frozen integer histograms and disjoint histogram moves give an expected
  regional continuation count with log2 at least249.25769045. This is not
  a bad-codeword count or a probability.
- Extended that subcount from the empty initial state to all4096 updates,
  with free terminal weight256. Its cut score is-6566 and its expected
  count exceeds2^63840. Thus the unpruned log2(b)=4 moment cannot close,
  even for one row set. These histories also exceed the surviving
  bad-codeword turnover cutoff. The result does not exclude other markers
  or proof methods; it identifies concrete overcounting in this expression.
- Proved additivity over selected interaction components. An inclusion-
  minimal maximum-score support has only positive-score components; the
  minimal maximizer is unique by supermodularity. This provides a selective
  target, but its ensemble tail remains open. Small positive components
  may accumulate; component probabilities are not asserted independent.
- Checks:49 directed coefficient cells,two root inequalities,4096 local
  masks,65536 component identities,2048 supermodularity inequalities, and
  exact replay of all frozen counts/moves. All runs serial,remote,file-free.
  Removed a missing optional SciPy dependency rather than installing it.
  Frozen prior helper SHA256 unchanged. No paper/library/commit/push changes.

Details: [partial-state audit](coded_state/QUARTET_TRANSITION_AUDIT.md).
Next: derive a quota-aware positive-component or constructive-flow deficit
bound that avoids this harmless subcount. Test that reduction before any
large numerical search. Full r=2100 closure and the other row counts remain
open; no new failure allowance has been added to the existing budget.

### 2026-10-27 — Whole-burst routing certificate and partial-state decomposition

The previous turn made progress: it proved a sufficient weighted-routing
condition at r=2100, while its first-moment ensemble test failed. This turn
keeps the full goal and construction unchanged and treats whole runs.

- Proved a directed ensemble routing certificate for all full-width/empty
  histories and ALL inactive-row sets of size2100. Exact complete-region
  hypergeometric moments retain the quota through16 updates. A two-state
  matrix sums256 regions with initial zero and free terminal. At slope1/4,
  death weight37/50, and frozen log2(b)=125/64, the offset is3702 with
  log2 tail upper -66035138609487/1099511627776<-60. Quartet partitions can
  be arbitrary. This is a restricted-family certificate, not whole distance.
- Proved a deterministic extension to partially active states. Relative to
  any full/empty template, coverage loss is at most holes plus4/5 of run-start
  holes plus1/5 of nonterminal run-end holes. Retaining excess deaths gives
  an exact two-state minimum correction. A correction<=14469 meets the
  sufficient routing threshold18171 on the burst good event, with no extra
  union over templates or deviations. Partial-family claims remain over J.
- Checked1048576 support/template inequalities and65536 template-DP cases
  in C++, plus1000 Python candidates,64 exact hypergeometric cases, and60
  independently enumerated regional entries. Directed replay is file-free.
- A rounded illustrative failed-moment trace has A=1537024, minimum D=67712,
  and minimum correction42884.96. It lies outside the new class. This is
  not a codeword or sampled state distribution. The residual issue includes
  partially active transition intervals, not pure full-width switching.
- Derived a second exact decomposition by averaging all512 support-weight
  thresholds: F(S)<=3702+sum_t[L_t max(u_t,v_t)/512-C_t]
  -(37/50)[D-sum_t max(u_t-v_t,0)]. It retains actual routing jointly with
  the reference and removes artificial hard-threshold discontinuities.
- The bounded general screens still fail. A hard-threshold reference gives
  offset1477389.99 and oscillates around the threshold. The averaged
  reference with three log-convex chord anchors gives329810.56. Neither is
  a certificate or evidence of bad code distance. Exact tests cover50 label
  tilts,625 reference paths,625 averaged-threshold paths,100 convexity cells.
- Added the burst exception to the ledger; all proved allowances still sum
  to less than2^-43. Independently replayed the r=2100 implication at B18171.
  Preserved frozen sources and earlier certificates. No construction, paper,
  library, commit, or push changes. All numerical jobs were serial and remote.
  One DNS failure was confirmed terminal and retried after connectivity
  recovered; a command-line size failure was resolved with a minimal source
  bundle. Temporary builds were removed; no raw artifacts were retained.

Next: analyze partial-state coverage and excess turnover using the exact
averaged-reference identity. First compare exact partial-state moment cells
with the chord envelope, to distinguish envelope loss from the combinatorial
count. Do not repeat the broad scalar sweep. The full r=2100 class and the
remaining row counts are still open; the global45%/40-bit goal stays active.

### 2026-10-27 — A sufficient weighted-routing condition for the hard class

The previous turn made proof progress by extending activity exclusions and
introducing a trajectory-sensitive prior. This turn retains the construction
and establishes a direct conditional closure argument, not full closure.

- Proved a joint normalized prior on support histories AND selected zero
  outputs. Its shared Dirichlet transition rows give one global 2^-74
  exceptional-event allowance, with simultaneous parameter choices. It
  covers all row sets, histories, selected zero sets, and free terminal states.
- Eliminated the required zero count to obtain a directed linear necessary
  condition involving A,C,D,w. Combined with A>=743225 at r=2100, it proves
  that C+0.74D>=lambda*A-B suffices for that distance class. Allowed integer
  B values are10739,18171,25604 at lambda0.24,0.25,0.26, respectively.
  This is an exact implication; the ensemble cut hypotheses are unproved.
- Generalized the min-cut and constructive-flow diagnostics to separate
  death and coverage costs. On the existing graph and ONE fixed row set,
  exact offsets are3301.30,3988.44,7832.88 at those slopes. Constructive
  offsets9918.38,14971.01,23175.07 also fit the sufficient thresholds.
  These are all-history fixed-geometry checks, not ensemble certificates.
- Bounded the attempted ensemble tail with the actual regional quotas and
  complete row-set union. Its FLOAT offset screen187529.67 is much larger
  than25604. Dominant artificial paths alternate nearly full/empty support.
- Proved canonical-cut rules by deleting a single state or all incoming
  states in one quartet. A minimum-cardinality maximizing cut has no
  singly covered inactive input when lambda+tau<=1, and obeys
  f<lambda*a+tau*h for a live quartet. These restrict geometric witnesses,
  NOT actual message histories. Exact toy and full-graph checks passed.
  Adding the rules changes the screen only to187364.56, still inadequate.
- Full field-transfer tests with joint weights also failed to improve the
  best screen. Stopped the local-weight sweep instead of treating these
  positive bounds as evidence of bad code distance.
- Preserved all checkpoints and the frozen helper. Tests include268 prior
  comparisons,1680 zero-count eliminations,16 terminal-marker paths,
  393216 toy cuts,61440 run erasures, and strict/tie canonical cases.
  Jobs were serial and remote; temporary builds were removed. No production,
  paper, construction, commit, or push changes. Exclusion budget remains<2^-43.

Next: analyze whole-run geometry or chronological-flow losses to obtain an
ensemble tail for the sufficient weighted cut, including every row set.
The r=2100 implication now closes algebraically once that hypothesis holds;
other row counts still need matching arguments. Do not claim the routing
lemma finishes the entire goal. Stop one-step retuning unless new evidence
shows it addresses the nearly full/empty trajectories.

### 2026-10-27 — All-row activity bounds and a trajectory-sensitive count

The preceding turn established an activity exclusion at one row count.
This continuation extends that exclusion and retains more of the field count.
The construction and the full 45%/40-bit objective remain unchanged.

- Replayed 127 exact intervals covering every unresolved inactive-row count,
  105 through 8191. Endpoint convexity and piecewise-affine field bounds
  cover the intermediate counts; this is not sampled interpolation.
  The combined geometry and field allowance is below 2^-60. Thresholds
  improve the original activity exclusion at 7712 row counts.
- Proved the all-row-count exclusion A-C_T<=128n-65794 with failure below
  2^-63, conditional on any fixed geometry. At n=6092, remaining histories
  must have A-C_T>=713983. No support-history restriction is imposed.
- Proved a normalized universal Markov prior for counting support histories,
  with a normalized row-set prior. It gives simultaneous linear necessary
  conditions indexed by birth/death parameters, using one 2^-74 allowance.
  These priors count histories; they do not model actual encoder states.
  Exact toy checks cover normalization, pointwise comparisons, and row priors.
- Added the resulting coverage, activity, turnover, and terminal weights
  to the full FLOAT diagnostic. The tested weights did not improve the best
  unconditional screen, +271673.74. Tests assuming the unproved small-offset
  routing cut also stayed positive; that cut is not established as sufficient.
- Checked that all proved exclusion allowances, including the stronger
  single-row checkpoint, sum to less than 2^-43. A remaining distance bound
  of 2^-41 would suffice, but that bound is NOT proved.
- Preserved the frozen helper hash and all earlier checkpoints. Numerical
  jobs were serial, remote, and file-free. No production or paper changes.

Progress is structural, not whole-code closure. The remaining obstruction
is the joint count of histories surviving the activity/coverage/turnover
exclusions while producing low-weight output. Scalar weights still let
artificial dense or alternating histories dominate. Next adapt the universal
trajectory prior to those histories and enforce the high-turnover cap jointly
before investing in another ensemble routing lemma. Use a bounded diagnostic,
not an unrestricted parameter sweep; preserve the existing construction.

### 2026-10-27 — A geometry-backed activity certificate

Previous turn was progress: its conditional tests ruled out treating the
fixed-setup cut inequality as the final missing lemma. This continuation
keeps the construction fixed and adds a directed structural exclusion.

- Diagnosed both new dominant bridges368->111->12 and12->120->368.
  Their artificial inactive counts3212 and1610 miss the regional quota2100.
  Their rescaled activity is also below the existing rigorous cutoff.
- Added a three-quarter-class operator with incoming maxima, a common
  complete-region p choice per path class, and free-terminal maxima.
  Twelve rational toy cells validate it. Bins8 improve+303966 to+297813.
  Targeted proof-p refinement and the existing zero-selection-free
  interpolation improve the full midpoint screen to+271673.74. This is
  still positive, FLOAT, and only r=2100; no distance closure claimed.
- Proved A<=Z_T+5C_T, where Z_T counts state coordinates touching no
  inactive input position. Its maximal zero-coverage support is explicit;
  the same fixed diagnostic geometry gives Z_T=475819. Checked196608 toy
  supports and zero-coverage maximality without changing field maps.
- Proved a route-setup tail using parity-separated adjacent-time pairs and
  Cauchy--Schwarz. Exact regional quota conditioning and the row-set union
  give Pr[any T of size2100 has Z_T>=599529]<=2^-60. Replay uses rational
  p=29/100,v=39/32 and directed logs, not numerical optimization output.
- Combined that tail with outer shortening and all-support entropy counting.
  Every nonzero J-message with6092 active rows and A<=747590 is excluded
  with total failure at most2^-59, including every row set, arbitrary
  histories, resets, and the free terminal state. Old cutoff was694488.
  Other row counts and larger-field restricted activities are NOT covered.
- Preserved the frozen helper hash and all previous certificates. Jobs were
  serial and remote; no local data, production/paper edits, commits, or pushes.

This is progress, not completion. Next extend the geometry-backed activity
bound across row-count intervals and combine its coverage relation with
the remaining output count. Keep the global45%/40-bit objective unchanged.

### 2026-10-27 — Exact support cuts and a conditional-closure diagnostic

The construction is unchanged; the whole-code45%/40-bit goal remains open.
This continuation made structural progress, not certificate closure.

- Added QUARTET_SUPPORT_CUT.md and a rowwise rank obstruction: with allowed
  state supports, A>F(S) and F(S)+D<W suffice for a low-output codeword,
  where F(S)=sum_i min(128,h_i). No such witness was found in these tests.
- Reduced min(C_T+D-lambda A) to an exact integer min-cut. On one full-size
  geometry and one fixed inactive set r=2100, verified offsets712.22,
  1063.64,7826.52 at slopes0.22,0.24,0.26 for EVERY support history.
  These are not ensemble or all-row-set certificates.
- Proved a short-run erasure lemma for minimizing cuts. The slope0.24
  optimum has no nonterminal run shorter than5, as independently checked.
- Added a chronological feasible-flow construction. Its exact offset is
  3526.70 at slope0.24 on the same geometry; work conservation is checked.
  A setup tail for this offset has not been proved.
- Marked bridge calculations expose inactive-label quota mismatch in the
  old sparse/dense artificial paths. Reoptimizing two proof parameters
  improves the full midpoint FLOAT screen to+303965.75, still positive.
- Tested the cut inequality's utility BEFORE pursuing an ensemble proof.
  Assuming R>=0.24A-5000 and applying marker2 gives+394432.61, worse.
  Activity marker1.5 without that assumption gives+462968.00. These tests
  do not show impossibility, but do not justify calling expansion the last
  missing lemma. The ensemble cut-moment screens themselves remain loose.
- Exact tests cover196608 toy cuts,36864 run erasures, three optimized
  flows, three conditional-cut inequalities, and the previous routing tests.
  All jobs were serial on Peach; temporary build files were removed.
  No local numeric data, production/paper changes, commits, or pushes.

Next: identify the remaining histories under the joint distance/expansion
count and their rowwise constraints. First demonstrate that a proposed
geometry lemma makes that count useful; only then invest in its ensemble
tail. Preserve all restricted certificates and the fixed construction.

### 2026-10-27 — Midpoint classes before the routing minimum

Previous turn classified as progress: it added an exact coverage/death
exclusion and identified a remaining-count window. This continuation keeps
the construction fixed and works on the full513-state routing calculation.

- Restored the zero-selection-free interpolation to the all-weight kernel.
  The bounded tau0.1 test is worse (+776562.53); no closure claimed.
- Inspected the actual artificial flow of the full regional envelope,
  rather than importing the old constant-p diagnostic. The endpoint-only
  envelope concentrates near weight272, not489, at regional boundaries.
- Derived a sharper operator: retain the midpoint weight w and sum over w
  only after minimizing the complete-region Bernoulli-p quota bound. Both
  halves use the same p before this minimum. No setup change or resampling.
  The FLOAT screen improves +572899.41 -> +404889.64 at r=2100.
- Retained complementary half-region count intervals, then quarter-boundary
  weights inside those intervals. Width64/16 half-quota bounds do not improve
  the full result. The width256 quarter refinement improves some entries,
  but not the dominant full-path bound. The joint coverage/death marker0.9
  is also worse (+420216.43) on the refined operator.
- The refined artificial boundary distribution puts about49.61% below128
  and50.39% above384, with17.85% at zero and a high peak near489. This is a
  distribution of the proof envelope, not the encoder or a bad codeword.
  The edge flow is almost entirely low-to-high then high-to-low. Traced
  representative edges0->489 and489->0: their largest midpoint classes are
  348 and332, with gradual growth/decay rather than abrupt reset. This
  locates a concrete hard family for further joint counting.
  Keeping all three quarter weights on those two edges saves only88.84 and
  73.49bits; it is not a complete-operator replacement. The remaining
  global gap is about1582bits per region, so finer checkpoints alone are
  not the most convincing next direction.
- Added QUARTET_PATH_ROUTING.md with explicit operator/conditioning proofs.
  CWC guided the distinction between proof choices and actual setup.
  Tests check27 exact midpoint components,54 exact split-quota components,
  54 exact quarter components,27 selected-edge components,4 exact bridge
  marginals,16 score telescopes, and all preceding tests.
  Retained free terminal states, every turnover, and all513 state weights.
- Jobs remained serial, remote, and file-free. Enlarged the explicit helper
  guard to120s/128MiB for retained quarter matrices; the longest run took74s.
  No local numeric job, disk cleanup, native benchmark, construction change,
  library/paper edit, commit, or push. A diagnostic variable shadowing an
  imported module was caught before numerics and corrected.

This turn is progress in the unrestricted routing bound, not a distance
certificate. Next resolve the remaining sparse/dense path family jointly
with outer constraints. Every row count and a directed final union remain
required; the original45%/40-bit goal stays active.

### 2026-10-27 — Exact coverage/death exclusion and remaining-count tests

Previous turn classified as progress: a uniform low-activity lemma removed
one global family. This continuation keeps the construction fixed.

- Proved a second all-row-count exclusion using R=C_T+D, the covered
  inactive input positions plus support deaths. Low-output histories with
  R>=max(0,W-128n+91522) have total setup probability below2^-63.
  The support/selected-zero automaton is [[1,1],[1,2]]; a rational potential
  bounds its count by (21/8)^N. The final union uses integers/rationals only.
- This covers every low-output history with n>=8088 active rows. It does
  not close the other classes. At r=2100, the remaining event has A>=694489
  and R<=255463. The combined exclusions cost less than2^-43 globally.
- The exact inactive quota gives C_T>=256r-U, with U uncovered inputs.
  Thus U-D>=282137 at r=2100, and consequently A+D<=1815015. This excludes
  some dense artificial histories without changing the construction.
- Added exact automaton replay, 5460 quota checks, nine short-path marker
  checks, and54 activation/reset cells. The full helper accepts independent
  death and inactive-coverage markers; defaults preserve prior behavior.
- Bounded full-state FLOAT tests of the remaining constraints did not
  improve the baseline +572899.41. Joint coverage/death marker0.9 gives
  +592233.64; free-minus-deaths marker1.1 gives +622120.49, and marker4
  gives +3227199.07. Activity tilts and a theta0.5 comparison are also worse.
  These positive upper bounds are not codeword counterexamples.
- Added and indexed QUARTET_HIGH_TURNOVER.md, quartet_high_turnover.py,
  and quartet_remaining_activity.py. CWC guided the proof's quantifiers and
  scope. Numerics stayed serial and file-free on Peach,90s/64MiB bounded.
  One lost job output was replayed only after confirming no remote job
  remained. No disk cleanup, setup, library, paper, benchmark, commit, or
  push changes. The frozen pair-defect helper hash is unchanged.

This turn is progress through an exact new exclusion, not closure. Next
retain the joint remaining counts inside routing, rather than continue
global-marker tuning. The all-history, all-row-count directed union and the
full45%/40-bit goal remain open and active.

### 2026-10-27 — Uniform low-activity certificate and joint quota diagnostics

Previous turn classified as progress: it isolated a quota mismatch and
built the full regional matrix. This turn keeps the construction fixed.

- Added count-bin regional composition, covering all weights and turnovers.
  It excludes the persistently low-label paths, but coarse interval costs
  leave positive screens: +785175 (state/count widths8/16), +708280 (1/16),
  and +647301 (2/8). These do not improve the previous best global screen.
  Thirty exact toy components check the interval and quota inequalities.
- Added exact three-axis coefficients for near-boundary fixed profiles.
  Weight489, arbitrary turnover, r=2100, first activation1 gives FLOAT
  -377912.25. Matched optimized Cauchy is already -376199.63, so the exact
  quota saves1712.62bits, not a new whole-code closure. Eight diagonal and
  eight artificial-p invariance checks pass.
- Proved a different, global-in-row-count low-activity lemma. The conditional
  pin/shortening count without selected output zeros is <=2*q^(A-128n).
  Summing all row sets and all supports with A<=b(n)n gives failure<2^-44.
  The dyadic thresholds b(n) rise from84 to118. This covers arbitrary
  varying states, resets, first activations, and free terminal states overJ.
  At r=2100 the excluded range is A<=694488. It is not an output-distance
  statement and restricted activity does not automatically lift to extensions.
- Added QUARTET_LOW_ACTIVITY.md and quartet_low_activity.py. The scalar
  replay uses exact rationals/integers, including all fourteen row-count
  ranges and the final union. No FLOAT evaluation enters this certificate.
- Preserved all checkpoints and indexed the result. Used CWC for the proof.
  Numerical work stayed serial and file-free on Peach; local disk remains
  below20GiB. Largest screens took63s under a120s watchdog and128MiB numeric
  budget. No setup, production, paper, benchmark, commit, or push changes.

Next: impose the certified A>b(n)n condition in the remaining full distance
count, rather than continue paying for excluded low-activity histories.
The low-activity event uses one2^-44 allocation across all row counts.
Remaining low-output histories and their directed union are still open;
the full45%/40-bit goal remains active.

### 2026-10-27 — Actual free positions and full regional matrices

Previous goal turn classified as progress: it built an all-weight kernel
and separated global overcounting from the restricted certificates.
Continued on the full objective with the construction fixed.

- Derived a regional cap using genuinely uncovered positions, with an
  affine marker envelope. Added quartet_covered_route.py and optional
  marked derivatives to the new full-weight helper. The all-weight FLOAT
  screen improves to about +622276, still uninformative as a probability.
- The marked artificial flow has mean state weight 488.55, almost no free
  positions, and only 20.63 inactive covered positions per update. A real
  region averages 131.25 inactive labels. This directly identifies loss
  from maximizing the missing quota independently of realized occupancy.
- Added quartet_regional_matrix.py: all 513 endpoint weights, no bins,
  complete 16-update matrices before entrywise minima over proof p.
  A nine-p screen gives +572899 at r=2100. Bounded theta=0.5 and ell=0.001
  comparisons are worse. These are FLOAT screens, not certificates.
- Checks: 275 exact rational cap comparisons; empty-quartet and marked
  identities; 567 direct matrix-product checks, including 110 underflow
  fallback cases; full-region vector comparisons within 3e-11 bits.
- Local disk crossed below20GiB. The guard stopped the next local run
  before computation. Used configured Peach through the connection skill,
  executing compressed source in memory, no remote files or caches. Runs
  stayed serial and bounded; no deletion, benchmark, production, paper,
  commit, push, or setup changes. CWC distinguishes inequalities from screens
  in coded_state/QUARTET_COVERED_ROUTING.md.

Next: enforce the actual inactive-label total jointly with the regional
state path, using coefficient extraction or a residual-quota bound. The
dense artificial flow supplies a concrete diagnostic that the stronger
bound must exclude. Rowwise RS feasibility, all-r coverage, and directed
closure remain open. Goal remains active; no external blocker to proof work.

### 2026-10-27 — Full-weight transfer and outer-coverage audit

Continued the active whole-code goal without changing the construction.
The preceding turn's joint-occupancy work was progress, not closure.

- Added quartet_weight_flow.py. Retaining all 17 live weights improves the
  matched restricted-band screen by about 5037 bits at r=2100, but leaves
  +7226.80 at the proposed tau=0.092430. Other tested slices are -1083.72
  and +1067.31. FLOAT only; not coverage of intermediate row counts.
- Added quartet_full_weight_kernel.py. A two-axis coefficient calculation
  includes all 513 incoming/outgoing weights and all deaths. Each build
  takes roughly 6--7 seconds. Independent activation/reset formulas and
  289 comparisons against the restricted three-axis kernel pass.
- Added a necessary RS coverage marker: every active row needs at least
  129 covered positions. A repeated 0->weight-one->0 history covers at
  most 10240 positions, versus 785868 required at r=2100. This rules out
  a concrete class of spurious histories; it is not a code obstruction.
- Full regional screens remain strongly positive. At ell=0.00377 the best
  six-iteration all-weight bound is about +685155, not near closure.
  Coverage marker 1.5 does not improve it. Constant-weight diagnostics
  identify a separate dense-state mismatch in the inherited proof-p grid.
  Adding p=0.94 improves that diagnostic but not the full regional bound.
- Tested 33-bin regional endpoint composition. Original and constant
  potentials give +827219 and +1042599 respectively. Park this coarse
  reduction; it is not a numerical certificate or an improvement.
- Added independent checks: 9 exact coefficient cells, 117 explicit paths,
  81 potential-conjugacy cells, and 81 coverage-monotonicity cells.
  Recorded scope, inequalities, failed screens, and commands in
  coded_state/QUARTET_FULL_WEIGHT_FLOW.md using CWC.
- All runs serial and source-only, with 60-second watchdog and >=20GiB free.
  The initial eight-p weight-flow run hit its watchdog; the matched run
  uses four p values. No generated dumps, benchmark, production, setup,
  paper, commit, or push changes. Frozen diagnostic hash unchanged.

Next: retain actual regional coverage and inactive occupancy jointly with
weight flow, rather than pricing every coordinate of a nonquiet update as
constrained. Recover rowwise RS feasibility, then extend to all r and
directed arithmetic. Restricted band closure remains useful but does not
show that the whole-code proof is close. Goal remains active.

### 2026-10-27 — Joint occupancy, classwise proof weights, and mixed minors

Previous goal turn classified as progress: it introduced the regional bound
and identified the separate-occupancy relaxation. Continued on the same
whole-code objective, without construction changes or new subagents.

- Retained total nonquiet count before selecting tau and summing histories.
  FLOAT bounds at r=1850,1974,2100 are approximately -17567,-2867,+14866.
  Extending regional p candidates up to0.76 did not change those results.
  No intermediate-row coverage or directed certificate is claimed.
- Implemented joint incoming-weight/inactive-occupancy coefficients and an
  exact regional-quota extraction in a two-state polynomial matrix. Baseline
  FLOAT runs take8--10s. Outgoing weight/death interval markers improve the
  bound modestly but do not close it; r=2100 is still about+21348 at tau0.1.
  These direct sums differ from the best-per-nonquiet-count calculation.
  A complete artificial-p invariance test agrees within4.66e-10bits.
- Identified that tau interpolation is dominated by using the smaller
  output marker ell'=ell*(q*ell)^(-tau) directly. Preserve the previous
  calculation, but do not treat tau as indispensable new machinery.
- Proved an actual-setup mixed-minor lemma:74 extra minor polynomials have
  total degree224. More than6 exceptional quartets has probability below
  2^-48.64778. Their designated-location union costs at most104.63079bits
  at the retained markers. All4096 local triples,5120 polynomial inequalities,
  and74 actual-field witness minors pass; directed exponents are frozen.
  The stronger local count improves this screen by less than1bit, so it is
  not worth pursuing here. No rejection rule was added.
- Added QUARTET_JOINT_REGIONS.md, QUARTET_MIXED_MINORS.md, and their source
  helpers. The CWC skill keeps the proved lemma separate from FLOAT screens.
  Runs remain serial, bounded, source-only, >=20GiB free; no production,
  paper, benchmark, commit, or push change.

Next: retain state-weight flow across updates/regions. The current scalar
band envelope can take different incoming-weight maxima at successive steps.
Test that relaxation before claiming a code obstruction or doing more local
minor refinements. Whole-code coverage and directed closure remain open;
the goal stays active.

### 2026-10-27 — Active whole-code goal; regional long-quiet diagnostic

User explicitly requested a goal to close. Registered the objective: the
unchanged quartet code, K=2^20, rate 1/2, 45% distance, at least 40 bits,
all messages, row counts, histories, and free terminal states. Keep it active
until the full coverage and directed union close. Restricted families are
progress checkpoints, not completion conditions.

Added coded_state/quartet_quiet_regions.py and QUARTET_QUIET_REGIONS.md.
No frozen helpers or construction parameters changed.

- Derived a regional cap by convolving quiet-position Bernoulli labels
  with the nonquiet labels, then maximizing the remaining binomial mass.
  Wholly quiet regions cancel their quota cost exactly.
- Combined the selected-zero and zero-selection-free conditional counts
  before averaging. Their interpolation retains the actual A+D score and
  gives a composable transfer without a cap on zero-run length.
- The cheap FLOAT screen at r=1850 improves from +23367 to -11255 bits.
  At r=1974 and 2100 it remains positive. This is not a finite certificate
  and does not cover untested row counts or other live-weight bands.
- A two-box comparison retains outgoing weights at the largest p values.
  It improves r=1850 to about -17093, but leaves the other gaps. Outgoing
  weight relaxation is not the sole loss. At tau=0.10 the three screens
  are approximately +8011, +10163, +16204.
- The remaining regional cap maximizes inactive-label occupancy separately
  from live transitions. Next retain that occupancy jointly, then convolve
  with the quiet binomial and impose each exact quota before composition.
- Passed 81 exact conditioning checks, 21 exact interpolation checks, and
  56 short-path coefficient checks. Replayed all three prior band/exception/
  reset bounds unchanged; the frozen pair-defect hash is unchanged too.
  All results are source-only, serial, with
  the 20-GiB guard. The initial eight-box comparison hit its 60-second guard;
  the bounded two-box comparison completed in about 31.5 seconds. No raw
  dumps, benchmark, production edit, commit, or push.

Whole-code closure remains open. Preserve the earlier band/exception/reset
certificates and continue on the joint regional count, not a mixer redesign.

### 2026-10-26 — Band exits and arbitrarily many brief reset episodes

Added coded_state/quartet_excursions.py and QUARTET_EXCURSIONS.md. The
construction, acceptance rule, and frozen band markers are unchanged.

- Derived a five-local-type envelope covering every state weight and
  turnover, including zero. Under the existing potential its log row bound
  is53.05571; regular band transitions have log row bound-271.80293.
  This avoids constructing a513-cubed array. The maximum local envelope
  occurs at empty incoming quartets; checked by exact fortieth powers.
- Noncommuting positive-operator bounds sum every placement of up to25
  arbitrary exceptional transitions. Initial/terminal weights are unrestricted.
  At45%, with1850--2100 inactive rows and every first activation time,
  the directed exponent is-3901196831485731/43980465111040 (~-88.70295).
  The26-exception envelope is inconclusive, not a code obstruction.
- Counted actual resets separately with a one-dimensional Bernstein DP.
  All local reset masks satisfy precisely |T|<a. A512-bit positive integer
  calculation bounds the band-to-zero row by2^-834.09978. The reverse
  row is at most2^-311.10508; a zero-to-zero edge has weight one.
- A four-state live/zero-age potential now covers arbitrarily many reset
  episodes, each with at most3 consecutive zero states. All nonzero weights
  stay216--232; nonzero-to-nonzero transitions have24--40 deaths. The full
  family union, including all activation times and free terminal zero ages,
  has exponent-13146227725709835/2199023255552 (~-5978.21223) at45%.
  This is separate from the25-arbitrary-exception family, not a claim that
  both freedoms can be combined without another proof.
- One shared outer message space and row-set union remain throughout.
  Quiet input equations stay in outer shortening; no independent outer
  resampling or restart at reactivation. Both results remain restricted.
- Frozen source-only replay asserts both exponents and the old band result.
  Exact tests:40 noncommuting-path bounds,125 occupancy products,72 Bernstein
  cells,64 reset masks,4 age-potential rows and32 terminal-age bounds.
  Old pair-defect hash unchanged. Serial~11s runs,60s watchdog,>=20GiB free,
  bounded memory, no outputs/builds/benchmarks/production edits/commits/pushes.

Next: longer quiet intervals need the actual regional quota and RS coverage
constraints, not automatic inheritance of live-band decay. Other nonzero
weights and row counts also remain. Short reset mechanisms are now covered
within this band; do not redesign the inner solely because a quiet loop
breaks a scalar continuation envelope.

### 2026-10-26 — Compose varying histories and every activation time at45%

Continued toward end-to-end coverage with the construction fixed. Added
coded_state/quartet_history_band.py and QUARTET_HISTORY_BAND.md.

- A positive transfer and telescoping potential sum all histories with
  live weights216--232 and24--40 deaths per transition. Births vary with
  weight changes; no fixed schedule or independent-state assumption.
- Exact positive coefficients are retained in a233-by57-by41 box, below
  16MiB including temporaries. Upward integer rounding, overflow guards,
  and rational logarithms give directed bounds without an FFT.
- Cover all inactive-row counts1850--2100 using binomial log-concavity,
  exact endpoint bounds, and an explicit251-count union. All row sets and
  messages are included in the conditional count, as in the prior proof.
- The original markers close first-step activation at45% below2^-33545,
  and activation in the first128 updates below2^-15639. At46%, the broader
  first-step band remains inconclusive (+72650); no counterexample claimed.
- Complete zero-prefix regions cancel their label-conditioning cost exactly.
  RS injectivity rules out first activation after step2048. A proof-marker
  proposal plus directed replay closes every possible activation time in
  the band, including the free terminal state. Frozen log2 upper exponent:
  -106287443107715149/10995116277760, about-9666.78664. No quota atom cap is
  needed for this final result. This is restricted coverage, not whole-code.
- A zero-selection-free fallback separately closes selected late-start bands.
  It helped distinguish a proof-weight issue from a construction obstruction;
  the final frozen bound uses one marker set for all2048 activation times.
- New exact checks include1536 labeled partition cases,9 varying-weight
  identities,16 DP cells,2 full-path potential inequalities,406 binomial
  inequalities,7 log checks,5 prefix quota identities,9 negative-zero checks.
  The prior224/32 frozen46% replay remains unchanged at-104.10742.

Final replay took about8.5seconds. Serial bounded file-free runs, >=20GiB
free, no benchmark/build/raw dump, no production/paper change, commit or push.
The frozen replay uses stored rational weights, not optimizer convergence.

Remaining: other row counts, state weights/turnover outside the band, exits
and returns, and full resets/reactivation. Next extend the composable transfer
to these cases while retaining the outer cost of quiet intervals. Do not
claim that strong margins inside this band cover paths leaving it.

### 2026-10-26 — Close the 224/32 gap with exact counts and regional quotas

User requested bounded two-update and rowwise analysis with the construction
fixed. Added coded_state/quartet_two_update.py and QUARTET_224_TWO_UPDATE.md.

- Proved support-only two-update composition equals the square of the exact
  one-step restricted row sum. Shared support labels alone do not improve
  that relaxed count. No code change is justified by that observation.
- Replaced three transition coefficient evaluations by positive integer
  coefficient extraction at (incoming,births,deaths)=(224,32,32). With the
  frozen rational theta=7/40, the directed exponent improves from +53740.2904
  to +1057.38794. The 35-term DP uses bounded uint64 arrays and explicit
  overflow guards; no FFT or floating-point result supplies the bound.
- Retained the exact count of 1974 inactive rows in each of the 256 regions.
  Conditional local label-count laws admit a certified rational dual bound
  on adjacent probability overlaps. A Fourier/Gaussian atom bound recovers
  another 1161.49536 bits while retaining the original conditioning cost.
- Final directed exponent: -4578692708060433/43980465111040, about -104.10742.
  This covers all GF32 messages/support histories in the fixed 224/32,
  1974-inactive-row family at 46%; not other profiles or a whole-code claim.
  No stronger mixed-diagonal acceptance or construction change is used.
- Derived the individual active-row coverage generating function and noted
  the cross-row dependencies. The individual 129-position floors are not
  needed for this closure; no unjustified product of row-tail bounds is used.
- --frozen replays stored markers/dual without optimizers, asserts the exact
  exponent, and checks integer arithmetic, roots, partition identities, and
  Fourier bounds. Previous frozen trajectory replays pass. Serial, file-free,
  roughly four-second runs with 60-second watchdog and >20 GiB free. No timing
  run, production change, raw dump, commit, or push.

Next: bounded neighboring-profile coverage using the same refinements, before
any inner redesign. Reserve must be summed explicitly across covered families.

### 2026-10-26 — Native quartets and varying-state trajectories

Added native/quartet_native.cpp, coded_state/quartet_trajectories.py, and
coded_state/QUARTET_NATIVE_TRAJECTORIES.md. Fixed outer/global route; no
promotion, paper change, or whole-code certificate.

- Pair/quartet arithmetic: 3.31/4.29 ms. Naive gathered/scattered recurrence:
  9.50/12.58 ms. Partition-ordered recurrence: 4.20--4.50/5.54--5.78 ms across
  two runs. Excludes outer/global route and setup. Sixteen accepted keys per
  family are reused in expanded coefficient streams; not a full independently
  sampled code. Scalar/SIMD checks pass. Outer-route layout integration remains.
- The native tower has its own checked common-minor witness and trace-one
  constant; no unsupported basis-change assumption is used.
- Directed 46% bounds close alternating 240/272 and 224/288 schedules and
  selected constant weights. Selected full-reset/high-turnover schedules are
  impossible by support counting, not every reset-containing history.
- New hard family: weight 224, 32 deaths/update, 1974 inactive rows. Directed
  46% log bound +57302.414 (inconclusive); 45% gives -57208.031. Optimizing
  theta at 46% still gives FLOAT +53739.703. An ideal uniform dense inner
  gives restricted exponent -11931328 with the same outer and route.
- Stronger mixed-diagonal acceptance improves the FLOAT gap by only 1.21 bits;
  an aggregate active-row coverage marker gives no gain. Stronger acceptance
  is not adopted. No actual low-distance quartet codeword was established.
- Previous 943 algebra checks, original rational replay, and 10740 pairing
  checks pass. All numerical runs serial/bounded; disk floor respected.
  No raw dumps/local build tree, and frozen helpers remain unchanged.

Next bounded test: two-update constraints with surviving coordinates and
rowwise outer occupancy at the 224/32 family. If that cannot recover the gap,
test modest cross-quartet mixing before another implementation. Do not infer
arbitrary-schedule coverage or resume a broad sweep from these family bounds.

### 2026-10-26 — Compare cheap inner mechanisms before tuning another bound

User requested fixed outer/route, structural obstruction checks, then proofs
for surviving cheap inners. Compared six-shear recursive quartets, concentrated
rank-one refresh, and a field-multiplication realization of the quartet control.
No new production construction, native timing, or whole-code certificate.

- Recursive forward/backward sweeps have tridiagonal inverses, including
  weight-at-most-two preimages of boundary unit outputs. An incoming boundary
  singleton also leaves one output coordinate fixed. Random ordering does
  not turn this into the uniform quartet law. These are local interface
  obstructions, not full-code counterexamples.
- A raw fresh rank-one update I+u v^T has exact fixed-input law consisting
  of a 1/q unchanged atom and a (1-1/q) uniform vector. Nevertheless, every
  update can be neutralized by one homogeneous equation on the message.
  With 240 of 512 fixed state coordinates allowed, the expected shortened
  outer dimension is 5429.210690 versus at most 4095 neutralization equations.
  Exact expectation/cap arithmetic gives probability >0.0027 of a codeword
  of relative weight <=46.875%, for every inner-key choice. Four low-rank
  directions also fail a 49.1%/40-bit target by a >0.0015 lower bound using
  weight <=49.0234375%. This is a real ensemble obstruction. The same test
  does not rule out a rank-one 46% profile.
- Constructed an explicit degree-four extension of GF(2^32). Uniform raw
  field multiplication gives an exactly uniform image for every fixed
  nonzero base-field quartet. An actual-field witness verifies all 138
  square minors of A and A+delta I are nonzero polynomials. Rejection thus
  has total density <=1024/989<2, the same allowance as the quartet control.
  Two nested quadratic multiplication circuits use nine variable base-field
  scalings and fifteen XORs per quartet, with reduction constants prepared
  in setup. This is cheaper than direct sixteen-product matrix arithmetic,
  but 50% more variable scalings than two accepted pair maps; no timing inferred.
- The directed hard-family replay transfers unchanged: log2 bound
  -2842457/64 for half-full histories, 32 swaps per transition, 1974 inactive
  rows, at 46%. All twelve neighboring FLOAT screens again pass, weakest
  about -99575.69. No unrestricted activity/turnover/reset coverage is claimed.

Added INNER_MECHANISMS.md and inner_mechanisms.py. Checks include 943 local
algebra fixtures, 64 exact rank-one law bins, exact rational obstruction
arithmetic, actual GF32 minor witnesses and subfield embedding, circuit/reference
agreement, and a counted nine-product runtime schedule. The polynomial-basis
generator was initially not of the required order; that assertion rejected
the attempt. A verified order-255 search and AES-root embedding now pass.
The final directed replay runs in under four seconds; the twelve-profile
FLOAT screen runs in under seven seconds. All runs are serial, bounded,
file-free, and under the 20-GiB disk guard. No raw data, build, benchmark,
paper change, promotion, commit, or push. Frozen pair-defect helper hash
remains 7238FEC9E5553E7B4D1CE3EDB3089CEC06ABCC34CB090C029D0F3D6B54BBAED3.

Next: broaden the field-quartet proof to varying activity and turnover, and
time this exact circuit against pairs before selecting it. Recursive sweeps
remain open under a different count. Identity-base low-rank refresh needs
a structural countermeasure before pursuing the near-half ladder.

### 2026-10-26 — Is random pairing sufficient, or should the mixer change?

User requested a design assessment, not another coverage sweep. Kept the
current code fixed and compared it analytically with uniformly partitioned
quartets whose matrices and shifted matrices are superregular. The quartet
point bound has total accepted-map density below 1024/989, hence below two.
This is a distinct comparison distribution, not a claim about two pair layers.

Derived a joint capped count retaining exact random-group occupancy and
necessary local zero-input/zero-output constraints. For the original
half-full, 32-swaps-per-transition, 1974-inactive-row family at 46% distance,
the pair FLOAT bound improves from +286860 to +282899 bits. This is still
inconclusive. The quartet counterpart gives -99576 bits. A separate exact
rational replay, with algebraic roots rounded upward, proves the restricted
quartet bound -2842457/64 in log2 units. It includes all messages, support
histories, and row sets in that specified family, not the whole code.

A bounded neighboring screen used swaps {16,32,64} and inactive-row counts
{1024,1974,2560,4096}. All twelve quartet FLOAT bounds are negative; five
pair bounds remain positive. Replacing 25%, 50%, and 75% of updates by
quartets gives +192759, +98523, and +960 for the original profile. These
shared-marker screens do not establish optimality, and do not support a
claim that rare stronger updates suffice.

No counterexample to pair-code distance was found. The recommendation is
nevertheless to reconsider one pair layer as the final target: broader
local mixing gives a substantial improvement under the same proof method,
whereas further constraints inside a pair give little. More coefficient
randomness alone cannot materially improve its existing q^-2 point bound.

Validation: 10740 local/partition/coefficient checks pass. One initial
FLOAT proposal underflowed and was rejected; exact positive recentering
fixed evaluation before subsequent screens. All selected screens then
finished within five seconds, serially and without data files, under the
60-second watchdog and 20-GiB free-space guard. No native timing, production
change, paper edit, commit, or promotion. The uniform quartet comparison
still needs broader proof coverage and an efficient actual realization.

Next: test varying activity and turnover for the quartet control, while
examining a cheap implementation distribution. Do not transfer its proof
to factorized maps without an image-law argument. See the
[sufficiency analysis](coded_state/RANDOM_PAIRING_SUFFICIENCY.md).

### 2026-10-26 — Conditional audit before more flow coverage

User requested reassessment and an audit. Kept the construction fixed and
granted the proposed flow conclusion hypothetically. Derived the general
conditional count for arbitrary state supports, deaths, resets, terminal
state, and exact active outer rows. The pin identity gives the exponent
E=d-J-D-z. When output-zero constraints are needed, it simplifies to
E=-(K-W)+128|T|-(C_T+D).

At |T|=1974 the selected coverage threshold 168785 gives only E<=0.
It does not pay for support histories or selected zero-output positions.
The simultaneous geometric conclusion removes a history union from the
geometry-failure event, not from the subsequent field count.

For full-length half-full histories with exactly 32 swaps per transition,
the exact counts cost between 1109590 and 1109591 bits for histories,
720465 and 720466 bits for zero selections, and 6519 and 6520 bits for target
row sets. The direct conditional union needs 57395 additional field-symbol
units of coverage to give a 2^-41 good-geometry bound. Its sufficient target
is 226180, requiring causal loss<=142320, or 6.788%, rather than 17.738%.
The retained mean-field diagnostic predicts 12.34%; it is not a theorem,
but exposes the mismatch between the numerical goal and this composition.

This does not exhibit low-weight words or invalidate restricted proofs.
The count may include many histories incompatible with the outer and
geometry. A sharper joint count, or a sharper geometric estimate, could
recover the loss. Low-activity regimes and other target-row counts also
remain uncovered. Whole-code distance is still open.

Added a bounded exact audit helper: 2304 geometry identities, 21168 exponent
checks, 169 pin-rank checks, and full-size integer combinatorial counts all
pass. The serial file-free run took 2.254 seconds under a 60-second watchdog,
with more than 20 GiB free and fixed working sizes below 256 MiB. No coverage
search, benchmark, construction change, paper edit, commit, or promotion.
Existing partial coverage and frozen helpers are preserved.

Next: derive a joint weighted count for the 32-swap family, retaining outer
compatibility and matching constraints before taking a geometry cutoff.
Require a plausible end-to-end negative exponent before resuming broad
histogram certification. This supersedes the numerical next step below.
See the [conditional audit](coded_state/RANDOM_PAIRING_CONDITIONAL_AUDIT.md).

### 2026-10-26 — Diagnose the exposure cell before more coverage search

The previous turn only answered a status question: no proof progress.
This continuation revalidated the pending work and found no live numerical
process. The construction, potential, and selected target-row count1974 stay
fixed. The new evidence changes the next proof action, not the code design.

The difficult155,373,250-histogram box splits exactly into three exposure
cells. Phase2 has105,823,530 histograms. Its degree-three bound improves
from9.646791297 with a deliberately forced nested-matrix correction to
9.530488820 with exact phase quotas. This is a small improvement, not the
main remaining loss. A center-based triangulation gives the same bound:
the worst control lies on the boundary h0+h1+h2=255.

The nearby integer histogram(53,88,114,86,171) passes at6.988341870 with
the same Gaussian center and priority. Its80-digit replay agrees; changing
only the center gives6.934139713. This diagnoses slack near the offending
control, not a uniform bound over the box or a whole-code property.

Added exact rational proof-mesh midpoints and scalar degree elevation
after integrating the common histogram controls. The latter preserves
their actual barycentric weights. At degree8 it lowers the phase2 bound
to9.101206642; subdividing two boundary edges gives9.080092587. Both fail
the7.679295 allowance. The worst control shifts to a neighboring triangle.
These degree-three comparisons do not beat the preserved degree-four
parent-box bound9.025866; their value is diagnosing the loss.
A preceding one-edge run ended at its55-second deadline without a bound.
The completed two-edge run used68.031seconds within its95-second budget.

The source checkpoint remains23 regions,645,724,890 certified histograms
and2,273,894,055 pending. No additional coverage is claimed. The original
80-digit point control still returns7.672286056280 and passes5021 checks;
the default coverage geometry suite passes13107. Separate scalar-elevation
tests pass1635 exact checks. Phase and rational-mesh tests are recorded in
the linked diagnostic. A fresh80-digit replay of checkpoint region22 agrees
at7.37675208771730657765, covering its32,314,965 histograms. Earlier checkpoint
regions were not all replayed in this continuation.

All numerical work is serial, bounded, file-free, with more than20GiB free.
The frozen fixed-pair helper hash is unchanged. No native benchmark,
library/paper changes, commit, push, or promotion. Goal remains active.

Next split the full h3/h4 exchange direction of the offending boundary
face, rather than chase individual control coefficients. If that still
fails, isolate the Jensen loss from mixing Gaussian control exponents.
Other histogram regions, target-row counts, low-activity histories, and
field-cancellation composition remain open. See the
[exposure-cell diagnostic](coded_state/RANDOM_PAIRING_FLOW_CONTROL_POINTS.md#exposure-cells-and-the-remaining-interpolation-loss).

### 2026-10-26 — Another intact region and an endpoint-aware flow interface

Previous turn made progress through elevated controls and bounded tails.
The next intact region ([0,79],[80,159],[234,287],[0,198],[0,198]), sum512,
passes at degree four with upper6.818549 and at degree three with7.376753.
Retained the cheaper degree-three certificate:275 controls,32,314,965
histograms. Its80-digit replay with13107 checks takes32.747s and returns
7.37675208771730657765. Checkpoint now23 regions,645,724,890 certified and
2,273,894,055 pending in four regions. Prior22 certificates were not
replayed here; their verifier/source paths are unchanged.

The adjacent155,373,250-histogram region still fails at9.025866, with
essentially all of that bound in the central contribution. Two degree-four
axis-split probes each hit their43-second deadline. Both remain pending;
timeouts are not failed points or construction obstructions. An initial
probe invoked split_region from the wrong module and terminated before
numerics; the corrected serial run completed with the two timeouts.

Also broadened the deterministic flow interface. For arbitrary code-state
supports H_1,...,H_L, set S=sum|H_t|, w_t=|H_t| and let D count deaths.
The exact telescoping identity2D=sum|H_t symmetric_difference H_(t+1)|
+w_1-w_L gives

    C_T+D >= (S+2(w_1-w_L)-random_loss-3n)/4.

The allocation is fixed independently of these supports, so the bound
holds simultaneously for varying-weight histories, including resets.
Adjusted activity at least nL/2 retains loss threshold371900. Aggregate
activity at least nL/2 with free endpoint weights instead requires370876.
The latter lowers the selected per-update allowance to7.654289121736758;
the checkpoint's existing worst upper7.651481649163582 remains below it.
Neither tail is uniform over all histograms yet. This is a geometric
extension, not a completed field-cancellation or minimum-distance proof.

The new deterministic tests exhaust147456 small histories across36 fixed
geometries, including135000 nonconstant-weight histories. All737354 checks
pass, as do the original4382 local checks, in0.905s. CWC separates the
arbitrary-width analytical inequality from these finite implementation
tests. The construction, potential and auxiliary transition kernel remain
unchanged. No native timing, library/paper edits, commit, push or promotion.
Goal active; all numerical work serial, bounded and file-free.

Next inspect central-envelope loss across the matching-exposure phase
boundaries before further axis subdivision. A phase-aligned partition may
retain the actual exposure counts more efficiently. Full auxiliary
coverage, other target-row counts, low-activity histories and the
field-cancellation composition remain open.

### 2026-10-26 — Central controls and bounded tails close a larger mixed region

Previous turn made progress: bounded sensitivity and two new regions.
This turn addresses central representation loss instead of only tuning
the tail. Extended triangulated controls to degrees3 and4. A retained
multi-index must have all active vertex pairs in the triangulation's edge
set. Every simplex-supported index remains, so the same barycentric
representation covers every histogram across every Gaussian annulus.
Extra retained cliques only weaken the bound. The degree-two path remains.

Added5477 exact blossom, reconstruction, simplex-support and exclusion
checks; the coverage suite passes13107 checks at40 and80 digits.
For the balanced64..128 region, degree-three restricted controls improve
the unrestricted7.659564 bound to7.628991, using136 rather than220 controls.
The previously failing degree-two bound was7.895777. This diagnoses
representation slack, not a construction failure.

The uncovered mixed region
([0,79],[80,159],[180,233],[72,252],[0,180]), sum512, now passes.
Degree-three central contribution7.838647 still fails; degree four reduces
it to7.597060108. The old tail remains too large, so integrated the previous
bounded switch result on radii16..32 into the same population controls,
retaining the global bound beyond32. Final upper7.597060108312 is below
allowance7.679295. The restricted basis uses670 controls instead of4845.
A40-digit run takes27.605s; an80-digit replay with tests takes57.299s and
agrees. Both are proof-check timings, not encoder performance.

Added this disjoint26,133,720-histogram region to the checkpoint. It now
contains22 regions totaling613,409,925, with2,306,209,020 pending in five
regions. Added explicit replay slices to keep serial runs bounded; each
slice reports its own verified count separately from the checkpoint total.
Fresh40-digit slices0:11 and11:21 replay the old21 regions in62.820and
55.795seconds. With the80-digit new-region slice, all22 regions have now
been rechecked. The original80-digit control passes5018 tests and retains
its7.672286056280 bound.
The selected target-row count remains1974. This is not a whole-code bound.

Construction, library and paper unchanged. Frozen helper hash unchanged;
disk above22GiB; numerical jobs serial and file-free. No benchmark,
commit, push or promotion. Goal active.

Next use elevated restricted controls plus the bounded annulus on the
remaining mixed regions. Compare central and tail contributions before
subdividing. Other target-row counts, arbitrary state histories, and
field-cancellation composition remain open.

### 2026-10-26 — Bounded sensitivity works; central controls still limit coverage

Previous turn was a status report, hence no new proof progress. Revalidated
the source checkpoint and resumed with a structural bound rather than
another blind subdivision run. The construction remains fixed.

Derived directed bounds on the switch gradient inside a Gaussian ball,
using actual outcome-probability ratios and their common tilt. Over all
centers and55 switches, squared-gradient upper bounds before beta^2 are
7.097414 at radius4 and14.073988 at radius16, versus the sharp unrestricted
constant32. The old constant remains on the infinite tail. The new helper
passes2358 rational/interval checks at40 and80 digits. An initial test was
stopped and repaired: rationalizing a rounded-zero Decimal endpoint was
creating a million-digit denominator. No numerical result from that
interrupted run was retained.

Added an opt-in bounded-switch scalar bound, including a16..32 annulus
before the original global tail. It certifies a previously failing region
of5,623,240 histograms, reducing its bound9.511594 to7.651482. Another
4,423,190-histogram region passes at7.379158 (the old scalar method also
passes that second region). Both new cases replay at80 digits. Added them
to the disjoint source checkpoint; a fresh111.220s40-digit replay verifies
all21 regions after7630 coverage checks. Total587,276,205 certified,
2,332,342,740 pending in six regions; worst upper7.651482 is below the
selected allowance7.679295. These are one-update bounds at row count1974.

The diagnostic also identifies where this refinement stops helping. On
the balanced64..128 region the scalar upper improves9.574624 to8.091269,
but the stronger degree-two matrix certificate already has central-shell
contribution7.895776209. Sharpening its tail cannot make it pass. Its
previous higher-degree certificate remains preserved. The neighboring
26,133,720-histogram region improves20.878262 to8.836943 and remains open.
These are proof-envelope failures, not code counterexamples.

CWC keeps local sensitivity, partial coverage, and distance claims
separate. The original80-digit control retains its7.672286056280 bound
and5018 passing checks. Frozen helper hash unchanged; disk remains above22GiB. Jobs
were serial, bounded and file-free. No encoder benchmark, library/paper
change, commit, push, or promotion. Goal remains active.

Next use the bounded-radius test for cheap coverage, but focus the hard
mixed populations on sharper central controls. Preserve the common
histogram dependence when combining shells. Other target-row counts,
arbitrary state histories and field-cancellation composition remain open;
do not present this coverage count as proximity to whole-code closure.

### 2026-10-26 — Continued coverage and a sharp sensitivity obstruction

Previous turn was progress: a replayable partial partition and a cheaper
matrix trigger. Continued from its complement at80 digits, certifying
78,927,954 additional histograms in a bounded100.002s search. No failed
singleton; one unfinished deadline remained pending. Added the nine
passing regions to the compact source checkpoint, preserving exact
nonoverlap and complement accounting.

Added opt-in40-digit directed arithmetic without changing the default or
80-digit mode. The7630 coverage checks and5018 original checks pass at40
digits. A serial fixed-case comparison takes9.595s at80 digits and5.835s
at40, with upper-bound agreement better than10^-25. A fresh104.881s
40-digit replay verifies all19 checkpoint regions:577,229,775 histograms
pass and2,342,389,170 remain pending. This is partial one-update moment
coverage at target-row count1974, not a full histogram or distance theorem.

Also tested a structural improvement instead of only subdividing boxes:
bound each matching switch using its common positive variables. A FLOAT
polynomial screen covers55 switches, with190 exact rational identity
checks. An exact leading-monomial witness shows that the existing global
Lipschitz constant4*sqrt(2)*beta is sharp. The switch
log F22-2 log F24+log F44 approaches gradient beta*(0,0,0,4,-4) along
x=(1,epsilon^3,1,epsilon^-1,1). This rejects a smaller unrestricted
constant, not the construction. The extreme Gaussian direction's mass
must remain in any sharper bound.

CWC separates the numerical screen, exact limiting argument and partial
certificate. Frozen helper hash unchanged; disk above22GiB. All jobs
serial, bounded and file-free. No native benchmark, construction/library/
paper changes, raw dumps, commits, pushes or promotion. Goal active.

Next continue the saved complement at40 digits and test radius-dependent
common-variable switch bounds against the local matrix envelopes. Do not
retry a smaller uniform global sensitivity constant. Other row counts,
arbitrary histories and field-cancellation composition remain open.

### 2026-10-26 — Larger corner regions and a replayable retuned partition

Previous turn was progress: triangulated bounds added coverage and
identified representational slack. Integrated them into the accounted
full-domain pilot. Its old-cap run certifies278,945,607 histograms,
including41,964,450 beyond rechecked caps. Timeouts remain pending.

Expanded the five disjoint corner regions to thresholds
(320,320,288,288,288), with polynomial degrees(2,2,4,3,3). All pass under
the same retuned weights; they cover448,226,845 histograms versus the
old236,981,157. Retained degree-two failures for the last three regions
and the coordinate2 degree-three failure. These degrees belong to the
proof representation, not the encoder.

The new-cap pilot initially added no interior coverage: its width160
matrix trigger spent the budget on oversized regions. Width96 instead
closes four boxes containing47,851,881 histograms in the same roughly
100-second work budget. A pending-only one-node continuation adds another
2,223,095. Its reporting separates prior certificates not replayed from
newly verified coverage.

Saved a compact source partition of the five corners and five interior
boxes. A fresh89.063s directed replay verifies all ten:498,301,821
histograms certified,2,421,317,124 pending in six exactly accounted regions.
Worst checked log upper bound is7.634662, below allowance7.679295. This
is a partial one-update certificate at target-row count1974, not a uniform
histogram theorem or whole-code failure bound. The checkpoint does not
silently include other earlier certified boxes.

The updated exact suite passes7630 checks. The original-point control
replays5018 checks and retains its7.672286057 bound. Frozen diagnosis hash unchanged;
disk remains above22GiB. All numerical jobs were serial, bounded and
file-free. CWC kept proof degree, population coverage and code claims
distinct. No construction/native performance/library/paper change, raw
dump, commit, push or promotion. Goal remains active.

Next use --retuned-pending-cover on the saved complement instead of
rediscovering passed regions. Prioritize large mixed-population regions;
only attempt matrix bounds after sufficient subdivision. Other target-row
counts, arbitrary histories and field-cancellation composition remain open.

### 2026-10-26 — Triangulated controls separate regional slack from point behavior

Previous turn was a status report, not new proof progress. The prior live
handle17021 is now missing, and no Python process was present. Its final
output was unavailable, so no final coverage is inferred from that run.

Preserved the construction and retuned potential. Added exact pulling
triangulations of the population polytope, retaining only vertex-pair
controls needed by its simplices. Shared controls remain fixed across
Gaussian shells and the tail. Added optional interior-point cones and
integral edge bisection, without changing the covered region.

- The merged8,327,716-histogram region passes below6.654506 using157
  controls instead of378; its check takes15.228s. A rectangular81-control
  bound finishes faster but fails at7.838282. An unrestricted replay hits
  its28s deadline; that timeout is not a mathematical failure.
- Extending h0<=43 to h0<=63 closes10,794,081 histograms below6.713346,
  adding2,466,365 beyond that merged region. No new full-domain total is
  claimed without replaying its entire partition.
- On the already degree-three-certified balanced band, degree-two pulling
  gives7.895777. Interior-point cones do not help. Two targeted edge splits
  improve this to7.764475, still above7.679295. A troublesome edge midpoint
  separately passes below6.999814. This identifies representational slack
  at that location, not a code obstruction or a uniform pointwise theorem.

Final exact suite:7627 checks, including exhaustive small lattice coverage,
independent slice-volume identities, polynomial reconstruction, degenerate
fallbacks and edge splits. The zero-dimensional star case was repaired
before replay. Original-point5018 checks and the old seven-region partition
still pass. Frozen helper hash unchanged; free space above22GiB. Source and
compact notes only; no native benchmark, production/paper change, commit,
push, raw output file or construction change. CWC keeps the scope explicit.

Next integrate triangulated controls into the accounted cover with bounded
edge refinement and ordinary subdivision. Do not keep polishing the already
closed balanced band. Uniform histogram coverage is unfinished, and other
row counts, general histories and field-cancellation composition remain.
The original full research goal remains active.

### 2026-10-26 — Population-aware radial controls and a complete-domain pilot

Previous turn was progress: the broader upper-integrand search supported
the fixed retuned weights. Implemented shared polar controls for the core
and gradient on the actual clipped population polytope. Convexity after
Gaussian integration gives a cheap scalar bound. Added an optional exact
degree-four product of the gradient polynomials to retain signed cross
terms. Both preserve the construction and fixed proof potential.

The low-h0 broad-region bound improves49.451373 to31.358239, still failing.
For the64..128 band, the scalar bound improves15.213360 to9.802749;
degree four gives9.765471, still failing. Splitting Gaussian radii into
shells gives9.574624 there but closes the existing593,589-histogram boundary
box at7.549269. This distinguishes useful domain tightening from the
remaining need for matrix/directional information. Do not claim scalar
bounds now solve the entire domain. Degree-four cases exceeding2000
controls are explicitly inconclusive; no guard was raised.

A fresh complete-domain pilot rechecks the five caps using the retuned
weights and subdivides their disjoint complement. The first short run
adds no coverage. The longer24-node run closes five new regions containing
15,372,391 histograms. It certifies252,353,548 including rechecked caps,
with2,667,265,397 pending in15 regions and no failed singleton. Counts sum
exactly to2,919,618,945. The run takes83.655s under a100s deadline, with
nonblocking output polls. No search tree or raw dump is saved. The new
regions are also disjoint from the previously retuned band and boundary
box, but the pilot does not silently import those separate certificates.

The checker now passes1060 exact/domain/Gaussian checks, including independent
rational expansion of every degree-four product coefficient in a test
case. The old-point80-digit verifier still passes5018 checks and gives
7.672286057. Frozen diagnosis hash unchanged; disk above22GiB; all jobs
serial and file-free. No native encoder timings, construction or paper
changes, commits or pushes. Whole-code proof remains open; goal active.

Next combine the cheap scalar pass with matrix verification for difficult
leaves and avoid repeated switch/envelope preparation across calls. Measure
coverage and verifier cost on a bounded continuation before a larger run.
Further scalar degree increases are not justified by this diagnostic.
Uniform histogram coverage is only one obligation: other row counts,
general support profiles and field cancellation remain after it.

### 2026-10-26 — Broader upper-integrand search supports the retuned potential

Previous turn was progress: new weights and directed regional checks.
Kept those weights fixed and searched the numerical upper integrand on
the coarse simplex grid, then locally on every boundary face. The bounded
run examines863 distinct histograms and takes16.580s. Its largest order3
value is7.059943 at(191,167,113,41,0); orders5 and7 give7.060040 and7.060025.
The allowance remains7.679295. Other face maxima are lower:6.850945 when
h0=0,6.589195 when h1=0,5.398263 when h2=0, and6.927041 when h3=0.
These are numerical searches, neither exhaustive nor directed integration.

No displaced numerical peak was found. This supports investing in a
uniform bound under the new potential, rather than changing the encoder.
Next examine the coarse region envelopes' use of rectangular populations:
some combinations are impossible under sum h=512. A candidate improvement
is a shared polar representation of the core and gradient on the actual
population polytope, with convexity used before their Gaussian integration.
That improvement has not been implemented or certified in this turn.

The source command is random_pairing_flow_quadratic.py --retuned-upper-search.
Disk guard passed, all runs serial and file-free. Construction unchanged;
no native timing, full histogram certificate, or whole-code closure. Goal
active. The user requested a status explanation before the next change.

### 2026-10-26 — Retune the proof potential before refining more boxes

The broad coverage diagnostic split21 pending regions in35.835s and closed
none. The broad radial envelopes lose too much across distinct population
patterns; more blind subdivision is not the preferred next step.

Accelerated the Gaussian-center proposal without changing its equation.
The mean update is2*kappa-Lipschitz, hence contractive at kappa=.08.
Direct iteration computes only first moments; second moments are evaluated
once afterward. The old damped solver remains as a reference. On72 points,
proposal probabilities differ by at most3.32e-13; the solve took0.166s
versus2.998s. The exact140 finite-population moment tests also pass.
These timings concern proof tooling, not encoding.

Used this speedup to explore the linear weights in the quadratic potential.
The construction, theta=.1, eta=.2, kappa=.08 and target-row count1974
remain fixed. A bounded1150-point lower-bound screen suggests the separate
candidate u=(0,-79200914,-154166955,-195740597,-227550536)/1e9.
It is not a uniform upper bound. Subsequent80-digit directed checks give:

- At(189,168,115,40,0), log M<=7.083655, allowance>7.679295226742.
- The previously failing64..128 population band closes at7.659564 using
  degree3 controls, covering8,250,385 histograms. Degree2 still fails.
- A593,589-histogram boundary box closes at7.464759 with degree2 controls.
- All five old caps reverify under the new weights, covering236,981,157
  histograms. The three sets of regions are pairwise disjoint.

These are selected one-update moment bounds under one new potential, not
whole-code closure. Do not combine old- and new-weight regions as one
telescoping certificate. Old weights remain the default, and each verifier
result now records its exact weights. The original partial partition is
replayed separately: all seven regions pass in37.006s. The original point
replay passes5018 checks at80 digits and remains below7.672286057. Six
additional interval checks verify invariance when the same constant is
added to every potential weight, including the endpoint allowance. The
frozen diagnosis hash is unchanged; free disk remains above22GiB.
Source and compact results only; no native timing,
encoder changes, promotion, commits or pushes.

Next validate the new weights against a broader upper-integrand search,
including all boundary faces, then build coverage under that single
potential. Further weight optimization may be cheaper than fine partitioning.
Full histogram coverage, other target-row counts, general support profiles
and field-cancellation accounting remain open. Goal active.

### 2026-10-26 — Global outcome geometry closes a band; explicit partial partition

Previous turn was progress: joint Gaussian controls and selected passing
regions. Added exact degree3/4 elevation of polytope controls. For the band
80<=h_i<=136, degree3 puts the finite-annulus bound below the allowance,
but the old tail still fails. Derived a global covariance matrix from each
pair's distinct outcomes: Cov(v)<=V V^T/2, using the categorical bound
diag(p)-p p^T<=I/2. This holds for every Gaussian tilt, not just near zero.
The mean matching curvature retains exact quadratic population counts.

Replacing the isotropic tail curvature by this matrix reduces its log
control bound from8.237293 to-1.180483. The complete band bound now passes
at7.672986463 for4,884,615 histograms. The larger64..128 band still fails:
degree3 gives7.754862, degree4 gives7.712463. An alternate exposure order
(4,3,2,0,1) gives7.715887 at degree3. Its remaining loss is local, not the
tail. Stop raising degree blindly; refine the population/exposure regions.

Added exact disjoint slab subtraction and a fixed partial partition with
five caps, the newly closed band, and B_4. Its80-digit replay verifies all
seven regions:241,866,501 histograms certified;2,677,752,444 pending in25
regions. Counts exactly cover the2,919,618,945-point domain. This is a
partial one-update moment certificate, not a full histogram or code proof.
No stored passing labels are trusted by the replay.

The helper passes1040 checks, including new rational covariance and
translated-matrix checks, elevated polar reconstruction and exact partition
membership. Combined replay37.630s under55s deadline; the single degree4
failed control took46.049s under75s. Deadline checks occur inside the control
loop, with at most2000 polar controls. All runs were serial and file-free,
above20GiB free. No construction, native timing, library, paper, commit or
push changes. Derivation and commands are in RANDOM_PAIRING_FLOW_CONTROL_POINTS.md.

Next refine the explicit pending regions with the appropriate bound and
exposure order, while retaining exact accounting. Whole-code closure also
needs other target-row counts, support profiles and field-cancellation
counting. Goal active; no construction obstruction or redesign decision
has been established by these failed upper bounds.

### 2026-10-26 — Joint integral controls preserve one histogram across shells

Previous turn was progress:236,981,157 histograms covered by five caps,
and a verified scalar-versus-matrix distinction. Tested a common lower
precision matrix over each population region. It certifies a421,750-point
interior box, but does not improve the tight boundary enough.

Derived a stronger interface from convexity of the Gaussian integral in
its exponent coefficients. Expand the constant, linear and precision terms
in one Bernstein basis; integrate each control exponent. Crucially retain
the same population weights across annuli: sum each control's contributions
before maximizing. Degree elevation improves the bound without subdivision.
The final80-digit replay closes the unsplit729-histogram boundary B_4 at
7.677136805, versus allowance>7.679168691101. Its separate-shell bound was
above the allowance. This removes a concrete correlation loss in the proof.

Added a quadratic polar representation on the actual population polytope,
avoiding impossible points from dropping the eliminated population's bounds.
It certifies the871,720-histogram band88<=h_i<=124, sum h=512, below7.485025.
Wider bands fail:80..136 gives8.618565;64..128 gives10.172558. These bounds
identify both local-envelope slack and displaced-center tail loss. They do
not exhibit low-distance codewords. Improved the region tail using the
previous global coupled bound and a joint-control tail; no local bound is
used outside its valid radius. One early broad run rejected its old tail
coefficient; no invalid result was accepted.

The helper passes778 exact domain, polynomial and Gaussian controls,
including polar reconstruction, degree elevation and analytical integrals.
Final boundary replay11.936s; passing polytope band9.199s. Runs remain
serial, bounded and file-free under the20GiB disk guard. New derivation and
failed controls are indexed in RANDOM_PAIRING_FLOW_CONTROL_POINTS.md.
No construction, library, paper, timing, commit or push changes.
The final regression preserves all five cap certificates and the original
matrix point bound; its5018 checks pass. Frozen diagnosis hash is unchanged.

Next construct an explicitly accounted interior partition using these
complementary bounds, and diagnose broad-region center displacement rather
than blindly increasing polynomial degree. These regions overlap; they
are not a replacement for full coverage. Other target-row counts, support
profiles and field cancellations remain required. Goal active.

### 2026-10-26 — Retain population dependence through Gaussian integration

Previous turn was progress: exact coverage accounting and coarse vertex
bounds. This turn diagnoses two distinct losses rather than enlarging a
blind point search. First, maximizing the Gaussian core and its gradient
separately loses cancellation. Integrating the linear Gaussian term first
retains that cancellation in one degree-four population polynomial.
Second, global switch Lipschitz bounds are too loose in the interior.
Added ball-local curvature bounds with an independent global tail, retaining
the same population dependence in that tail too. The encoder is unchanged.

Five disjoint caps, with dominant thresholds(352,352,320,320,320), now
certify236,981,157 auxiliary histograms. Their largest log bound is below
5.896333, versus allowance>7.679168691101. The80-digit bounded global replay
verifies those caps and leaves2,682,637,788 histograms pending, with no
unclosed singleton, after12 attempts/18.196s. This is about8.12% coverage,
not a uniform moment or whole-code certificate. Exact accounting is intact.

A three-shell scalar bound also certifies the16,730-histogram interior
region96<=h_i<=108, sum h=512, below7.665537. It fails at the known difficult
single histogram(185,188,83,56,0):7.851141. The existing anisotropic matrix
bound still certifies that point at7.672287. Thus subdivision alone cannot
make this scalar envelope sufficient everywhere. Preserve the matrix
method near the tight boundary and use cheaper bounds where they suffice.

The cover now passes711 exact domain, polynomial and Gaussian-integral
checks, including an exact core/gradient cancellation example and central
Gaussian tail controls. The separate80-digit matrix regression passes5018
checks and reproduces7.672286056279. All runs were serial, bounded, and
file-free with more than20GiB locally free. No native timing, production
change, paper edit, commit or push. Frozen diagnosis helper is unchanged.

Next use the resulting hybrid proof interface to classify the remaining
interior, retaining matrix directions in tight cases. A promising targeted
improvement is to bound its matrix variation without losing population
correlations, rather than shrinking every region to singleton histograms.
Other target-row counts, support profiles and field cancellations still
remain outside this selected one-update calculation. Goal active.

### 2026-10-26 — Accounted simplex coverage and a coarse radial bound

Generalized the verifier to five inclusive population bounds constrained by
sum h=512. Exact inclusion-exclusion counts each region; disjoint splits
preserve the count. The full domain has2,919,618,945 auxiliary histograms.
These populations describe the causal-flow proof, not field-symbol values.
The new bounded cover checks certified+pending+unclosed=total throughout.
Its696 domain, coordinate and polynomial checks pass.

Derived a cheap all-Gaussian radial envelope retaining rare-color exposure
counts. A simplex Bernstein bound avoids the impossible corners introduced
by a rectangular relaxation near a vertex. Five disjoint vertex regions
now certify16,499,285 histograms for the selected one-update moment.
An80-digit replay takes3.190s and reproduces the passing caps and failing
larger controls. A bounded global run certifies16,505,404, leaves
2,903,113,541 pending, and reports no failed singleton. This is only about
0.57% coverage, not a uniform moment or distance certificate.

Also checked a489-histogram face with an absent auxiliary state and two
70-histogram vertex neighborhoods. Directed positive-curvature witnesses
reject joint concavity of the natural Gaussian core, including a fixed-fifth
slice. They reject that proof shortcut, not the code or desired moment bound.

Local disk recovered above20GiB. Numerical jobs remained serial and bounded;
one earlier local box replay completed its checks but exceeded its25s final
wall guard, so it was not accepted as a successful run. Later bounded runs
completed. No raw output files, construction changes, native benchmarks,
commits, library edits or paper edits. Frozen diagnosis hash unchanged.

Next inspect the interior regions where the coarse bound loses slack,
retaining joint histogram/Gaussian dependence where needed. Do not replace
this diagnosis with blind enumeration of billions of points. Even full
histogram coverage here would address only the selected target-row count;
other support profiles and field-cancellation counting remain open.

### 2026-10-26 — Local boundary covered; generalized phase-independent regions

Previous turn was progress: directed histogram boxes, with one unsplit
zero-fifth-color box failing. Covered that125-histogram box by eight radius1
boxes; all pass. The initial uniform bound over those pieces is7.678944.
No extra union factor is needed for a cover of conditioning histograms.

Derived and verified a sharper shell bound by subtracting a lower bound on
Gaussian mass beyond each upper radius. Both spectral bounds are directed.
The selected point now has log moment<7.672287, versus allowance>7.679168691101.
Its hypothetical40.66-bit slack is NOT a certificate: this point bound is
not uniform over all histograms. The full unsplit125-point box still fails
at7.679669; its subdivided cover remains valid.

Generalized histogram centers and exposure priorities. For regions crossing
an exposure-phase boundary, derived a nested-matrix affine majorant using
monotone squared ranges and capped prefix populations. Its extension is
valid across the boundary even when its formal final quota becomes negative.
The proof does not replace bounds inside negative actual quotas. Interval
LDL checks the nested matrices. Added exact phase and Gaussian-shell tests.

Four off-neighborhood region checks pass, including a729-histogram box
crossing exposure phases (log bound4.992448). The80-digit phase replay takes
10.331s and passes5018 tests. The final boundary replay takes17.820s and
certifies its125-point cover below7.676047. Updated the existing verifier and mathematical
note, not the construction. Frozen diagnosis helper hash unchanged.

All numerical checks remain serial, file-free on Peach under30s/256MiB caps;
/tmp has63GiB free. Another account build was present; no encoder timing or
benchmark was run or inferred from these proof-check runtimes. Local space
remains below20GiB. No raw files, commits, library or paper changes.

Next build an explicitly accounted complete histogram cover, including
simplex faces; selected boxes are not enough. Other target-row counts,
support profiles and field-cancellation counting still follow. Goal active:
this turn closes the previous local boundary and removes a structural
restriction of the region proof, but does not claim whole-code closure.

### 2026-10-26 — Certified histogram boxes; boundary coverage remains open

Previous turn was progress: a directed local integral, including its infinite
tail. This turn checks all16 single-coordinate neighbors; all pass, including
introduction of the fifth auxiliary color. The worst neighbor bound is below
7.675895. Generalized the same verifier to whole constrained histogram boxes.

Plain interval arithmetic loses cancellations and fails tiny boxes. Rewrote
the exact matching mean, gradient and curvature as quadratic polynomials in
four independent histogram coordinates. Retained the affine exposure quotas
jointly, then combined the potential with the Gaussian shift before bounding.
Tensor Bernstein coefficients bound the resulting degree4 polynomial. This
is region coverage, not an extrapolation from sampled points.

Directed80-digit results: boxes B(r,l,u), defined in the quadratic note,
(1,0,0):27 histograms, bound7.678556; (2,2,2):125,7.660477;
(4,4,8):3645,7.650837; (8,8,16):44217,7.668881. All are below the directed
allowance>7.679168691101. Control B(2,0,0),125 histograms, still fails the
bound at7.682042. Boxes overlap and do not cover the simplex. No whole-code
or structural-tail certificate is claimed. Recentring the auxiliary Gaussian
representation was essential for the largest box; the encoder is unchanged.

Box replay10.448s; final combined neighbor/box replay21.653s,399
arithmetic/polynomial checks, all expected passes reproduced. Exact rational
tests check the quadratic expansion and its square; interval LDL verifies
all selected matrices, solves and determinants. Fixed overly sharp shell
eigenvalue candidates for interval boxes, and a zero-population guard that
had rejected harmless outward-rounding below zero. No invalid check was
accepted. Numerical jobs were serial/file-free on Peach under30s/256MiB
limits,63GiB scratch free. Local disk remains below the guard. No raw files,
native benchmarks, commits, production or paper changes. Frozen helper
SHA256 remains7238FEC9E5553E7B4D1CE3EDB3089CEC06ABCC34CB090C029D0F3D6B54BBAED3.

Next focus on the zero-fifth-color boundary, where the same box method still
loses the small slack, and extend the proof to a complete histogram cover.
The full objective additionally requires other row counts/support profiles,
field-cancellation counting and eventual implementation timing. Goal remains
active; this turn adds directed region coverage and a discriminating failed
control, rather than merely sampling more points.

### 2026-10-26 — Directed Gaussian envelope passes the selected histogram

Previous turn supplied numerical slack, not a certified integral. This turn
derives a finite analytical envelope and verifies it with outward-rounded
Decimal arithmetic. At h=(185,188,83,56,0), the actual uniform-matching/IID-tag
one-update log moment is below7.675753. The directed allowance is above
7.679168691101. This is a local statement, not uniform histogram coverage or
a distance certificate. If uniform, this bound would leave20.18 bits in the
selected structural calculation; that premise remains unproved.

The first radial envelope failed (>7.74). Retaining affine switch directions
in verified ellipsoids, and bounding tilted covariance by its individual
centered contributions, removes enough loss. Gaussian quadratic integration
covers a radius4 ball and five shells through16. The prior global bound
covers the infinite tail. Matrix positivity, linear solves, determinants,
elementary functions and accounting costs are checked with intervals. No
quadrature agreement is used as an error bound.65/80-digit replays agree;
the latest replay passes183 arithmetic checks and takes0.616s.

Added two small helpers and extended RANDOM_PAIRING_FLOW_QUADRATIC.md.
All numerical work ran serially, file-free on Peach (63GiB free, account
idle before starting), with30s/256MiB caps. Local space remains below20GiB.
No native timing, construction change, raw outputs, library/paper edits,
commit or promotion. Frozen diagnosis helper hash is unchanged.

Next test and bound neighboring histogram regions, especially a nonzero
fifth color, rather than infer a global maximum from the numerical screen.
The broader remaining tasks are other row counts/support profiles and
field-cancellation counting. Goal active; this turn made a concrete proof
bound advance rather than merely extending a parameter sweep.

### 2026-10-26 — Quadratic candidate gains numerical slack; certify next

Previous turn was progress: full-width block interface and diversity
candidate. Tested mixed histograms, not just homogeneous starts. Exact
finite-population first/second moments pass140 cases. The initial kappa.01
potential fails an exact Jensen witness (>35 nats after conservative union
costs). Larger kappa survives; kappa.08 is the current candidate, with fixed
rational u=(0,-84450914,-150416955,-196990597,-222300536)/1e9, theta.1,eta.2.

Derived an actual-matching upper integrand using a four-dimensional Gaussian
representation and Doob/Hoeffding exposure. Exposing rare colors first
improves the conditional-range budget;2520 exact path checks pass. The
analytical bound is valid, but numerical quadrature is not yet certified.
592 grid/refined histograms select h=(185,188,83,56,0). Order5/7/9 log
integrals converge near7.648277, below7.679169 after endpoint, initial-batch,
row-union, and regional-conditioning costs. If made uniform and directed,
this would leave about182.5 bits for the selected structural tail. It is
NOT a certificate: integration error and untested histograms remain open.

Derived a global radial Gaussian tail envelope; at the selected point its
floating log absolute tails are-8.8275 (R16),-39.0334 (R24). Next certify the
central integral and cover histogram regions, retaining the exact matching.
Whole-code closure still additionally needs other rows/profiles and field
cancellation. No native timing or construction change.

Added RANDOM_PAIRING_FLOW_QUADRATIC.md and one bounded helper. Local25.19GiB
became insufficient; guard stopped execution. Used file-free Peach with
63GiB /tmp, idle account,30s/256MiB caps,Python-B. Frozen-parameter full screen
14.027s, serial, no raw files/packages. Corrected an overly strict odd-moment
roundoff test before the successful three-order run. Old proofs preserved;
no library/paper/commit/promotion. Goal active with new discriminating evidence.

### 2026-10-26 — Full-width two-update envelope; retain histogram diversity

Previous turn was progress: verified separable-tail barrier and exact
three-update mean identity. Derived a valid all-histogram two-update
envelope via pairing polynomials, an entrywise PSD majorant, Gaussian
quadratic moments and Holder. Its evaluation uses5-by-5 powers through
degree256, not full histogram enumeration. Exact width4 checks cover all70
histograms;15 fourth-moment identities also pass.

The512-wide floating screen fails: selected union exponent+7689 versus
one-step+6026 at the chosen markers. Even exact homogeneous-start moments,
with endpoint weights tuned over9 marker pairs, give best sampled optimistic
exponent+3664. This is not a global impossibility theorem. It diagnoses
loss from allowing homogeneous histograms at every block boundary rather
than carrying diversity across blocks. Do not just polish the PSD shift.

Proposed exp(kappa sum h_a^2/n) times the product potential. Derived and
checked its exact homogeneous-start regional mean change; kappa.01 gives
about-2.71/-2.49 nats. Mean drift is not an exponential tail. Next test
mixed histograms and the actual conditional moment, retaining the mixture
normalization instead of taking separate worst-case product bounds.

Added RANDOM_PAIRING_FLOW_BLOCK.md and bounded helper. Native timing and
whole-code closure remain open. Local disk20.83GiB before runs; serial,
file-free Python-B. Full screen1.636s; nine-point rerun15.263s; final exact
suite0.094s with5 additional regional collision checks. Fixed a misplaced
reporting call after an initial ladder run; successful rerun covers all
reported ladder data. No construction/library/paper/commit/promotion change.
Goal active: concrete proof-method evidence, not a completion or blocker.

### 2026-10-26 — Verified separable-tail barrier; exact three-update gain

Previous turn was progress: capacity-sharing and exact transition kernel.
This turn tests a tail template rather than doing more local allocation
tuning. A rational stationary-flow witness proves any time-independent,
separable, all-histogram one-step envelope followed by the separate target-
row union misses by >1019 bits, even before regional conditioning. The
witness has exact mean tag count and sufficient loss; KL cost is bounded
by an exact chi-square expression <5500 bits, versus row union >6519 bits.
This excludes that proof template, not the construction or all tail methods.

Derived an exact delayed benefit from uniform re-pairing. The first two
loss means ignore pair diversity; the third retains a(a+b-1)/(n-1), with
a,b counting auxiliary states-2,-1. Extended to actual remaining regional
counts through six-position falling-factorial probabilities. This is a
mean identity, not a tail or minimum-distance certificate. Next use a
multi-update envelope retaining histogram structure and regional counts.

Added RANDOM_PAIRING_FLOW_TAIL.md and a bounded spectral-screen/exact-witness
helper; extended the existing kernel test. Exact barrier check0.005s;
kernel suite0.895s, including5 IID and9 regional mean checks. Local space
22.56GiB restored; Python-B serial/file-free runs, no native benchmarking.
No construction, library, paper, commit, or promotion changes. Goal active;
whole-code closure and implementation cost remain unproved/unmeasured.

### 2026-10-26 — Capacity sharing and exact finite-population flow kernel

Progress, not whole-code closure. Demand-dependent sharing improves the
selected regional 32-update flow lower bound from1280.75 to1483.75.
Derived a causal five-state auxiliary flow rule, exact loss conservation,
and non-amplification of one-unit perturbations. The R256 diagnostic now
asks for a tail on random-update loss at371901; this is a sufficient-flow
failure threshold, not a codeword criterion.

IID mean-field loss is about12.34%; fixed-partner IID control is15.827%
with exact rational stationarity. Neither is the actual finite ensemble.
Derived the actual uniform-matching histogram transition with regional
remaining-label counts. Exact tests pass36 matching profiles and8 complete
small transition kernels against labeled brute force, plus2916 threshold,
1215 allocation,81 cut,170 stability checks. The global tail remains open.

Recorded RANDOM_PAIRING_FLOW.md and two compact standalone helpers.
Tests ran serially and file-free on Peach (63GiB free, account idle),
0.293s under30s/256MiB bounds. Local disk was below20GiB guard. No new
native timing, construction changes, paper edits, commits, or promotion.
Next: rigorous moment envelope for the finite-population histogram process,
starting at1974 target rows and paying the all-target-set union. Preserve
the restricted certificates and do not enumerate the full histogram space.

### 2026-10-26 — Scope extension and a joint-dual barrier

Previous turn was progress: verified fixed-T cut bound g_T<=-47942.
Derived deterministic sensitivity g_T'-g_T<=128|T triangle T'| atR256.
This certifies g<=-70 for every T' within374 membership changes (187 swaps
at equal size). Extreme target sizes<=655 or>=7537 have g<=-47 independently
of setup. These exclude only the sufficient rank-gap obstruction, not bad
words in general. Exact integer counting shows a same-size neighborhood
cover at the current margin would require at least2^4425 centers.

Proved a limitation of the obvious joint graph-cut extension. Relaxing
both half-full layers and target size through linear prices admits a
mixture of empty/full states with anticorrelated target sizes. Any such
dual lower bound is <=-128 min(s,8192-s), or-252672 at s1974, too weak
against the-83887 threshold. Per-layer prices and death reweighting do
not remove this barrier. It is a proof-relaxation obstruction, not a code
counterexample. Do not build a large joint-cut optimizer expecting tuning
alone to fix it.

Screened an explicit equal-split flow as an analytically simpler control:
each target gate assigns half service to each incoming coordinate timeline.
Corrected an initial DP normalization bug, then passed486 exhaustive
timeline checks. At supply1/4, IID extrapolated means164116 (window32) and
177284 (window256) are below the adaptive cut's216727; no actual setup
tail or distance claim follows. The ability to share gate capacity appears
important; next investigate a demand-dependent local allocation interface.

Added RANDOM_PAIRING_CUT_SCOPE.md, exact scope test and service diagnostic.
Scope tests pass4096 sensitivity,360 mixture,375 dual-price checks, with
toy hard minima0 versus negative relaxed barrier. Disk22.35GiB initially
dropped to12.37GiB; local guard prevented tests. Used connect-to-peach skill
and file-free stdin after checking63GiB /tmp and no numerical account job.
Serial standard-library jobs,30s timeout,256MiB cap,Python-B,remote20GiB
guard; tests0.10s, corrected screen0.04s. No remote packages/files or local
raw data, no native timing, code construction, paper, commit, or promotion.
Goal active: analytical progress and remote test fallback, not blocked.

### 2026-10-26 — Verified fixed-set cut bound; preserve half-fullness

Previous turn was progress: bounded rank search and exact joint geometry
moment. Implemented the proposed directed-cut lower bound for coverage
of fixed target rows plus deaths. OR nodes encode live-pair/next-state
coverage; unary prices relax half-full activity. Verified integer flows
for capacity, conservation, antisymmetry, and residual cut equality.
Exhaustive small comparisons pass72 cases across both death formulations.

At seed20261027, choose the1974 lowest-coverage rows of the constant-mask
history. Complete128-window replay proves C_T+D>=194833.5 for ALL half-full
histories for this fixed T and setup, giving g_T<=-26048.5. The minimizing
cuts exploit sparse early/dense late states. Replace directed deaths by
half symmetric differences, exactly equal on half-full histories. Complete
replay strengthens to216727 and g_T<=-47942. Runs15.60/23.68s, serial.
This is deterministic fixed-set geometry coverage, NOT all-target-set,
setup-probability, field-failure, or whole-code distance closure.

64-update attempt at17/64 price verified windows0..27, sum97239.875,
then hit25s budget. Continuation from28 exited with allocation failure;
do not use any unreported results from that process. Later resource check:
disk17.55GiB (<20GiB), physical8.38GiB, virtual11.45GiB free. Cause of
allocation error not established. Stopped numerical work; no deletions or
system/process changes. Added per-flow disk guard and graceful partial
reporting on MemoryError, not numerically rerun under the low-disk guard.

Added RANDOM_PAIRING_CUT_BOUND.md and cut/test/diagnostic sources. Frozen
helper hash unchanged. No raw data, native timing, construction/library/
paper changes, commits, or promotion. Goal active: substantive fixed-set
bound established, with resources stopping the optional larger-window run.

Next address the union over target rows rather than merely testing more
selected sets. A weighted death price also permits a fixed-B lower bound
C_T+D>=L_eta+(1-eta)B and may tighten the critical profile. Need a valid
family/setup-probability interface and field count, not just more cut
instances. Recheck the20GiB guard before any further numerical work.

### 2026-10-26 — Actual pairing rank search and joint geometry check

Previous turn was progress: birth negative result and a structural rank
criterion. Implemented a bounded full-size route-and-pairing search at
seed20261027, keeping every state half-full. Fixed target-row coverage plus
deaths is optimized by exact blockwise swaps, then target rows are updated.
Verified2240 swap deltas and140 descent steps before full-size execution.

Constant-start block64/16/4 run18.34s found no improvement over rank gap
-83887. Independent block4 start reached89408 swaps, J412590, d501998,
z173295, gap-173295 at the25s cap. All row coverages>=135, zero parity
deficit; lowest1974 rows mean155.096. No positive full-size gap found.
Bounded one-seed local search is not a certificate or global optimum.

Derived joint local polynomial (x+u)^256(h+u)^k(1+u)^(256-k), retaining
deaths and outside crossings with the actual matching law. Conserved
target-label routing composes through h=(1-p+ptv)/(1-p+pt). Exact tests
pass8 polynomial and30 routing cases. SelectedR256,B131072,s1974 structural
event bound remains+178752 FLOAT; corresponding field term+309953. One
smaller profile's structural screen is negative, but that is neither a
directed certificate nor a distance guarantee. Do not replace the prior
proof with these selected checks.

Added RANDOM_PAIRING_RANK_SEARCH.md, search/test and joint-moment/test
helpers. No raw files, concurrent numerical runs, native encoder timing,
construction/library/paper edits, commits or promotion. Initial26.68GiB
free. Previous restricted certificates preserved. Goal remains active.

Next implement a small-window directed-cut lower bound on target-row
coverage plus deaths. OR coverage and directed death costs are submodular;
a Lagrange price for half-full state supplies a bound for all histories,
not merely found histories. Start with exact small checks, then fixed
target sets from the full-size sample. A fixed-set bound alone is not the
all-target-set or whole-code proof; it is the next obstruction diagnostic.

### 2026-10-26 — Birth constraints are insufficient; derive a rank-obstruction test

Previous turn was progress: directed16384-swap restricted coverage. This
turn kept the code fixed and tested the promised empty-pair birth condition.
With exact S, an inactive-row birth must be a crossing endpoint. Retained
it through three-category route conditioning and two conservation markers.
Exact moment checks pass. Selected R256,B131072 bound improves only297971
to296463 FLOAT bits; added conditioning even worsens smaller B profiles.
Old bounds remain available. An optimistic, NOT-code-bound model forcing
all births and deaths onto S still gives282518 in that case. Do not claim
an optimality theorem from this bounded numerical proposal.

The critical old moment is concentrated near6218 active rows, with tilted
U occupancies .703 on active and .120 on inactive rows. It already avoids
most inactive-row union positions. Derived and tested an alternative
forward support-transfer count with conserved source-label mass Rs.
It is weaker on the critical class (~+1.81M atR256,s6218,B131072).
This prevents mistaking a new counting formulation for an improvement.

Derived a structural obstruction criterion for fixed route AND pairings:
use true V_t=live-pair coordinates union next-state mask, shortened RS
dimension d, outside crossings J, deaths D, selected output zeros z.
Potential recurrence equations = state variables +J+D. If d-J-D-z>0,
a nonzero message exists regardless of gate coefficients. This is a
sufficient actual obstruction, unlike a positive probabilistic upper bound.
Initial padding must be removed for this rank test. No full-size positive
gap found or searched yet. Next bounded search should optimize this gap
using actual pairings, not merely the route-only projection objective.

Added RANDOM_PAIRING_BIRTH_DIAGNOSTIC.md and source helpers. Checks pass12
birth moments,10 uniform-image transfers,3888 rank-count identities,669
toy positive-gap witnesses. Toy gates are not the accepted SL2 ensemble;
they test the distribution-independent rank argument. Serial runs<2s,
26.77GiB free at initial check, no raw files. No native benchmark,
construction/library/paper changes, commits, or promotion. Earlier directed
certificates and frozen helpers preserved. Goal active; progress this turn
is a tested negative refinement and a sharper obstruction diagnostic.

### 2026-10-26 — Exact active rows; restricted coverage to16384 swaps

Kept construction fixed. Retained exact active outer-row set S inside the
conditional cap. Outside crossings J must land on S. On compatible rows,
d_S-|J|=sum_{i in S} b_i-128|S|. After a conserved routing tilt, summing
all S factors into (r^-128*a+b)^8192. The coverage>=129 condition can be
dropped safely; selected optimized markers did not exploit it.

Directed replay covers every B<=16384 in both R192/R256 half-full suffixes,
at46% overGF2^32. Nine intervals each, fixed r/t8/4 and16/12; worst directed
log bounds -404016/-35799. Quoted joint restricted bound <=2^-113.
This doubles the preceding turnover coverage, not whole-code closure.
Varying weights, resets, other suffixes, GF128, and performance stay open.

Important negative result: seven sampled active-row counts pass atR256,
B32768, but the full row-set sum gives+54115 FLOAT bits. Do not infer
coverage from the samples. AtB131072, full sums remain+208414/+297971.
These failed upper bounds are not counterexamples. At very large turnover,
forced output weight makes the selected screens pass again.

Added RANDOM_PAIRING_ACTIVE_ROWS.md, directed replay and exact tests;
extended existing capped diagnostic with --active and --active-sum.
Exact tests pass175 dimension,12 row-moment,12 conserved-route checks.
Initial512-wide replay exceeded25s after completingR192;2048-wide directed
replay completed both profiles in7.80s. Tests0.008s, FLOAT screens<1s each.
Local free space26.78GiB; serial guarded runs, no data files. No library,
paper, native timing, commits, or promotion. Prior controls preserved.

Next retain birth/death nonzero-input requirements jointly with S. Current
count drops these and allows unsupported activations on inactive rows.
Test one remaining R256 intermediate-turnover profile before widening
coverage or redesigning. Goal active; progress is restricted proof coverage.

### 2026-10-26 — Cap shared pairing events; extend to8192 swaps

Previous turn was progress: conserved-route count and256-swap certificate.
Revalidated sources and local space (23.97GiB initially). Derived a different
count: condition on route/pairing geometry, union messages, cap probability
at one, then average geometry. Many messages can share the same rare
crossing event; first-moment multiplication was unnecessarily amplifying it.

For fixed support history, conditional bound is min(1,2*C*q^(d-|J|-B-z)).
Here J are crossing endpoints outside consecutive-state unions, d is the
outer dimension shortened to U union J, B counts deaths, z selected output
zeros. Identity d-|J|=|U|-128M+sum_i(128-covered_i)_+ converts it to a deficit
moment. Actual endpoint sets are pointwise dominated by fair coins with
rho=2^256*m_256, ~4.82645 bits per update. Fractional moments and exact-column
Bernoulli conditioning avoid a simultaneous bad-route union over histories.

Directed replay now covers all half-full histories with<=8192 swaps in
R192 orR256 suffixes. Bound for joint family at46% overGF2^32 <=2^-119;
actual directed log contributions -182846/-57995. Includes all masks,
placements and messages; not varying weights, resets, whole-code, or GF128.
Earlier certificates preserved. FLOAT B4096 now -244657/-137794; B16384
also passes selected checks, while B32768 remains+41208/+106989.

Local space dropped to4.15GiB before the run; local guard correctly stopped
Python. No benchmark/compiler/Python process was found running. Read-only
inspection found allocated pagefile45413MiB, not a verified cause of the
space loss. No deletions/system changes. Used connect-to-peach skill for
the approved host: /tmp63GiB free, no numerical job underprindal. SciPy was
absent; implemented standalone standard-library helpers. Streamed source
via host-key-pinned Plink stdin, no remote files/packages. All jobs serial,
30s timeout,256MiB address-space cap,Python-B,remote20GiB disk check.

New diagnose_random_pairing_capped.py, certify_random_pairing_capped.py,
test_random_pairing_capped.py and RANDOM_PAIRING_CAPPED.md. Exact tests pass
652 endpoint,280 capped-average,295 deficit-identity,30 route-coverage cases.
Screen0.10s,replay0.28s,tests0.10s. No raw outputs, native benchmarks,
production/paper changes, commits or promotion. Goal active, not blocked:
remote execution provides a safe path despite local disk pressure.

Next retain exact active outer-row sets inside the cap: J must lie in S,
and each active row needs>=129 covered positions. Test intermediate turnover
before expanding the support-weight range or resetting histories.

### 2026-10-26 — Space restored; changing-mask certificate and route diagnosis

Previous turn was no progress after revalidating the disk blocker. Current
read-only check found25.33GiB free, so the blocker is cleared. Ran the
pending screen serially under the20GiB guard; no raw files were generated.

The fixed-weight turnover count passes low B but fails at moderate B.
R256,w256,B4096 gives+4700129 FLOAT bits; most loss is the simultaneous
route cap8.165M bits. Derived an occupancy-conserving moment using an
independent Bernoulli surrogate conditioned on exact regional sizes.
Its scalar tilt is handled before conditioning, not through invalid
negative association of nonmonotone functions. Convex smoothing handles
unequal regional occupancies. The same case improves to+932193 with a
4.397M-bit route cap but still fails. R192,B4096 improves+381995 to+236648.
At R256,B450000, old+3897867 becomes-5728231. These are FLOAT screens.

Directed replay now covers ALL histories with <=256 total support swaps
at constant weight256, for R192 orR256. Union masks are padded to equal
total area, history counts use a rational positive marker, and route
moments use rational tilts. Good-route bounds -306404 and-817359; route
failure <=2^-120 per family after union over histories. Each restricted
failure <=2^-119; combined <=2^-118 at46% over GF2^32. No changing-weight,
reset, whole-code, or extension-field closure is claimed. Earlier controls
are unchanged. A first r9/8 replay proposal failed; r3/2 passed.

Fixed-seed route-adapted search (no inner maps) found histories with swaps
378/1919/3944/16256/32552 and shortened dimensions32216/45455/55269/86289/
117894. Their actual projection-envelope bounds remain strongly negative.
This does not bound the worst history or prove absence of an obstruction;
it indicates slack in the current uniform route event. No codeword found
or claimed. Added source-only diagnostic and compact note.

Exact tests pass12044 support histories,12 coefficient identities,405
convex smoothing inequalities,48 conserved-conditioning cases. Screen
0.06s,conserved screen0.11s,directed replay0.72s,tests0.10s,route search
4.20s (serial). New certify_random_pairing_turnover.py and
test_random_pairing_turnover.py; updated notes/index using CWC. No native
benchmarks, production/paper edits, commits, or promotion. Goal active.

Next target the joint route/history count at B4096 before more tail-weight
tuning. Then varying weights and repeated activation/reset boundaries.
The moderate-turnover gap is currently a proof-summary obstruction, not
an exhibited construction failure or a request for a redesign yet.

### 2026-10-26 — Disk blocker revalidated, second consecutive occurrence

Read-only check now reports1.19GiB free on C:, still below20GiB. No live
Python/benchmark/frontend/compiler process was found. No numerical job
started, no data deleted, and no new proof result is claimed this turn.
The preceding turn made analytical progress but encountered this same
blocker. The pending discriminating screen remains unrun; restore space
or obtain a scoped cleanup direction before running it. Goal stays active.

### 2026-10-26 — Composed turnover bound; disk guard stops numerical work

Previous turn was progress: two directed confined-support certificates
and exact local changing-support checks. Revalidated those sources.
This turn derived a composed fixed-weight changing-mask bound. Histories
are counted by C(512,w)[x^B](sum_b C(w,b)C(512-w,b)x^b)^(16R-1).
Union area Lw+B, intersection area (L-1)w-B, deaths B and required output
zeros are retained jointly. Time-varying regional occupancies use a convex
binomial envelope, with history entropy included in the simultaneous route
event. This is an analytical bound, not new numerical evidence.

Added diagnose_random_pairing_turnover.py and RANDOM_PAIRING_TURNOVER.md.
The intended serial FLOAT screen was stopped before Python launched by
the20GiB disk guard. Read-only check found C: free1.49GiB; no Python,
benchmark/frontend, or compiler process was found running. No restart,
cleanup of unrelated data, raw output, native benchmark or numerical
result. The new helper is untested. Existing certificates are unchanged.

First occurrence of this external blocker; goal remains active, not marked
blocked. User informed. Next recheck free space; once>=20GiB, run the
bounded turnover screen and compare B0 with the preserved control. Then
diagnose high-turnover profiles before adding varying weights and resets.

### 2026-10-26 — RS projection closes two random-pairing families

Previous turn was progress: exact joint transfer and toy composition.
This turn tested actual RS conditioning. The initial constant-support
pin bound fails at R192,w256,s4096 (+6041139.9 FLOAT bits). Replaced that
summary, not the construction, with a forward outside-projection count.

For a confined state, crossing endpoints must equal the nonzero outside
input positions. Their exact probability is P_w(k)/C(512-w,k), whose
maximum m_w gives a per-step matching cost. RS punctured enumerators
sum all messages together: sum_x q^-outside_weight <= product_i
(1+2^(R-4096+31*b_i)). Independent regional uniform subsets give negatively
associated row counts; a simultaneous route event bounds this product.
Output zeros add independent coordinate constraints on live pairs.

Directed integer/rational replay certifies R192,w256 and R256,w256,
each including all masks and messages, at46% over GF2^32. State is zero
before the region-aligned suffix, then has one constant support through
the free terminal state. Good-route log bounds -638009 and-434936;
route failure after mask union <=2^-97 each. Each restricted failure is
<=2^-96; joint two-family contribution <=2^-95. Not a whole-code cert,
not a reset-family cert, and not a GF128 extension claim.

FLOAT neighbors: all1024 weights for R192/256 pass; worst -27384/-36369,
both w2. Additional30 selected profiles across R129/160/192/224/256 pass.
Only the two half-support profiles have directed replay so far.

Derived changing-support local bound: b births release at most a factor
2^b in the exact endpoint probability; d deaths add d field equations.
Outside union(H,H_next) input weight and selected internal output zeros
are charged jointly. Exact GF8 tests pass64 RS enumerators,12 regional
mgf comparisons,540 confined and3510 changing-support constraints with
actual accepted SL2 maps. This local lemma is not yet composed globally.

Added diagnose_random_pairing_confined.py (retains failed pin screen),
diagnose_random_pairing_projection.py, certify_random_pairing_confined.py,
test_random_pairing_projection.py and RANDOM_PAIRING_PROJECTION.md.
CWC used to distinguish actual setup, comparison bounds, and coverage.
Serial guarded runs: selected screen0.03s,1024-neighbor screen1.34s,
directed replay0.49s,exact tests4.11s. No raw files or native benchmarks.
Existing proofs/library/paper untouched; no commits or promotion.

Next compose changing-support counts with the actual RS outside projection.
Pay history multiplicities jointly with births/deaths and output zeros;
extend the simultaneous route event to time-varying union masks. Then
handle repeated activations and reset boundaries. Goal remains active.

### 2026-10-26 — Joint matching/pin transfer and toy outer validation

Construction unchanged. Extended the temporal pin identity to arbitrary
changing matchings: pins start distinct zero-output path components.
Derived a weight-indexed positive transfer polynomial jointly retaining
live coordinates, next-support births, fixed-zero inputs and output zeros.
It counts relaxed amplitudes, not physical transition probabilities.
Documented the RS interpolation interface and required exact-region
conditioning; no independence claim for RS evaluations is made.

Exact checks pass 30 polynomial rows against literal matching/support
enumeration and 128943 temporal pin histories. A separate GF4 two-step
repetition-outer experiment checks every output-weight coefficient for
one and two active rows (3456 and114048 enumerated images). At weight<=3,
exact expected counts69/32 and83/128 are bounded by115/16 and617/48.
This uses uniform pair images as a comparison, not actual accepted SL2.

Diagnostic n512,w256,q2^32,p1/2,theta0: without inactive births the weighted
crossing mean remains~1.42e-14; with all256 inactive coordinates born it
becomes~171.001 and zero-cross fraction~8.30e-123. Retaining pins changes
which matchings dominate. Confined-support histories remain unresolved.

Added diagnose_random_pairing_pins.py, diagnose_random_pairing_outer_toy.py,
and RANDOM_PAIRING_PINS.md; indexed continuation. Tests/report ran serially
in about0.94/0.08/0.29s under the20GiB disk guard (38.73GiB free). No raw
dumps, kernel changes, native timings, paper edits, commits or promotion.
CWC kept comparison counts separate from actual probabilities. Old proofs
remain intact. Goal active; no new restricted or whole-code certificate.

Next: evaluate actual RS shortening and region conditioning for confined
support and activation/reset histories. Account for outer dimension and
exclude the all-zero message before extrapolating the local transfer.

### 2026-10-26 — New goal: independently random state matchings

User requested a goal to explore random pairings. Created an active goal
covering exact law, whole-code proof integration, discriminating tests and
eventual implementation cost. Each update gets an independent uniform
perfect matching at setup; reuse it for every message. Outer,regional
route,512-state,accepted SL2 arithmetic remain fixed controls. Old cyclic
eight-defect certificate does not transfer to this new distribution.

Derived exact crossing law and conditional uniform external partner set.
Derived hypergeometric overlap with marked inactive coordinates. At w256,
mean live mass384.2505,zero-cross probability~2^-255.5007. At w2,zero-cross
probability1/511. Tests exhaust11463 matching/support cases,160 marked
profiles,and256 affine parity probability checks; all pass.

Key proof interface: fix output y first; Q_next=delta Q+y makes every
state independent of matching setup. Conditional on matchings and route,
inverse inputs are fixed on empty pairs and pointwise bounded by uniform
on live-pair coordinates with total density<2. Probability of outer-code
membership<=2 q^-rank(H_F) times the affine compatibility indicator.
Sum over low-weight outputs and average actual matching/route jointly.
This addresses the post-setup message quantifier without a false uniform
good-matching claim for every adaptively selected support.

Discriminating negative result: E[q^-cross] at q2^32,w256 has only255.5007
bits of extra inverse exponent,not~4096 from average128 crossings. The
no-cross term dominates (>0.99999999999999). Hence average spreading alone
does not close proof. Must retain next-support/outer input compatibility.
Typical random matchings also make the old rare-defective-batch partition
inappropriate; singleton pairs are common and useful spreading events.

Added diagnose_random_pairing.py and RANDOM_PAIRING.md; README indexed.
Serial bounded Python-B tests~0.05s,diagnostic~0.15s,no data files. Disk
38.75GiB before runs,guard>=20GiB. No implementation,benchmark,paper,
library,commit,push,promotion,or frozen proof changes. Prior controls kept.

Next retain the pin identity, crossing law, and marked next-state births
in a joint transfer bound. Test half-support persistence and a complete
activation-to-reset episode with actual RS support constraints before
native implementation work. Goal active,local identities only,no cert.

### 2026-10-23 — Stop at the requested redesign checkpoint

Previous turn made structural progress. This pass tested separated larger
budgets instead of extending h one unit at a time. At T1967,a1, cases
(h,c)=(8,8),(12,12),(16,16),(12,6),(16,8),(16,4) distinguish repeated
activations from many internal defects. The h12c12 bound fails at+101404;
h16c16 at+244421 after retaining the reset bound even when the localized
cap exceeds512. An ideal dense-image diagnostic is comfortably negative.

Root cause in the present summary: pair mixing permits eight clean proper
states, a tight support-geometry bound. At h12c12T1967 the105 proper-state
allowance leaves random-coordinate floor953554 below W964689. Even
retaining clean-run minimum weights yields959976<W. Routing rank tuning
alone cannot fix that forward-volume summary; more joint support
information or faster local spreading is required. This is a proof
obstruction, not an actual bad word or a proof the current code is weak.

Compared independent uniform superregular4x4 and8x8 local maps with the
same outer, route, and recurrence. Explicit rejection bounds give global
density losses<=2 and<=2^10 respectively. Their proper clean-run bounds
are4 and2. At h12c12T1967 they give-264685 and-631118; at h16c16T1967
they give+85505 and-435116. A ten-profile window over T1936,1952,1967,
1984,2000 at h12/h16 has1 pair pass,6 uniform4 passes,10 uniform8 passes.
Selected profiles only,not a new code certificate. More mixing does not
monotonically transfer the old result: h12T1936 is excluded by pair growth
but the current uniform4 bound fails because it can grow faster.

Recommendation: discuss a second crossing pair layer inside each update.
Do not assume it realizes the uniform4 control. Exact accepted-pair tests
show input(1,1,0,0) produces only two nonzeros after two crossing layers
with probability2/7 in GF8 and2/15 in GF16; input(1,0,0,0) always spreads
to four. This identifies the cancellation event that its proof must count.
Mixing-only multiplication counts are1.5/current,3/two-layer,4/dense4,
8/dense8 per coordinate; no new performance measurement or claim.

Added diagnose_pair_scaling.py and UNIFORM_PAIR_REDESIGN_CHECKPOINT.md;
README indexed. Tests pass8704 closure cases,20212 clean-run-mass
transitions,and exact small-field cancellation counts. Final ten-profile
replay31.14s under45s cap. Runs serial,Python-B,guarded>=20GiB free,
no numerical files. Existing proofs,control,encoder,outer,route unchanged.
No commit,push,promotion,benchmark,or paper edit.

Decision: the exploratory goal reaches its explicitly requested stop
condition: enough evidence for an incremental redesign discussion.
This does not complete the mathematical whole-code proof. Stop this
proof push and report; no redesigned encoder is authorized or implemented
by this checkpoint. Next user-directed task should establish the actual
two-layer local law and cost before attempting a new global certificate.

### 2026-10-23 — Growth precursors close the full eight-defect family

Previous turn made progress; construction unchanged. On the per-input cap
S<512, activation is proper and each entry into full support costs an
internal defect. Full-state runs<=e=h-c+a, rather than h+1. Also proper
incoming states<=c*floor(log2 S)+9e. Growth to full uses g defects, leaving
at most e-g for the last-full-to-zero decay. These bounds reduced the341
old full-reset residuals to five (versus69 from the prior refinement).

Exact-binomial output volume resolves four of those five. The last is
h8,a1,c7,T1967,S133,e2. Before its first full state, backwards growth
requires512->190->29, supplying380+58=438 random inverse coordinates in
two proper batches. Their union with selected full batches touches at
most124 regions, so RS parity columns remain independent. Rank rises
972800->973238; exact exponent+12851->-1165. No construction change or
independent-state assumption is used.

New diagnose_pair_budgeted_reset.py --integrate recomputes the complete
old full-reset cover:341 old residuals,341 covered,zero remaining;72
retain the old decay bound,four use exact binomial volume. Good-route
maximum-108 (a0c5T1963). Three route exceptions at-80 give the full-reset
family<=2^-78. With the preserved no-full-reset result, all histories
with at most eight defects close at46%,restricted failure<=2^-77.
Higher defects, whole-code certification, and extension-field closure
remain open. This is not a claim about the construction's actual distance.

Tests pass31239 support transitions and526 first-full precursor cases.
Selected check and complete replay were serial,Python-B,source-only.
Successful replay46.15s under60s cap. An earlier exact-binomial replay
timed out and was not counted; avoiding exact work on already-closed
old-window cases resolved it. Frozen diagnose_pair_defects.py hash checked.
Added UNIFORM_PAIR_BUDGETED_RESET.md and README index. No benchmark,
library/paper edit,delete,commit,push,promotion. Control preserved.

Disk briefly fell to14.49GiB after an earlier35.32GiB check. Held the
large replay; Peach had62GiB free but its old transfer client could not
load the configured key, so no source transfer occurred. Local space
recovered to32.35GiB; guarded local jobs require>=20GiB before starting.
No numerical output files were created.

Decision: meaningful structural proof progress; no redesign requirement
established. Next test the combined full-run/precursor bound on separated
larger defect budgets, rather than extending h one unit at a time. Locate
whether the next loss is activation, full-run fragmentation, or omitted
random coordinates in proper states. Retain actual pair support before
concluding an inner change is needed. Goal remains active.

### 2026-10-23 — Coupled terminal growth closes all h8 no-full-reset histories

Previous turn was progress. Tested growth-cost coupling first: requiring
u+g<=e does not improve worst a1c7T1955(+14000),because its old maximizing
cell already hasu0. This discriminates the actual loss: global local-input
cap ignores which batches can feed growth in the final interval.

Added restricted regional moments with terminal complete-region lower
bound and touched-region upper bound. Terminal-specific caps resolve28
of45 remaining h8 no-full-reset profiles;17 remain,all a1c6,e3. Split the
initial terminal region from complete terminal regions and retain order:
initial weight<=P,k<=r-1 defects there,then complete-region inputs<=U.
State upper 2^(e-k)(2^(k+1)-1)P+(2^(e-k)-1)U,using e=h-c+a-u.
The coarse endpoint union still permits516 at worst T1955(P33,U36,r15).

Retained fixed terminal end and initial start when a1. Endpoint union now
13*(2c-a-1)+ceil_log(2(h+1))+12+1+80; route exception still2^-80.
All45 profiles pass,worst new exponent-1909 at a1c6T1959. Complete replay
checks69551 profiles:706 old residuals,661 old localized rescues,45 new,
zero remaining. Overall good-route-230; three route exceptions-80 give
complete h<=8 no-full-reset family at46% <=2^-78 (48.3s).

Replayed preserved h8 full-reset refinement (41.1s):341 old residuals,
69 remain. Groups:a0c6(count5,T1965..1969),a0c7(2,T1967..1968),
a1c5(5,T1983..1987),a1c6(14,T1970..1986),a1c7(24,T1959..1985),
a1c8(19,T1963..1984),with gaps possible. Worst+37396. These are failed
bounds,not bad words. Complete combined certificate remains h<=7 at-77.

Added diagnose_pair_growth_budget.py,UNIFORM_PAIR_GROWTH_BUDGET.md.
Tests:3071 terminal geometry,1375 ordered growth,375 restricted one-row
moments pass; joint-row proof reuses the positive-moment argument. Runs
serial,Python-B,bounded caches,30/45/60s caps,no result files. Disk checks
21.23/23.15/35.33GiB before numerics. No construction/library/paper edit,
benchmark,delete,commit,push,promotion. Old helpers/control preserved.

Decision: concrete progress; continue proof. Next full-reset h8: retain
the actual interval that reaches full support and later resets, tighten
its endpoint/region bounds, and couple its growth/decay defect budget to
the other intervals. Do not reuse the no-full-reset terminal anchor there.
Higher defects and whole-code/extension-field closure remain open. No
redesign requirement established. Goal active.

### 2026-10-23 — Local occupancy closes all h<=7 histories; h8 limitation mapped

Previous turn was concrete progress. Preserved construction and old proofs.
Mapped all25 h6 full-reset residuals. Last-full to first-zero window has
length<=9(e+1),e=h-c+a. Replace L^512 death-time union by L*[9(e+1)]^512;
this closes16 cells. Exact initial/terminal alignment bounds partial
regions by2c+1-a and uses mass(T+c+1)/16 without eight-batch padding.
c6T1968 dimension609->358 closes. Remaining eight c7 cells have e0:
full state requires full activation, but eligible row cap<512. Complete
h6 full-reset good-route-229,family<=2^-78; combined h6<=2^-68.

Then generalized route count to eligible rows whose coordinate lies in
one chosen batch. Its joint moment is(1/16)^s times the other regions'
128-hit moment. Union over interval endpoints and chosen batch leaves
route exception<=2^-80. Actual input batch weight<=S on this event.
Clean updates cannot grow support; defective updates satisfy w_next<=2w+S.
After e internal defects,w<=(2^(e+1)-1)S. This is a deterministic necessary
condition, not a uniform-state assumption.

At h7, window/alignment leaves10 of132 full-reset residuals. Local growth
excludes all10; hardest c7T1968 hasS156,e1,statecap468<512. Full-reset
good-route-263,total<=2^-78 (30.4s). New no-full-reset cover checks61376
cells,excludes all379 old residuals,remaining good-route-310,total<=2^-78
(23.6s). Thus all histories with<=7 defective incoming batches at46%
close at<=2^-77. Still not whole-code or extension-field closure.

Bounded h8 no-full-reset diagnostic:69551 cells,706 old residuals,661
rescued,45 remain (39.6s). Groups:a0c5T1950..1951(e3,Smax36),
a1c5T1940..1960(e4,Smax39),a1c6T1943..1959(e3,Smax59),
a1c7T1955..1959(e2,Smax87). Worst old exponent+14000. First example
S35,e3 gives525>512,so the growth exclusion alone no longer applies.

Added diagnose_pair_full_reset_frontier.py,diagnose_pair_localized_frontier.py,
UNIFORM_PAIR_LOCALIZED_GROWTH.md. Tests:54241 decay,2047 anchored geometry,
400 exact localized-row moment checks pass. Serial Python-B, bounded caches,
45s numerical/20s test caps,no numerical files. Disk checks20.73/21.54/
22.78/22.93GiB before jobs. No benchmark,construction/library/paper edit,
delete,commit,push,promotion. Certified control preserved.

Decision: broader proof progress, not a redesign obstruction. Next couple
the defects required to grow to512 (from S) with the early closed-episode
budget. The separate bounds currently spend those defects twice in their
worst-case envelopes. Test the remaining h8 classes before further tuning
or redesign. Full-reset h8 and higher histories remain untested here.
Goal active.

### 2026-10-23 — Single internal defect closes; all h<=6 no-full-reset covered

Previous turn was progress. Kept the construction and all old helpers
fixed. Derived a single-batch history count: births occur only at j,
with coordinate choices binom(512,births) and death times in the remaining
r steps. Proper clean prefix/suffix force w>=2^j and w+e>=2^r, j,r<=8.
Sum position/charge bands while retaining the envelope q^(-w-e) tradeoff.
Vandermonde and a denominator endpoint bound give a cheap two-fragment
factor. Actual first alignment retained, whether first episode is clean
or exceptional. Three-fragment cases all excluded by the eligible-row cap.

Added diagnose_pair_single_defect.py,UNIFORM_PAIR_SINGLE_DEFECT.md.
Tests:24232 support histories,6588 cumulative history caps,1006 denominator minima,
30180 Vandermonde allocations pass. T1941 all lengths:6162 count cells,
4256 excluded,worst new good-route exponent-103. Integrated all34 residual
T1908..1941:64128 exceptional count cells,worst-103; clean-early subclass
replays-74. Old otherT-6244,old route-80,two new route exceptions-100.
Complete h<=6,a1,c6 no-full-reset branch at46% <=2^-71 (24.7s locally).

Fresh audit of ten remaining branches a0/1,c1..5:zero residuals,worst
good-route-1267 at a1c4 (7.1s). Together with preserved a0c6/a1c7 bounds
at-78, all h<=6 no-full-reset histories close at46% <=2^-69. This is not
a whole-code certificate. Replayed frozen h6 full-reset cover:25 residuals,
worst+9147 at a1c7T1953,first a1c5T1969 (18.9s). No redesign conclusion.

Numerics serial,Python-B,bounded caches,20/30/40/60s internal caps; no data
files or cache output. Started22.86GiB free. No benchmark,construction,
library,paper,delete,commit,push,promotion. Certified control preserved.
Frozen defect helper SHA256 remains7238fec9...54bbaed3. Disk recheck
before final tests remained above20GiB; unrelated disk usage increased.

Decision: concrete progress; continue proof. Next map the25 full-reset
residuals and retain defect budget in the full-to-zero decay window and
route cost of full activation. Old reset cost allows arbitrary whole-stream
death times per coordinate; test whether that is the lost information.
Then higher defects and extension-field closure. Goal active.

### 2026-10-23 — Unaligned branch closes; internal defect isolated

Previous turn was progress, closing h6/a1/c7. Kept that proof and source
unchanged. Added diagnose_pair_unaligned_frontier.py and
UNIFORM_PAIR_UNALIGNED_FRONTIER.md, using five clean early episodes,
arbitrary single/crossing placement, actual RS zeros, and shared-region
correction. Six interval endpoints cost156bits; cap requests256bits.

Fresh old-verifier replay for h6/a0/c6:40 residualT=1893..1932,
worst+3003 at1931,remaining good-route upper-6990. Integrated all40:
5711 count cells,new worst-139 at1925. Old route exception-80,new-100;
complete no-full-reset branch at46% <=2^-78. Not wholecode.

Inspected neighbors from source. a0/c5 has no residuals,worst-76454.
a1/c6 has34 residualT=1908..1941,worst+4091 at1940,remaining-6244.
For the latter's clean-early subclass, the one internal defective batch
may be in the unrestricted terminal interval. Reuse the five-factor
count, conservatively ignoring forced first alignment. All34 close:
9822 count cells,worst new good-route exponent-74 at1941,combined<=2^-72.
This is a subclass, not a complete a1/c6 branch certificate.

576 explicit five-episode length-envelope checks pass. The preserved
a1/c7 integration replays unchanged at-78 (13.1s locally); frozen defect
helper SHA256 remains7238fec9...54bbaed3. Numerical jobs
serial,Python-B,bounded caches,40s internal caps; integrations about8s
locally. Started with22.86GiB free; no dumps or caches. An initial inline
command failed Python quoting before executing; replaced by retained source.
No construction/library/paper edit,benchmark,delete,commit,push,promotion.

Decision: concrete proof progress; no redesign requirement established.
Next address a1/c6 with the single internal defective batch in one closed
early episode. Four other closed episodes are clean. The exceptional
episode can have length<=17 and span three regions. Use its single-batch
defect budget, rather than allowing births everywhere. Test T1940/1941
first, then integrate. Full-reset/general defect-rich/extension-field
closure remains open. Goal active.

### 2026-10-23 — Complete aligned h6/a1/c7 no-full-reset branch closes

Previous turn was progress: all-charge selected geometry plus coverage audit.
Acted on that audit, targeting the retained worst global profile rather
than another isolated charge slice. Kept the construction fixed.

Added the one-region activation/reset factor, valid for the aligned first
episode and noncrossing episodes. For h6,a1,c7, six later activations consume
all defects; six early episodes are clean/proper with lengths1..8. Enumerate
maximum length d,total early length E,crossing count p<=5,terminal-region
count,and duplicate-region count g. Coefficientwise local length envelopes
relax only the local product; actualE still controls routing. Generalized
the finite-population correction to any number of fragments sharing a region.
Retain actual RS zeros and virtual deficit Z+g*s.

The T1938 screen covers1048 count classes,656 excluded by eligible caps;
new good-route exponent -127. Recomputed the old branch from source: exactly
63 residualT=1878..1940,worst+6491 at1938. Integrated all63,23763 count cells;
worst new good-route exponent -127,first at1926,new route exception -100.
Old verifier covers every otherT with good-route exponent<=-3806 and route
exception -80. Combined complete h6,a1,c7 no-full-reset branch at46% <=2^-78.
The new exclusion itself uses no output-weight bound. Not allh6 or wholecode.

Added diagnose_pair_aligned_frontier.py,UNIFORM_PAIR_ALIGNED_FRONTIER.md.
66 exact one-region image checks,4746 multifragment ratios,105 allocation
checks pass. Integrated replay4.7s on Peach under256MiB/45s limits; local
test/screen began>22GiB free. Jobs serial,Python-B,only source/compact notes.
Old defect helper hash unchanged. No benchmark,production/paper edit,delete,
commit,push,promotion. Certified control and earlier proof checkpoints preserved.

Decision: progress; continue proof. Next target the actual a0,c6 clean
activation-maximal branch with five early episodes and no forced alignment;
recompute its residuals before integration. Then use all-charge envelopes
for branches with internal defective batches. Full-reset branch and general
defect-rich/extension-field closure remain open. No structural redesign
requirement established. Goal active.

### 2026-10-23 — All charge bands join sharing; audit identifies geometry gap

Previous turn was progress: defect-envelope identity and isolated e1/e2
profiles. This pass combined that interface with shared regional routing
and removed the small-coordinate-charge restriction for one closed episode.

For charge band[a,d], use probability q^(-w-a), envelope<=w+d, initial-mask
count at d, and history count T(l,w,d). The history count covers the whole
band; no extra charge-cardinality union is needed. Cap it by2^(b(l-1)), the
number of all intermediate support masks for fixed endpoints. Retain actual
RS zero count and exact shared-region correction from prior work.

Five proper crossing episodes of length3,R131,one shared region,four clean
and one arbitrary closed episode: allsource counts8..48 included. Bands
1--2,3--8,9--32,33--128,129--1536 give good-route exponents -503,-447,-352,
-480,-3292. Clean case is -623. Contiguous replay verifies full charge
coverage; six-way union<=2^-349, plus the common route exception2^-100 once,
gives restricted family<=2^-99. No output-weight restriction. Arbitrary
internal activity in one episode is covered, not arbitrary whole-code history.

Coverage audit against actual h6 partition: a0,c6 has five early episodes
but noncrossing/multiple-sharing cases remain; a1,c7 has six early episodes
and an aligned first activation. The retained worst(a1,c7,T1938) has no
internal defect budget, so new charge refinement alone cannot close it.
One-internal-batch branches(a0,c5)/(a1,c6) have different episode counts or
alignment and variable lengths. No claim of global residual reduction.

Added diagnose_pair_defect_bands.py,UNIFORM_PAIR_DEFECT_BANDS.md.72 band
inequality tests pass; complete six-class deterministic replay5.3s on Peach,
256MiB limit,40s timeout,35s per-profile cap. Local C: below20GiB so no local
numeric jobs; remote/tmp62GiB free. Serial runs,Python-B,source-only transfers
to existing temp directory; no raw files,caches,benchmarks,deletions,commits,
pushes,library/paper edits. Old construction and certified control preserved.

Decision: progress and continue proof. Next implement the single-region/
aligned activation factor, then target a1,c7,T1938 with correct geometry
coverage. This reconnects the new interface to the real frontier instead
of accumulating isolated examples. General defect-rich and extension-field
closure remain open. No structural redesign requirement established. Goal active.

### 2026-10-23 — General defect envelope and first internal-defect profiles

Previous turn was progress: one shared region closed for short clean episodes.
This pass targeted the next required proof interface, rather than another
clean-length sweep. Construction, old proofs, and certified control preserved.

Derived an exact closed-episode identity. With initial weight w, H births
from empty pairs, P one-sided pair occurrences, and e=H+P, the number C of
zero next-state coordinates in live pairs is C=w+e. Proof: deaths=w+allbirths;
live-pair births consume exactly the corresponding singleton zero constraint.
All inputs lie in the union of initial support, empty-pair birth coordinates,
and singleton partners, whose size is<=C. This replaces the false assumption
that defective inputs remain in initial support; no state-uniformity claim.

For fixed input/support history, probability<=q^(-w-e) in the uniform-image
comparison, with one global accepted-map density factor two. Fixed-U history
count<=l^w sum_(B<=e)binom(bl,B)l^B via birth labels and component death times.
Retained initial singleton-count mask enumerator and enlarged envelope in
the actual regional subset count. This supplies a joint defective-episode
route bound compatible with the previous actual-RS-zero DP.

Selected five length-three crossing episodes, disjoint regions,R131, live
tail: four clean, one of exact charge1 or2. All21 source counts8..28 and
allfive choices of defective episode included. Good-route exponents -747
and -750; route exception -100; each restricted family<=2^-99. Clean control
with the same actual-zero count is -813 before route exception. No output
weight restriction. These are separate families, not whole-code closure.

Added diagnose_pair_defect_envelope.py and UNIFORM_PAIR_DEFECT_ENVELOPE.md.
115893 exhaustive support-history identities,245 history-count bounds,66
initial-mask counts pass. Profile runs<2s each, tests~1s. Local space checked
separately before runs (>23GiB); serial jobs, bounded caches,20/35s source
caps,Python-B, no raw output/cache. No old helper,library,paper edit,benchmark,
deletion,commit,push,or promotion.

Decision: progress; continue proof. The inner already charges defects that
enlarge input support. Next combine this interface with shared regions and
test a retained internal-defect residual. Map coordinate charge e back to
the certificate partition carefully: e is not the defective-batch count,
which can hide many coordinate events. Full defect-rich and extension-field
closure remain open. No structural reason for redesign yet. Goal active.

### 2026-10-23 — One shared region closes after retaining actual RS zeros

Previous turn was progress: complete disjoint clean-crossing family at131
regions. This pass tested dependence when two intervals share a region.
Construction, old helpers, and certified control remain unchanged.

Derived a virtual-region comparison with explicit correction. For s active
rows, each shared region costs at most(s+1)*2^s*(M/(M-s+1))^s: denominator
ratio for finite-population sampling plus a sum over random allocations.
The actual message dimension remains k=R-128, not the virtual k+g.
Eligible-row factorial moments are bounded by splitting weights only after
the real-region same-set moment bound. No independence of actual fragments.

The relaxed count closed the selected all-length-three profile but failed
the56-profile length1--4 family, worst+694 at five length-four episodes.
Root cause: early deficits were free up to their maximum, while every
message was still charged q^(ks). Retained actual RS zero count Z instead:
message count <=binom(Rs,Z)q^(ks-Z), virtual early deficit<=Z+g*s.
Joint route/reset DP now tracks that budget before summing Z.

For five clean proper crossing episodes, lengths1--4, an unrestricted live
tail,131 actual regions, exactly one shared region, all56 length multisets
close. Good-route upper2^-106, route exception2^-100, full restricted
family2^-99. Worst five length-one episodes, s2,Z0. All row counts and
endpoints included, no output-weight restriction. Not whole-code closure.
Longer episodes, more sharing, noncrossing intervals, and internal defects
remain open; no inference of their closure from this result.

Added diagnose_pair_shared_regions.py and UNIFORM_PAIR_SHARED_REGIONS.md.
3860 finite-population ratio tests,225 exact image/allocation tests,995
RS-zero/splitting checks pass. Corrected56-profile run7.3s on Peach under
256MiB and45s limits. Earlier local/remote all-length sweeps hit40s caps
and are incomplete. Runs serial, no raw data/caches. Local disk initially
>23GiB, later about0.1GiB; final tiny test shared a shell call with the space
check and ran before that low-space result was inspected. Stop further local
numerics pending revalidation. Existing small remote source directory reused.
Old defect-helper hash unchanged. No deletion,benchmark,commit,push,promotion.

Decision: progress; continue proof. Extend the joint zero/activation count
to longer clean episodes with bounded profile groups, then births/singletons.
The latter can violate initial-support containment and need a new explicit
interface. Current evidence still points to counting losses, not a structural
reason for redesign. Goal stays active.

### 2026-10-23 — Joint activation containment closes131-region clean family

Previous goal turn was progress: extra-RS-zero coverage at130 touched
regions. Continued the fixed construction and tested the omitted activation
constraint on its hardest retained131-region profile.

Audit correction: available shortened dimension does not determine actual
word weight. The previous reset-price coupling omitted additional zeros
inside available positions. Repaired with the weighted RS enumerator:
each extra zero costs32 message bits, has<=256 locations, saves<=16 reset
bits. All extra-zero classes inflate the weighted row bound by<=256/255,
or<=ceil(s/128) bits for s rows. Replayed all792 length multisets for129
and130 regions; reported restricted bounds remain2^-99 and2^-70. No gap
is silently carried forward; note and source now distinguish both notions.

Derived a stronger joint route formula for a clean episode. If initial
support has even size w, choose w/2 matching pairs, then place every later
input into positions using those same coordinates. Multiply exact subset
counts across the two regions and weight by(l/q)^w. No separate independent
pair-closure assumption. Across five disjoint episodes, optimize counts
with the RS zero budget sum deficits<=(R-129)*s. Use q^((R-128)s) messages
before route exposure, avoiding the repaired coupling altogether.

The hardest(3,3,3,3,3) profile atR131 now has good-route bound2^-429;
with the route exception,2^-99. Extended to all792 length multisets:
707 excluded by eligible-row caps; worst is five length-one episodes;
good-route upper2^-90, route exception2^-100, full restricted family2^-89.
All active-row counts/endpoints included. No output-weight restriction.
Still not whole-code, full-six-defect, or extension-field closure.

Added diagnose_pair_activation_route.py and UNIFORM_PAIR_ACTIVATION_ROUTE.md.
80 exhaustive width-four regional-image tests match the exact weighted
formula. Selected profile<0.1s, complete family2.3s on Peach. Local C: had
only0.23GiB free, so no local numerical work. Remote/tmp had62GiB free;
jobs serial,256MiB address limits,10/45s timeout,Python-B. Existing tiny
remote task directory reused; only small source files transferred, no
raw numerical data. No deletion,benchmark,commit,push,library/paper edit.

Decision: progress and continue proof; clean-episode failures were missing
joint route constraints, not established structural defects. Next extend
this count to shared regions, retaining without-replacement dependence;
then address internal births and one-sided pairs. Do not multiply episode
route factors as independent when they share a region. Goal remains active.

### 2026-10-23 — Extra RS zero closes via eligible-row dimension counts

Previous turn was progress: a reusable clean crossing family closed at
129 touched regions. This turn kept the construction fixed and tested the
first extra-zero extension, rather than another broad defect sweep.

For fixed interval geometry, defined eligible rows by at least129 routed
coordinates inside the allowed input support. Positive factorial moments
bound threshold counts at129,130,... jointly without row independence.
The nested counts control shortened RS message dimension. More message
degrees also force fewer free RS zeros, hence stronger reset charges.
Retained this dependence and combined its good-route bound with an
unconditional route/reset first moment, without counting route rarity twice.

For five clean proper closed episodes, each crossing a regional boundary,
and an unrestricted live tail, with disjoint touched regions, closed all
792 length multisets at130 total regions. Good-route contribution <=2^-71;
route exception <=2^-100 after endpoints/thresholds; combined <=2^-70.
All active-row counts included. No output-weight restriction used. At129
regions this proof gives2^-99; preserve the stronger prior2^-142 family.
Neither result is whole-code, all-six-defect, or an extension-field closure.

At131 regions the current bound remains positive, worst for five length-three
episodes with eligible-row caps(28,18,13), reset prices(56,63,70). This is
an unresolved counting bound, not a codeword or construction counterexample.
Next test actual activation-pair closure and subsequent matching constraints
on this fixed family; the current route eligibility count ignores both.
Shared regions and internal defects remain open. No redesign justified yet.

Added diagnose_pair_eligible_rows.py and UNIFORM_PAIR_ELIGIBLE_ROWS.md.
Tests pass810 factorial majorants,512 zero-price checks,100 balancing checks,
738 nested-dimension checks. Complete three-region-count replay took26s,
with a35s cap and small bounded tables. Final source differs only by the
diagnostic column name and added tests from the measured proof replay.
Numerical runs serial; no raw outputs/caches. Revalidated >22GiB free before
work and >23GiB after. Old defect helper hash remains unchanged. No benchmark,
library/paper edit, commit, push, promotion, or deletion. Goal remains active.

### 2026-10-23 — A reusable joint reset count; clean crossing family closes

Previous turn was progress (five-defect certificate and a six-defect
restricted profile). This pass tested a reusable proof interface rather
than merely incrementing the defect cutoff. Construction remains fixed.

Derived a positive weighted RS/route/reset enumerator. For fixed intervals,
activation positions carry weight l/q, ordinary allowed positions weight1,
and forbidden positions weight0. The one-row nonzero enumerator is bounded
by sum_z q^(128-z) e_(256-z)(regional mean weights). Positive expansions
and weighted sampling without replacement control all rows jointly without
independence. Integer fragment polynomials retain early lengths, defect
costs, and region counts. The activation-only version reduced the old137
coarse six-defect no-full-reset residuals to135, but did not close them.

Identified why it misses constraints: it permits rows to avoid activation
and inject later without charging their effect on the state. In a clean
proper episode, supports decrease and are pair-closed. Every input lies
in the initial coordinate support U, with |U|>=2^l by backward closure.
If s row inputs occupy u batches in a region, then |U|>=ceil(s/u).
These deterministic facts supply a stronger joint source/reset bound.

Closed a family with m=5..64 clean proper closed episodes, each crossing
one regional boundary, followed by an unrestricted live terminal interval.
Different intervals touch disjoint regions; total touched regions=129.
Every active RS row is then nonzero in every touched region and has one
message degree of freedom. Joint reset and containment costs dominate
that count for all source sizes. Exact rational checks cover480 cases
(m, maximum early length), all interval placements and all row choices.
Worst profile upper is -151; full family union is below2^-142. There is
no output-weight restriction in this family exclusion, and no restriction
on internal births/one-sided pairs in the terminal live interval. This is
still not whole-code or complete six-defect closure.

The theorem covers the whole h6,a0,c6 no-full-reset band T1893..1908, and
joint cells whose terminal interval reaches119 regions. A completed replay
with the broad-band gate left120 coarse residuals (upper+6491); the full
reset branch was not changed. A later replay after adding the joint-cell
gate hit its50-second cap and was rejected as incomplete. Do not report
a new complete frontier count from that timed-out run.

Added diagnose_pair_reset_weights.py and UNIFORM_PAIR_RESET_WEIGHTS.md.
Tests pass625 exact GF4 weighted-RS checks,2500 weighted route inequalities,
128 fragment envelopes,69900 coordinate containments,1062 backward-closure
checks, and480 exact family cases. The completed diagnostic used about21s;
the family test used <1s. Runs were serial with bounded tables/caches and
20/50-second limits. No output files or Python caches were generated.
Prior defect and frozen generic helper hashes remain unchanged.

Local C: began above22GiB free but dropped to about0.6GiB during the final
replay. Stopped local numerical work at that safety-floor breach; no
unrelated files deleted. Future runs need the configured remote host or
revalidated free space. Small source/proof notes only were added here.

Decision: continue proof, no structural need for redesign established.
The next unresolved interface is joint RS zero-pattern and shared-region
accounting: beyond129 regions an active row can vanish in a region that
the clean crossing bound charges. Test that extension on a retained
residual before another sweep. More within-episode defects also remain
open. Goal stays active; no benchmarks, production/paper edits, commits,
pushes, or promotion.

### 2026-10-23 — Five defects close; partial resets explain a six-defect residual

Previous turn was progress: four-defect restricted certificate and exact
checks. Continued with the fixed three-shear model, not an encoder redesign.

Derived a simultaneous actual-route rank bound for unions of full batches.
Only endpoint regions are partial. Positive factorial-moment expansions
handle joint row events using the regional falling-factorial law, without
row independence. Integer/rational tail caps bound each layer of the RS
shortened dimension and its complement. Uniform endpoint/threshold unions
permit setup-adaptive interval choices. The selected h5,T2058,c6 frontier
improves from +427855 to -1402737 on its good-route event.

Coupled early interval lengths to their necessary extra defects. Length l
of a proper-only early interval needs >=ceil((l-8)_+/9) internal defects;
these cannot also enlarge the final partial tail. A joint (length,cost)
DP replaces separately maximized regional cover and proper-state count.

Complete h<=5 coverage at 46%: no-full-reset branch checks 45020 base cells
and 29991 refinements, with good-route upper -1253 and route failure -80.
Full-reset branch checks 45020 base cells and 1085 refinements, with
good-route upper -1599 and route failure -80. Both have zero residual.
Each branch totals <2^-79, and the full restricted family totals <2^-78.
This is not a whole-code certificate or a claim about defect-rich histories.

Tested h6 to assess scalability rather than extrapolating. The broad
no-full-reset replay leaves 137 cells (upper +6491); full-reset leaves 25
(upper +9147). The first no-full-reset residual h6,a0,c6,T1893 forces five
length-one early episodes, 129 occupied regions and eleven partial regions
of occupancy1/16. A one-row message would activate a forbidden one-sided
pair. For s>=2 active rows, route containment costs 2^(-44s), message count
is <=q^s, and five short resets cost q^(-5s). The exact summed bound gives
this entire profile <2^-131 including the broad endpoint union. Other h6
cells remain open; the full coverage replay intentionally still reports them.

Added diagnose_pair_region_rank.py and UNIFORM_PAIR_REGIONAL_RANK.md.
Tests: 1500 exact joint factorial checks, 18 exact directed tail-cap checks,
21 exhaustive geometry tables, 32768 sparse-reset length vectors. Fresh
h5 replays agree (6.4--8.2s each); h6 diagnostics take 12--15s. All are
serial standard-library runs with 20/30/45-second caps and bounded caches.
Local free space stays >23GiB. No raw files, caches on disk, benchmarks,
production/paper changes, commits, pushes, or promotions. Previous defect
helper and both frozen generic-helper hashes match their checkpoints.

Decision: continue proof, no redesign condition met. The next reusable
interface counts activation support in every closed partial episode jointly
with active outer rows and routing. Extend it to neighboring h6 cells, then
test whether it connects to the pin-transfer bound for defect-rich histories.
Do not spend further turns merely incrementing a defect cutoff without a
credible path to broader coverage. Goal remains active.

### 2026-10-23 — Four-defect closure with the actual route and reset cost

Previous goal turn: progress, with a proved zero-defect family and recorded
checks. Continued the same fixed three-shear construction; no redesign.

Derived P<=9(h+a)-c for proper incoming states, where a indicates initial
activity and c counts nonzero intervals. Births from zero consume defects
but cannot themselves be proper states. Full states form <=h+1 intervals.
Retaining outer parity-check rank on their union improves the all-route
count. Complete integer coverage closes h<=3 at 46% with log2 upper -161687.
This uses 28656 count cells and includes every placement and amplitude.

Separated histories with a full incoming state followed by a zero state
before the free terminal state. Without this full reset, the last interval
continues to the end; earlier intervals are short/proper. A small regional
cover DP plus forward/inverse counting closes all h<=4 no-full-reset
histories with log2 upper -3723. No independence shortcut is used for
adaptive random-coordinate counts; padding defines the comparison draws.

For general h<=4, exact Hamming volumes leave only 30 count cells. Their
allowed inputs occupy <=994816 positions in <=130 regions, with <=11
partial regions. Positive joint-row expansions and the actual regional
falling-factorial law bound the probability that any such support has
32 eligible RS rows by 2^-151, including the full interval-placement union.
On the complementary route event, every such input space has dimension<=62.
A full reset forces an active-to-zero transition on each of 512 coordinates;
time-ordered exposure costs at most 2*L^512*q^-512 after choosing times.
This closes the remaining cells. Contributions are -151 for route failure,
-293 for full resets on good routes, and -3723 without full resets.
Their sum is below 2^-149. This is the ENTIRE h<=4 restricted family,
not whole-code closure or a certificate for defect-rich histories.

Added coded_state/diagnose_pair_defects.py and UNIFORM_PAIR_DEFECTS.md.
Checks: 5934 support states/121552 transitions; 32767 interval geometries;
16384 terminal-live geometries; 20099 binomial inequalities; 2108 exact
volumes; 6912 adaptive-draw checks; 1500 exact route-tail checks; and 256
coordinate-death checks. Complete h<=4 replay checks 36839 count cells,
17707 deterministic exclusions, 34 exact-volume refinements and 30 route
rescues, with no residual cell. Runs take seconds with 20/30-second caps.
Local free space recovered to >24GiB; serial local standard-library work
created no raw data or caches. Previous proof helpers are unchanged.

Five defects remain open: the no-full-reset coarse frontier is near
T2058,c6, where retaining actual partial-region occupancies should improve
the outer rank. Next test that frontier, then reconnect this small-defect
argument to the pin-transfer method for widely distributed defects.
Decision: real proof progress, no evidence yet requiring redesign. Goal
remains active. No benchmark, library/paper change, commit, push, promotion.

### 2026-10-23 — Active goal: episode-aware proof before another redesign

Created the requested application goal. Keep the three-shear model fixed;
develop the episode-aware proof, and stop/report if evidence instead calls
for a specific incremental redesign. Prior proofs and controls are preserved.

Closed the entire no-internal-birth/no-one-sided family, including every
aligned regional suffix and all activity counts, at46%. The proof confines
states to a full prefix, at most eight partly active incoming states, and
a zero suffix. It includes the unconstrained terminal input/output. This
confines the random inverse inputs to a prefix, giving d(F)-D<=-K+24b
for every route unless the RS outer excludes the message outright.
The direct output union gives log2 family bound-194007 at46%.

An independent interval replay covers394240 count intervals across all
R129..256, with259072 deterministic exclusions; worst class bound-217567.
The formerly failing R256,V1048320 profile has sharper exact-volume bounds
-227102 at49% and-171302 at49.1%. Its49.5% bound remains positive. None
of these results is whole-code closure or a general distance limitation.

Derived the next deterministic interface: if h incoming batches contain
one-sided pairs or a following internal birth, there are at most9h+8
partly active incoming states. Clean runs have at most eight each. An
exhaustive width4/length20 support DP verifies the general d*h+d-1 form
with4140 final states and411280 transitions. Also checked972 exact volume
majorants and274 exhaustive proper-support starts.

Decision: useful proof progress; no redesign stop condition has been met.
Next retain the intervals created by one defective batch in the actual
regional route count. Then split clustered defects from widely distributed
defects and combine the latter with the pin-transfer bound. Do not assume
that uniform-feedback reset probabilities transfer to partially live states.

Added diagnose_pair_episodes.py and UNIFORM_PAIR_EPISODES.md; no previous
source changed. Local C: had only0.5GiB free, so all numerical work ran
serially on Peach with62GiB free,256MiB memory cap,30/60s timeouts, and no
raw outputs or caches. Runs finished in under two seconds. No benchmark,
library/paper edit, commit, push, or promotion. Goal remains active.

### 2026-10-23 — Incremental pair redesign removes a proof complication

User requested root-cause analysis followed by incremental inner changes.
Used the previously proposed rejected-uniform SL2 pair model; kept the
outer, regional route, state512, reset, old candidate, and certified control.
No native encoder change or new timing. See coded_state/UNIFORM_PAIR_PINS.md.

Derived the exact identity rank(L)=H for fixed inverse-input constraints:
H counts initial active coordinates and births from empty pairs. Such
constraints pin distinct temporal components at their first vertex. With
amplitude dimension a=A-Z, the auxiliary kernel has dimension a-H. The
old fixed-coordinate links and exceptional labels disappear. A two-state
support transfer therefore retains zero-output rank fully in both counting
branches, rather than tracking later pinning of unclosed components.

Independent GF4 ranks check all28561 four-coordinate, three-vertex
support/zero-output histories. Exact literal sums check64 two-update
contractions. The actual SL2 distribution/factorizations retain their
earlier exact local checks. All setup remains fixed across messages.

Closed the stronger-model profile R192,z393088,f294816,B98272 over all8192
active-row counts and all one-sided orientations. Its matching-envelope
dimension cap is64000. Four frozen proposals replay with rational inputs
and directed integer logarithms. After the source-count union, good-route
log2 bounds are -5716/-5717 at log scales2^24/2^28, both for40% and49%.
This is a restricted base-field profile certificate, not whole-code closure.

R160,z=f327552,B81888 at s4096 improves from a +168436 one-update proposal
to +23047 with two-update/no-singleton counting and to directed -218 at
m27. Only this slice is checked; no full-profile claim. R256,z=f524160,B0
still leaves about+3798111 in the analogous two-update test.

Proved why the latter example needs temporal geometry: without births or
one-sided pairs, a proper support cannot remain nonempty through all nine
distinct pairing phases. Backward matching closure would force the initial
support to contain all512 coordinates. Thus a full-state prefix is followed
by at most eight partially active incoming states and an inactive suffix.
The free terminal update remains separate. Exhaustive small-width support
tests pass at2,4,8,16. This is not a bad-codeword construction.

Next: keep the three-shear model fixed and combine that full/partial/dead
interval structure with actual regional routing. Then charge births and
one-sided pairs as defects. Only add a crossing layer or periodic stronger
refresh if a remaining partly-live family demonstrates what additional
rank/probability property is needed. Do not use more randomness to patch a
bound that discarded the placement of zeros. The whole-code GF32 control
still reaches46%; the stronger-pair whole-code certificate remains open.

Preserved all79 prior coded_state Python sources and both frozen helpers.
Added only diagnose_pair_pins.py and the research note/index entries.
Checks and bounded proposals ran serially locally with >20GiB free, small
arrays, Python -B, no raw files, no native benchmarks, and no commits/pushes.
CWC was used to separate the proved identities, restricted certificates,
failed bounds, and next construction decisions.

### 2026-10-23 — Bounded global diagnosis and ideal-inner comparisons

User requested a bounded diagnostic before deciding whether to finish the
current proof or change the design. Preserved every existing coded_state
Python source and both frozen theory helpers. No encoder or setup changes.
See coded_state/LAYERED_BOUNDED_DIAGNOSIS.md for the selected grid, derivations,
compact results, checkpoint hash, and limitations.

Screened 8480 selected aggregate tuples across ten suffix lengths; 2132 fail
necessary tests. With three bounded weight proposals and complementary
forward/inverse bounds, 4506 retained tuples pass and 1842 remain unresolved.
Each retained tuple includes all active-row counts. The gaps between grid
points are NOT covered, and retained tuples need not represent feasible
histories. No percentage of whole-domain coverage is claimed.

The frozen bank badly overstated some failures: reweighting two tested
slices changed positive bounds into -5568812 and -10091250. A dedicated
R200 central-profile probe still leaves +53315 at s8192. Larger suffix
length is not monotonically harder: the same central counts pass at R224
and R256. Stop treating one more suffix as a measure of global progress.

Replayed the complete GF32 uniform-feedback running-sum control at 46%,
at 80 and 120 digits: same 459 rectangles, 4196352 episode-count pairs,
partition hash, and conservative 42-bit margin. Same outer, route, state
size, and reset-support ranks. Thus these components are compatible with
more than 40%, without certifying the cheaper local randomness.

Derived an intermediate rejected-uniform SL(2,q) pair model. It preserves
both pair superregularity conditions and has point cap 1/(q^2-7q+2) on any
nonzero pair. Total density inflation at q=2^32 is below 4096/4089. Accepted
matrices factor into three shears, versus two currently; no native timing.
Exact GF4/8/16 enumeration checks 333 image laws and 2763 factorizations.
The simple comparison bound passes 1385 of the 1842 old failed tuples,
including the R200 central profile. The remaining 457 are mostly empty-heavy;
different old e-values can project to the same ideal class. No full pair
certificate is claimed.

Recommendation: change the inner research model incrementally to stronger
pair randomness, retaining the old candidate as a checkpoint. Next separate
zero whole-state intervals from holes within a live state in the ideal-pair
proof. Reuse the uniform control's episode accounting for the former; test
a crossing layer or periodic stronger mixing only if the latter remains
the obstruction. Do not launch another unbounded neighborhood sweep.

Validation: six exact point-replay comparisons and all three old suffix
regions pass. Local C: free space was only about 10.5 GiB, so numerical runs
were moved to Peach with >60 GiB free. Runs were serial, time bounded,
1.5 GiB process-memory capped, and used Python -B. One diagnostic timed out
due to repeated binomial calculation and was rejected; caching constants
gave a complete selected screen in about five seconds. Only source and
compact notes retained. No benchmark, paper edit, commit, push, or promotion.

### 2026-10-23 — Wider suffix 193 and actual-moment coverage at 194–195

Reoptimized proof weights under the sharper route bound, keeping every
construction component and its setup distribution fixed. A new fractional
weight closes the suffix-193 z=49500 boundary slice at -3095 and the full
z=0..49500/radius-128 region at -368. At suffix 194, continuous reweighting
still left a +1444 slice and a +4726 complete-profile upper bound.

The residual routing loss comes partly from binomial domination of matching
moments. Added the actual finite factorial-moment series, evaluated with
outward integer rounding at 96-bit scale and fixed tilt 3/2. The same
subset/suffix union gives the same 2^-56 matching failure budget. The new
route event implies the preceding event; the head event is unchanged.

Central dimension caps drop from 155878/159381/162877 to
154310/157740/161180 at R=193/194/195. Complete-profile good-route upper
bounds are -51731/-35636/-17034. Certified count regions, each unioning all
active sets, s, integer counts, histories, and subdivision boxes:

- R193: z=0..51000, independent radius 512 in the other four counts;
  12 boxes, log2 upper -398.
- R194: z=0..50000, radius 256; five boxes, upper -882.
- R195: z=0..49264, radius 128; six boxes, upper -1084.

The remaining counts are centered at (24568,589632,368520,49136). These
are GF32 history-count bounds at the exploratory 40% threshold. Full
base-field coverage, and hence the intended scalar-extension distance
guarantee, remains open. The GF32 sandwich is still the certified control.

Explicit failures: R193 z through 52000 at radius 512 retains failed leaves;
reweighting its difficult s=8024 slice gives +2698. Central retained-weight
bounds at R196/197/198 are +6294/+32620/+51482. These are failed bounds,
not codeword witnesses or evidence of optimality.

Next: apply the actual-moment envelope across the broader count-domain
cover and target uncovered boundary cells. Compare dedicated weights with
additional nonzero-input relations where refinement stalls. Do not change
the construction or substitute a sampled screen for complete coverage.
See coded_state/LAYERED_MATCHING_MOMENT_COVER.md for derivation and replay.
Validation: exact replay reproduces all three gapless regions. Passed 140
enumerated matching identities, 315 rational tails, five precision checks,
14 integer-weight equivalences, 1368 local fractional-weight checks, and
480 direct-expression majorants. Existing location, component, rational-log,
and binomial-matching regressions also pass. Both frozen theory helper hashes
are unchanged. Disk free space remains about 43.2 GiB; no raw dumps retained.
No benchmark, raw output file, production/paper edit, commit, push, or promotion.

### 2026-10-23 — Fractional crossover weights and suffix-193 coverage

Kept outer, regional routing, layered inner, and setup fixed. Warm-started
proof-weight refinement closed the first active-row crossover. At the next
crossover, the optimized integer interpolation exponents 24 and 25 still
gave +307 and +327. Generalized the checker to rational tau=2^m, retaining
exact local weights and directed global exponent bounds. A weight near
m=24.5 gives -1257 for that slice. All 243 selected neighbors now close;
the weakest complete-profile bound is -564, or -556 after their union.

The full R=192 region z=0..49264 and independent +/-128 variation in the
other four counts closes with 29 disjoint boxes and good-route upper -111.
This covers interiors, not just the selected neighbors.

Continuous interpolation alone did not close R=193: a bounded refinement
at s=8061 gave +5685. Traced the loss to the matching-edge envelope. Used
the already-proved binomial moment bound directly with four rational
Chernoff tilts, instead of only its Hoeffding/Bernstein simplifications.
This refines the route event with the same 2^-56 matching failure budget.
The head event is unchanged. The new event implies the old event.

At R=193 the central dimension cap drops from 156715 to 155878 and the
complete-profile upper becomes -10674. The full region z=0..49264, with
independent +/-128 variation in the remaining counts, closes in 19 boxes
with count-union upper -261. At R=192, z=0..50000 with +/-256 variation
closes in six boxes, upper -758. These are good-route contributions at
40%, not whole-code margins or a near-half certificate.

Unresolved: central R=194 and R=195 have retained-weight bounds +9195 and
+29245. Wider R=193 attempts (radius 256 or z through 49500 at radius 128)
still have failed leaves at a 256-test cap. No partial region is labeled
complete. The next test is reweighting under the refined route cap, then
expanding the cover; nonzero temporal relations remain a possible fallback.

See coded_state/LAYERED_CROSSOVER_COVER.md for the proof, exact coverage,
and serial reproduction. Retained only source, compact weights, and notes.
Regression checks pass: 14 exact integer/fractional equivalences, 684 local
rows, 160 direct bounds under both route envelopes, the matching moment
tests, 4044 count-box inequalities, and the existing location/all-pair tests.
The old single-profile exact replay remains -5153. Both frozen helper
hashes are unchanged. The new three-region replay has no failed cells.
C: has about 43.1 GiB free; unrelated dirty work is preserved.
No benchmark, construction change, paper edit, commit, push, or promotion.

### 2026-10-23 — Neighboring empty-pair profiles and larger count regions

Kept the construction fixed. Extended the location-aware count to boxes by
retaining the shared dependence of the route dimension and affine count
terms. The max-min inequality gives an integer upper bound over each entire
box. Convex active-row bounds cover every s=1..8192, not just test values.

Certified gapless regions at R=190 and R=191 with z=0..49168 and independent
+/-128 variation in (r,f,e,B) around (24568,589632,368520,49136). Their
good-route count-union log2 bounds are -824 and -3933. At R=192, certified
z=0..49200 with +/-64 variation (-707), and z=0..48000 with +/-128 (-417).
The four certificates use 4,6,9,9 boxes respectively. They include integer
count unions; no passing sample is substituted for region coverage.

Replaced loose integer-power rounding in the new verifier with directed
fixed-point rational log2 bounds. This removes thousands of bits of
accumulated numerical slack without changing the inequality or construction.
The old single-profile verifier and frozen certificate helpers are preserved.

The fixed-weight neighbor screen now passes 237/243 selected profiles. The
worst remaining point has upper +2306. Per-s selection among eight marker
sets improves it to +1760 but does not close it. R=193 and R=194 at the
central counts still fail, with +12704 and +37637. These are failed bounds,
not low-weight witnesses. No complete layered-code certificate is claimed.

New checks cover 4044 shared-dimension inequalities, 640 direct expression
comparisons, directed rational logarithms, complete interval replay, and
honest reporting of an uncertified cell after the subdivision budget ends.
See coded_state/LAYERED_COUNT_REGIONS.md for derivation and reproduction.

Final serial validation passes the new logarithm and region tests, the
location-aware algebra tests, and the exhaustive all-pair component tests.
The original single-profile exact verifier still returns -5153 with no gaps;
all four new regions replay with no failed cells. The two frozen certificate
helper hashes are unchanged. No raw outputs were retained; C: has about
46.6 GiB free. Unrelated dirty work is preserved.

Next: refine proof weights at the active-row crossover, cover the six failed
neighbors and R=193, then expand the count domain. If weight refinement
stalls, retain additional nonzero temporal relations jointly with outer
constraints. Keep the GF32 sandwich as the certified whole-code control.
No benchmark, production/paper change, commit, push, or promotion.

### 2026-10-23 — Location-aware count closes the empty-pair obstruction

User requested joint forced-zero locations, surviving components, and outer
constraints, testing the hardest remaining profile with the construction fixed.
Kept RS[256,128], the regional permutations, and the layered inner unchanged.

The profile R=192,z=49136,r=24568,f=589632,e=368520,B=49136 at W=838860 now
closes across every active-outer-row count. The exact good-route class upper
bound is 2^-5153, replacing a failed log2 upper bound of +491591. This is one
complete aggregate history class at 40%, not a whole-code certificate.

The new proof retains forced-zero AND forced-nonzero fixed inputs. It
stratifies by the active outer-row set, charges zeros only when routed into
active rows, and requires nonzero inputs to route there. Independent Bernoulli
colours conditioned on each region's exact row count reproduce the actual
permutation membership distribution. The normalization is explicit; a tilt
of those colours is only a proof parameter, not a setup change.

An ambiguous fixed input is split into zero/nonzero cases. Its zero case
adds one independent temporal equation to the actual state space, but not a
second reduction of the auxiliary kernel dimension. This tradeoff stays
inside the same component/exception transfer. Both one-update and two-update
operators are used, with exact rational inputs and row maxima.

Seven compact marker proposals suffice; six are used by the interval cover.
The earlier count handles active-row counts 1..5900 with log2 upper -22137.
Six contiguous intervals cover 5901..8192. Discrete convexity in the active
row count certifies interiors from exact endpoint bounds. The worst endpoint
upper is -5166 at precision 2048; a 13-bit union allowance gives -5153.
No sampling or optimizer convergence flag is used as certificate evidence.

The exact quotient identity kappa+K_R-rank(H G_fix) is checked in 128 GF4
fixtures. In 63, relations among nonzero fixed inputs add rank beyond the
forced-zero constraints. Other checks: 972 exact colour bounds, 528 actual
local classifications, 1188 exact/float rows, and 11 independent input-zero
equations that leave kappa unchanged. Full derivation and reproduction are
in coded_state/LAYERED_OUTER_LOCATIONS.md.

Final serial validation reran the affine/fixed-rank, inner-law, mixed-history,
old component, all-pair component, and new location tests successfully. The
gapless exact active-row cover replayed successfully. Both frozen certificate
helper hashes are unchanged. Runs used Python -B, retained no raw outputs,
and left approximately 51.8 GB free on C:. Existing unrelated work is intact;
this research workstream remains untracked.

Next: apply the location-aware count to nearby empty-pair profiles, then
certify count boxes toward full coverage. Keep nonzero temporal relations
as a possible next refinement, not an assumed generic-rank gain. Preserve
the GF32 sandwich as the whole-code control. No production/paper edit,
benchmark, commit, push, or promotion is part of this result.

### 2026-10-23 — All old failures close; all pair types and count boxes

User requested remaining failures, other pair types, and progress toward
full coverage, with the construction fixed. Reproduced the original 55,161
selected profiles, including 9,194 necessary exclusions and 40 failures among
45,967 retained profiles. The original marker set closes 24 failures. Added
all 01 transitions and one rational set closing the 13 r-positive failures.
Three r-zero failures at R=176 require three-step label lookahead; a third
rational set closes them. All 40 now pass exact replay, worst upper −47,363.

Extended the component operator to 00 as well: 100 internal transition terms
and 13 terminal branch terms on nine incoming states. Removed forced-active
restrictions whenever positive 00/01 counts are allowed. Exact tests enumerate
2,207 histories, check path-component equality for every one, and verify a
complete weighted partition bound. No construction or setup change.

Added a rigorous box majorant: the route dimension decreases with z,e and
increases with r,f; all other log terms are affine in the five counts.
Verified 40 radius-128 neighborhoods, keeping z=0 and keeping r=0 where needed.
Worst per-class bound −36,212; after its integer-count union, −36,187. All 424
distinct corners agree with direct exact replays; the analytic majorant
covers the interiors. These boxes are not full parameter-space coverage.

With all pair types enabled, the R=192,r=24568,f=589632,e=368520,B=49136 profile
closes every z in 0..8192 through 33 verified gapless segments. Its count-union
upper is −10,577, conditional on the existing good-route event. The z=49136
profile still fails a bounded one-update reweighting: float +586,317, exact
+594,094. This is not a low-weight witness.

Composing two actual successive updates retains correlations between pairs.
Exact 81-by-81 local composition validates this two-update contraction. Its
retained bound is +556,423: improvement, not closure.
A deterministic forced-input-zero count J0=max(0,2z-B-b) gives K_eff=508080
instead of K_R=524288, improving the exact bound to +491,591. Checked 24,768
shortening inequalities and the zero floor on the exhaustive histories.
Even the optimistic independent rank of all 2z empty-pair positions would
not close this bound at the retained interpolation weights. The next step
needs joint outer/history information, not just a larger local lookahead.

Added LAYERED_ALL_PAIR_TYPES.md, all-type source/tests, a compact rational
marker bank, and an integer box verifier. CWC keeps route conditioning,
class closure, interval coverage, and whole-code obligations separate.
No benchmark, production/paper edit, commit, push, or promotion.

Final serial regressions passed: affine/fixed-rank tests, mixed-history tests,
the earlier hard-class and 27-neighbor replay, the 2,048-defect restricted
cover, all 40 new failure replays, the count boxes and empty-pair interval,
and the exact two-update composition. The restricted defect cover retains
its earlier 52-bit result; it is not a whole-code certificate. SHA-256 hashes
of persistent_mds_certificate.py and stratified_certificate.py match the
frozen controls. Runs used Python -B and wrote no raw data or bytecode caches.
C: has approximately 52.5 GB free after validation. The workstream remains
untracked; unrelated working-tree files were left untouched.

Next: retain forced-zero input placements jointly with components and outer
constraints for the substantial 00 profile. Test the ideal possible gain
before choosing another relaxation. The 40 count boxes and z=0..8192 interval
are saved controls; keep the construction fixed and do not extrapolate them
to full coverage or near-half distance.

### 2026-10-23 — Hardest mixed class closes through the joint count

User requested surviving components counted jointly with exceptional updates,
construction fixed, hardest class first. Kept the layered encoder, route,
GF32 sandwich, and frozen helpers unchanged. No native benchmark or promotion.

Built a three-state temporal component transfer. A live component is charged
only when finalized; subsequent pins may otherwise kill it. The first attempt
retained the exact nine-bit schedule but relaxed agreement of the two pair
endpoints' branch tags. That bound remained too loose. Kept a bounded source
diagnostic, not a certificate or claimed worst-case component improvement.

Enforced shared pair branches and interpolated the two inverse amplitude
bounds on the SAME history. This gives weighted component factor2^5 and
zero-rank factor2^-27 at theta=27/32. Positive potentials depending on current
and next label bits retain two-step pairing lookahead. Pointwise local row
bounds compose over disjoint pairs; no temporal independence is assumed.
The free terminal state and first output are included explicitly.

Exact rational replay closes R192,f589632,r0,e368520,B49136,W838860.
The route dimension bound is recomputed as163923. Conditional on the existing
simultaneous good-route event, integer-power log bounds give-38885 at p256
and-44123 at p512, versus the earlier+734997 failed bound. This closes one
whole count class, including all its placements and amplitudes, not the full
code. The same rational markers pass27 selected neighbors (f/e offsets1024,
B offsets512); worst p256 bound-14603. No interval coverage is claimed.

Exact tests enumerate650 width-four histories, compare transfer components
with the path oracle, and check their complete weighted sum against the
positive-potential bound. Added LAYERED_COMPONENT_COUNT.md, bounded proposal
source, and an exact replay source; indexed them. CWC distinguishes the
fixed-history lemma, class sum, good-route conditioning, and global gap.
The old desired uniform kappa<=407155 cap is not proved or needed here.
Final serial regressions preserve the256 affine-consistency fixtures,
864 kernel patterns,9728 inverse-fiber checks, mixed-history exact checks,
and the full100892800-class2048-defect restricted cover (worst-1663).
Both frozen helper hashes match. More than52GiB remains free; no raw dumps,
builds, paper edits, production edits, commits, or pushes were made.

Next: apply the new count to the other retained failing profiles, then extend
the local tables to01/00 pairs and obtain verified coverage. Keep the40%
exploratory closure distinct from the near-half goal. Do not modify the
construction or claim a whole-code certificate from selected passing classes.

### 2026-10-22 — Joint route test and fixed-coordinate consistency

Kept the layered construction and GF32 sandwich control fixed. The new joint
head/internal-edge envelope conditions on actual head counts and bounds the
remaining hypergeometric matching count. It passes42 exact range inequalities
and783 exact conditional probability checks. On the retained hard profile,
however, the dimension bound remains163923 and the best log2 upper bound
remains+734997. No broad coverage run was warranted.

Instead retained the fixed input coordinates in the inverse argument.
For a fixed history, let a be its state-space dimension, D its free-input
count, e its exceptional count, d the outer shortening dimension, and kappa
the kernel dimension of the fixed-input map on that state space. The new
amplitude exponent is min(a+d-D-e, kappa+K_R-D-e). The second branch does not
require substituting an upper bound for d into a claimed rank saving.
This is a proved fixed-history counting lemma, not a whole-code result.

The kernel is exactly a count of unpinned components along temporal coordinate
paths. Implemented both component union and a two-state scan. Checked256 GF4
routed-RS fixtures and864 exhaustive path patterns against Gaussian rank,
including a negative control where multiple fixed equations are redundant.
Two small fixtures show a positive saving; that frequency is not a K20 claim.

For the hard profile, the retained pattern count and new exponent would close
its -100 target if kappa<=407155 uniformly. The conditional bound is-107;
the required kernel cap has NOT been proved. Added LAYERED_FIXED_CONSISTENCY.md,
layered_joint_route.py, and layered_fixed_rank.py and indexed them. CWC keeps
the exact lemma, toy checks, failed route bound, and open numerical target
separate. No construction, native benchmark, production/paper change, commit,
push, promotion, or raw data output. All tests run serially and bounded.
Final regressions reproduce9728 affine-fiber checks, the mixed-history exact
tests, and the full2048-defect restricted cover:100892800 classes,1019 inverse
intervals, worst-1663. Both frozen helper hashes match. Free disk exceeds53GiB.

Next: sum or bound unpinned temporal components jointly with exceptional
locations and support changes, starting with the retained hard profile.
The two-state scan gives a concrete transfer interface. Repeated pins along
one component cannot be counted as independent constraints. If that count
fails, identify the surviving long components before another bound refinement.

### 2026-10-22 — Mixed histories: transition accounting and a defect extension

Construction and GF32 sandwich remain unchanged. Used CWC to record the
completed restricted extension separately from failed broader diagnostics
in LAYERED_MIXED.md. No native timing, library/paper edit, or promotion.

Extended the one-sided40%/52-bit result to at most2048 arbitrary00/11 pair
locations after the first batch. Repairing defects to type10 preserves the
01 count; temporal rank loses at most2H. Forward cancellation roots lose at
most2H. Counted all defect locations, support types, and exceptional masks.
Complete integer coverage checks100892800 classes, using38 forward suffixes
and1019 inverse intervals. Largest interval bound:-1663. At4096 defects the
same bound fails its -80 target at R206,r416587 (upper bound-63).

For unrestricted mixtures of10/01/11, introduced the number B of newly active
coordinates. Few activations force next-state cancellations; many force output
weight and temporal zero rank. Counted support changes and exceptional masks
with integer binomial upper bounds. Full nine-step wiring supplies necessary
orientation inequalities; a small exact min-plus dual strengthens them.
This rejects aggregate profiles that maintain a large static support with
only the favorable orientation, which is impossible under changing bits.

Replaced the crude subset-union head bound in the diagnostic by an additional
ordered-head envelope. Excess head degrees above each threshold have an exact
binomial expectation and a without-replacement martingale bound. Union over
thresholds/suffixes costs <2^-63. Existing controls and certificate helpers
are untouched. The old head and matching events are retained too.

The final selected screen has55161 profiles;9194 violate necessary wiring
constraints, and40 of the45967 remaining profiles still fail. Worst retained:
R192,f589632,r0,e368520,B49136, dimension upper163923, rank lower585828,
best log2 upper+734997. This is a failed upper bound, not a bad-codeword
witness. Most remaining failures involve slowly moving holes and many full
exceptional updates; a broad coverage run is premature.

Exact tests:4096 defect supports,729 small-cycle orientations,400 local root
cases,87552 nine-bit intervals,4000 joint transition/rank/root fixtures,
5049 binomial bounds,18 exhaustive occupation-dual comparisons,1080 ordered
head identities,576 swap-range checks. Runs were serial and bounded. Only
source and compact notes retained, with free disk above53GiB at final check.
Final regressions preserve the original one-sided closure (the added V-W
rank term improves its largest interval bound to-9631), the9728 affine
checks, and all old matching-graph thresholds. Frozen helper hashes match.

Next: bound head bias and internal matching edges jointly for the SAME
outer-block subset, instead of letting separate worst cases combine. If
insufficient, use the fixed-coordinate outer consistency at exceptional
full pairs. Keep the construction fixed; no evidence yet requires a change.

### 2026-10-22 — Sustained one-sided histories: complete restricted coverage

Kept the construction, GF32 sandwich control, and frozen helpers unchanged.
Retained the orientation of one-sided pairs: type10 has one free inverse
coordinate, whereas type01 has two and cannot take the exceptional branch.
The earlier shorthand treating every one-sided pair alike was conservative.
Added LAYERED_ONE_SIDED.md and layered_one_sided_gate.py; CWC separates the
restricted theorem, exact checks, and remaining whole-code obligation.

The actual route yields uniform head-subsets as well as uniform matchings.
Proved a simultaneous head-count envelope and combined it with the existing
edge envelope to bound every head-plus-selected-tail free set. On the forward
side, maintaining one nonzero coordinate per pair under changing pair bits
forces b/4 fresh cancellation roots per internal transition. On the inverse
side, support intersections limit potential output zeros; each such zero
imposes an independent edge equation along a temporal coordinate path.

Complete integer coverage closes40% (W838860), with failure below2^-52,
for all J-valued outputs supported on a regional suffix whose incoming
pairs after its first batch are all one-sided. This includes arbitrarily
long half-populated histories, not just the previous short warm-up budget.
Covered100892800 (R,r) classes:39 suffixes by the forward bound, the rest
by919 verified gapless intervals. Largest accepted interval log2 bound is
-7487. Exact tests pass371 head factorial checks,6144 oriented graph cases,
768 wiring-mask cases, and290 GF4 temporal-rank cases. A repeated final run
reproduced these receipts. Runs at41% and45% fail the per-class -80 target;
those failed bounds are not low-distance witnesses.

Also derived a stronger dense-family defect measure: only2*n00+n10 weight
units are lost, rather than twice every nonfull pair. Type01 is free of that
penalty. The old W854892 threshold and52-bit restricted margin apply whenever
W+2*n00+n10<=854892. For arbitrary mixed histories, amplitude counting gives
V-(D+e)=e-n01; what remains is joint geometry/rank and pattern counting.

No whole-code certificate, scalar-extension shortcut, native timing, library
or paper edit, commit, push, or promotion. No raw dumps or large artifacts.
All numerical runs were serial and bounded; source and compact receipts only.
Final regressions reproduced the9728 affine-fiber checks and both earlier
restricted receipts, plus all100892800 old matching-graph classes. Both
frozen certificate-helper SHA256 hashes match the retained values. Free
local disk remains above55GiB.

Next: combine the near-dense and sustained one-sided arguments for mixed
00/10/01/11 histories, retaining actual temporal rank and free-set geometry.
First test whether sparse/full transitions force enough roots or output
weight, before spending effort on a large numerical coverage search.

### 2026-10-22 — Layered inner: obstruction checks and two restricted closures

Kept the pair-shear construction and GF32 sandwich control unchanged. Used
CWC for LAYERED_INNER_PROOF.md. No implementation, timing, paper, or library
changes; no commits or raw artifacts.

Found a genuine limitation of the singleton-diffusion shortcut: once activated,
u*e_0 can persist with one input and one output symbol per batch for every key
setup. The prescribed input amplitudes depend on keys and need not be outer
codewords. Also proved a local weight-two output reset has probability1/(q-1),
not q^-512. A full local-wiring cycle is not a uniform refresh. None of these
facts is a whole-code distance counterexample.

Derived the actual affine branch bound q^d(F)*(q-4)^(-D-e): regular pairs have
two free inverse-input coordinates, one-sided pairs have one, exceptional
p=0 pairs have none and cost one root. This retains the fixed coordinates
that the density-only shortcut discards. Added a forward exceptional-count
bound and a sequential common-quotient rank interface using k conditioned
on a (uniform on S/a), without double-charging a shared coefficient.

Closed two restricted classes analytically, with exact integer receipts:
- Over GF32, all outputs of weight<=floor(.45N) whose incoming histories have
  at most8192 zero symbols in total have failure probability<2^-78655, for
  every route. This is not a dense-state theorem for arbitrary GF128 histories;
  full scalar-extension transfer awaits coverage of all GF32 classes.
- Over GF128 or GF32, the route excludes every nonzero output confined to at
  most23 fixed coordinate indices with probability>1-2^-51, for any pair keys.
  The same coarse bound for24 coordinates fails the40-bit target.

test_layered_inner.py passes9728 exact GF16 affine-fiber checks, the exhaustive
255-input/81-setup GF4 width-four cycle census, sparse-path/reset checks, and
both integer receipts. Tests are serial, bounded, source-only. Free disk was
56.43GiB before runs. The80-digit GF32 sandwich46% control replay reproduces
all459 rectangles/4196352 pairs and partition faf9d54b...b00bef6e, with the
same worst log2 bound below-65.77481065587761.

Remaining obstruction: the worst-case shortening bound d(F)<=|F|/2 is too
loose for intermediate active histories. For a dense suffix of length~.75N
with~.15N exceptional pairs, both current coarse forward and inverse bounds
are positive. Actual regional occupancy or retained consistency rank could
close that gap; there is no bad-codeword witness. Extended zero prefixes,
changing sparse supports, and multiple episodes are not covered by the two
restricted certificates.

Continued within this turn rather than stopping at that gap. Recognized that
free sets for full incoming pairs select edges in independent uniform perfect
matchings of outer blocks. Proved a factorial-moment domination by a binomial
law, then a simultaneous internal-edge envelope for every vertex subset and
regional suffix, except with probability<=2^-56. This bounds shortened outer
dimension for ALL free-pair masks at once. The former R192/e314572 diagnostic
now closes at45%, with conditional log2 bound below-1018495.

Added LAYERED_ROUTE_GRAPH.md and layered_route_gate.py. The graph tests pass203
exact factorial-moment checks and6144 selected-edge cases. Tested a heap-based
integer sweep, then replaced the large sweep with an exact piecewise-affine
minimization, checked against exhaustive small sweeps and all four retained
large receipts. All128 regional suffixes and100892800 exception-count classes
are covered at W854892 (~40.7644%), for fully populated post-activation states.
The combined restricted-class failure is below2^-52, including the route.

Extended this analytically to at most6144 nonfull pair locations after the
first batch of the first active region. Condition their keys for each fixed
output; uniform graph bounds avoid an extra mask union. Each ignored pair
costs at most two weight units, giving W842604 (~40.1785%) with the same52-bit
restricted margin. This includes arbitrary starting batches within a region
and a bounded warm-up/holes budget. The statements remain over GF32; full
scalar extension awaits complete base-field coverage. No whole-code claim.

Next discriminating test: retain the useful free coordinate in one-sided
pairs instead of discarding both, and combine this with coordinate-excursion
rank. Address long sparse periods and multiple episodes; the current gap is
histories exceeding the6144 nonfull-pair budget, not a need to change the
construction or tune its implementation. Keep the GF32 sandwich fixed.

### 2026-10-22 — First-principles proof audit and a cheaper certified field

User requested analysis of actual proof requirements and simpler mechanisms,
not more micro-optimization. Used CWC; no native benchmark or production edit.

Separated state dimension, output/reset counting, and anti-concentration
against the routed outer. Almost-uniform feedback on every nonzero state
is a sufficient interface, not the distance requirement itself. Found that
the old running-sum proof discards most reset rank on short excursions.
For r closed episodes among v batches, support rank R obeys
w<=2R+b(v-2r); unrestricted reset rank is br. Added the resulting two output
count bounds as an opt-in refinement, preserving old defaults and hashes.

Whole-code GF32 certificates now close45% and46% for RS256/state512 at K20.
Directed80/120-digit checks and independent coverage controls pass289/459
rectangles, respectively, covering4196352 pairs each. Worst log2 class bounds
are below-63.65608870/-65.77481066. After the class union and actual sandwich
density inflation <2^(1/1024), both retain40 bits (46% retains >42 bits).
Proved exact scalar-extension distance equality: GF32 coefficients give the
same symbol distance over GF128, without a four-lane union penalty. No new
timing is assigned, and this does not assert uniform GF128 lane mixing.

Rejected permanently independent tiny state groups structurally: one RS word
supported on the last129 regions gives distance <=8385*g/131072 for group
size g at state512. This excludes g1/g2/g4 for45%, regardless of coefficients.

Specified a cheaper layered pair-shear candidate. Each batch pairs one index
bit, cycling over nine bits. Two shears form A; output uses A+delta*I and
state update is equivalently delta*Q+Y. Simple setup restrictions make both
pair maps superregular, at one variable product and one fixed GF256 scaling
per output, with no large FFT. A singleton spreads across all512 state
coordinates after nine empty-input batches; this does not cover arbitrary
states. Derived exact regular and exceptional feedback point laws. The open
obligation is counting dependent zero-state history constraints and fresh
key exceptions together with outer membership. No whole-code pair claim.

Added INNER_REQUIREMENTS_REVIEW.md and test_inner_requirements.py. Exact tests
pass68 reset-count cases, four extension examples,255 GF16 feedback laws over
180 key pairs, and160 diffusion fixtures. Existing wrapper and sandwich tests
pass; the old GF128 receipt reproduces its original partition and bound.
All tests/screens serial and bounded; no raw files, commits, or promotion.

Next: retain the GF32 sandwich as a certified cheaper control, and derive
a multi-batch rank/output count for the pair-shear family before tuning it.
Start with repeated sparse states and reset schedules; escalate mixing only
when a concrete obstruction or accounting loss identifies the need.

### 2026-10-22 — Native gate for certified running-sum feedback

Implemented the actual two-diagonal/Cauchy feedback and running-sum reset
in native/sandwich_native.cpp, with the retained optimized RS256 outer and
regional route. No production or paper changes. The ideal-setup45%/40-bit
certificate is unchanged. The other RS512/state256 profile is not yet timed
end to end; both256/512 fixed maps pass algebra checks.

Used a full GF128 quadratic tower over the existing GF32 arithmetic so
the GFNI outer and small-field FFT coefficients require no basis conversion.
An independent scalar oracle checks1024 products and two complete512-symbol
inner batches. Trace, Cauchy columns, involution, inverse, specialized/control
map equality, and routed/staged/materialized equivalence pass. Optimized
and UBSan checks pass on Peach. Full-field keys remain nonzero and independent
in the ideal model; the reproducible timing fixture is not a crypto sampler.

Serial native iteration: first exact full encoder43.0ms; specialize fixed
FFT constants to reduce its map from9.7 to6.9ms; test constant-reduction
specialization and outlining (no meaningful additional benefit); separate
random gathers into batch staging. Final31-sample K20 medians: outer3.206ms,
outer+route6.681ms, complete contiguous inner27.260ms, direct fused full38.209ms,
materialized full34.493ms, batch-staged full34.730ms. Staging saves the32MiB
routed buffer but is not measurably faster than materializing. All full
paths agree exactly. This is slower than the retained16.6ms RS reference.

The principal obstruction is now measured arithmetic cost, not a missing
45% proof: one flat tower diagonal takes10.946ms, and two plus D take27.244ms.
A separate sparse-byte CLMUL representation diagnostic takes6.648ms for one
diagonal; it cannot be substituted without adapting maps or paying conversion
costs. No extrapolated full-encoder win is claimed. SANDWICH_NATIVE.md records
scope, commands, hashes, controls, and compact timings under CWC.

Runs were serial, time/memory bounded, with >20GiB local and63GiB remote
free. Only source and compact notes were retained locally. No raw dumps,
commits, pushes, library changes, or promotion. Existing proof helpers and
the frozen reset/outer controls were preserved.

Next: evaluate full-field representation and fixed-MDS cost together,
including conversions. If still too expensive, reduce full-field products
through a new feedback distribution and prove that actual distribution.
Do not replace full-field keys with small keys under the existing lemma.

### 2026-10-22 — Close45%, remove the reset transforms, cheapen feedback

The user requested completion of45%, continued feedback reduction, and
distance treated as a dependent variable. Closed group32 for both profiles
using exact regional moments and outward-rounded80/120-digit replays.
Independent integer-slice coverage checks and removed/duplicated/moved-box
controls pass. Group32 receipts have300/309 rectangles and conservative
44/42-bit margins. Frozen proof-helper hashes remain unchanged.

Iteration revealed that the same proof closes d=1: Q'=Q+Y, retaining all
b=256/512 state symbols. Both directed precisions pass814/462 rectangles,
covering16781312/4196352 feasible pairs. Worst class bounds are below
-66.46197877725/-65.77481065587, giving conservative41/42-bit margins.
This eliminates reset FFTs and translations without reducing state dimension.

Proved feedback R=diag(u) D diag(v), with independent nonzero full-field
diagonals and fixed superregular D. For every fixed z!=0, every point has
probability <=(q-1)^-b. Counting zero output coordinates using Cauchy minors
proves this directly. Its total trajectory inflation is at most
(1-1/q)^-N<2; both forward and inverse counting transfer, so the running-sum
variant retains45%/40 bits. Exact GF4/GF8 tests cover4110 first-scaled
vectors,252 support laws,74 minors, and shortcut negative controls.

The new feedback costs2N full-field products and one fixed MDS map per
batch. It removes one map relative to the retained exact-uniform two-map
construction. Raw keys use64MiB rather than8/16GiB for literal dense
feedback. No native performance claim is made: the current2.8/2ms figures
are older reset-only group32/group8 timings, not new full-encoder times.

Distance diagnostics: group32 has complete FLOAT covers at47% for both
profiles; the coarse48/49% per-class tests fail. Running sum closes45%
but its coarse47% test fails. Group8/16 pass45% FLOAT screens but are
superseded for this checkpoint by the certified running sum. State64/128
running-sum45% screens fail; these are bounds, not structural exclusions.

Added GROUPED_45_CERTIFICATE.md, grouped_reset_certificate.py,
test_grouped_certificate.py, and test_sandwich_feedback.py. All numerical
runs were serial and bounded; >20GiB free. No native benchmarks, remote
jobs, generated dumps, production/paper edits, commits, or promotion.
CWC kept exact local facts, complete certificates, FLOAT diagnostics,
and untimed operation savings separate.

Next: implement and time the two-diagonal MDS feedback with the simple
running sum, then compare the full inner and fused outer/route. Preserve
the certified controls. Small-field random diagonals do NOT inherit the
lemma; a cheaper coefficient distribution needs explicit root/rank accounting.

Exploratory update (2026-10-22): the user allows45% as an intermediate
target. Without changing group32, its simpler two-dimensional bound now
has complete FLOAT covers for both outer/state profiles:300 rectangles
for RS512/state256 and309 for RS256/state512. Worst class bounds are
-69.79359 and-65.77481, against thresholds-65 and-63 that budget a40-bit
whole-code union. This is not yet a directed certificate. The next gate
is independent coverage verification and directed replay; no additional
mixing or larger reset groups are needed for this proposal. Reset-only
timings remain about2.8ms and still exclude uniform feedback, outer, and
route. Serial bounded runs used under8 seconds total, with over20GiB free;
no generated data or native benchmarks. Reproduction is recorded in
coded_state/GROUPED_RESET_JOINT.md. Preserve the49.1% work unchanged.

The user resumed the purpose-built inner investigation. Preserve the fixed
outer/route profiles and seek a cheap inner with >=40-bit setup-failure
margin at K20. Distance is now a dependent variable:45% is a completed
intermediate checkpoint, while49% remains desirable for exploration.
The [coded-state index](coded_state/DESIGN.md) and latest ledger entries
record the current candidates. The [ideal reset ladder](coded_state/IDEAL_RESET_LADDER.md)
now separates outer concentration from state size and reset algebra. It
has uniform-feedback whole-code proof controls, not a selected fast kernel.
The64/128-state ensemble is excluded at49.1%; the256-state running sum
remains open. Preserve the certified window and all earlier partial proofs.

## Earlier goal: near-50% distance, preserving the finite window design

Updated 2026-10-21 by explicit user direction. Save and park the small-state
35% investigation, then drive the certified full-field window design closer
to 50%. Completing or benchmarking the 35% proposal is not the next task.

Keep K=2^20, rate 1/2, GF(2^128) symbols, fixed RS[256,128], one stratified
regional permutation, and independently sampled uniform Toeplitz feedback.
Start from memory L=256 and batch b=64, certified at >45%/72 bits. The working
near-half target is 49% with at least 40 bits of setup-failure margin.
First exploit proof slack without changing the encoder; then consider L/b
changes within the same design. Preserve all existing certificates. No
production, paper, commit, or push changes are implied.

The [new finite checkpoint](NEAR_HALF_WINDOW.md) now closes >49.1% with at
least 40 bits at (L,b)=(1024,256), retaining the outer and route exactly.
Both directed precisions and independent coverage checks pass. The
near-half proof milestone is achieved; performance for this profile is
unmeasured. Next measure and optimize the same certified Toeplitz design.

The [parked checkpoint](PARKED_SMALL_STATE.md) records the small-state evidence
and its current bound's approximately 37.86% ceiling. That bound is not a
candidate for the new high-distance task. The prior stages below are history.

## Historical hill-climbing goal (superseded)

The user requested hill climbing from the reviewed uniform-model proof.
Replace one expensive component at a time, preserving a whole-code bound
for the actual distribution. The [finite checkpoint](STRATIFIED_CLOSED_CERTIFICATE.md)
now closes fixed RS plus one stratified route and uniform Toeplitz feedback:
45% distance/72-bit margin at (B,L,b)=(256,256,64), and
35%/40 bits at (64,128,64), for K=2^20 and rate 1/2 over GF(2^128).

The goal remains active. The native exact-uniform window control now takes
43--44 ms, far above the target. A new [persistent MDS state](PERSISTENT_MDS.md)
has complete finite checks at 35% over GF(2^32), using state64 with B64 or
state32 with B128. Its additive-coset implementation now takes about23--24 ms
for B64/state64, after removing separate moment and state-multiplication work.
Both independent native oracles and UBSan pass. Packet-route certificates
now include four-record B64/state128 at35%/59 bits. Moving packets halves
standalone route cost, but total encoding remains about23--25 ms because
the certified profiles require a larger outer or state. A full Cantor
schedule now reduces the best scalar profile to21.5--21.8 ms and the
B64 packet profile to22.0--22.7 ms. Their certificates are unchanged.
The compact-key A/B is complete: no reliable encoding win. The experimental
one-diagonal profile now takes about 16--17 ms, but has no whole-code
distance certificate. Partial-packet moments and delayed syndrome rank
are proved and checked; completing their support count is the next gate.
Preserve the window control and its timings.
The target remains competitive with the preserved 3.798 ms GF256 baseline,
not merely another proof of 35% distance with more encoder work.

The current structural priority is the historical 64-input/16-state geometry
with full-GF32 nonzero DHD diagonals. Its [sparse gate](SMALL_STATE_PROOF_GATE.md)
covers every one-/two-group message at 35%, below 2^-405 combined.
The [direct paired-route transfer](SMALL_STATE_ROUTE_GATE.md) extends
coverage to one through six groups, below 2^-44 combined. A new
[whole-state operator](SMALL_STATE_RANK_BRIDGE.md) addresses every actual
source size. Directed gates now cover all messages on at most 128 groups,
below 2^-44 combined, and the separate intervals 2048..2073 and 4096..4108,
each below 2^-100. Two-batch fractional grouping and direct actual-size
conditioning close the former interval without adding encoder work.
Counts 129..2047, 2074..4095, and 4109..8192
remain open. Continue with complete intermediate anchors and directed
interval extensions; selected-profile scores are not coverage.
This is a partial proof with no new timing. The complete one-diagonal proof
remains open and preserved; it is not the performance destination.

The earlier coupled-scan objective and checkpoints below remain preserved.
Uniform full-field operations are allowed as reference components in this
new ladder; removing their cost is part of the open goal.

## Historical target and acceptance criteria (superseded)

Explore the coupled scan with one cache-local structured global route and
small coefficients over GF(2^128). Initially retain the old grouped outer and
route to isolate the inner change. These choices may change if a structural
obstruction justifies it; a new choice must receive its own specification.

- First distance checkpoint: K=2^20, rate 1/2, delta=0.35, and ideal setup
  failure at most 2^-40 over all nonzero messages chosen after setup.
- Later distance objective: move toward 50%, without demanding the exact
  rate-half limit.
- Performance reference: approximately 3.8 ms K20 precomputed forward encoding
  of 128-bit records. Compare matched implementations; do not infer timing
  from operation counts or mix historical machines and harnesses.
- No general large-field multiplication in the proposed hot loop. Preserve
  fused route reads and one final output write unless measured evidence
  supports a change.
- Keep original SPIN, the certified GF256 reference, and previous cheap-lane
  proofs intact. No production or paper promotion is implied.

A complete certificate must cover every nonzero message class, use the actual
sampling distributions and terminal convention, and have independently checked
proof accounting. Experiments and partial covers are not substitutes.

## Historical work stages (superseded)

1. Specify the construction and audit its local algebra.
2. Test coefficient-independent reset and support obstructions at widths 2, 4,
   and 8, including unions of outer groups.
3. Build the smallest sound rank continuation bound for a width that passes
   those checks. Width two has been rejected; widths four and eight are next.
   Include unions of active outer groups from the outset.
4. Extend to all support classes only if the sparse gate is promising.
5. Implement a matched native kernel after structural checks pass; benchmark
   serially and adjust proof and construction together.

When progress stalls, identify whether the failure is a real bad codeword,
a proven ceiling, a lossy bound, or a numerical limitation. Do not expand
coefficient searches without a specific explanation of what they can resolve.

## Progress ledger

### 2026-10-20 Initial specification

- Created a separate workstream and active goal; prior work is preserved by
  links rather than copied proof trees.
- Specified exact field, outer, regional pair distribution, inner coefficients,
  state persistence, inverse, and free terminal state.
- Independently reviewed invertibility for w<16, whole-vector persistence,
  width-two empty-input support bounds, and the fixed-mask rank lemma.
- Derived the complete-batch reset obstruction and the full-support condition
  for extinction by a single allowed input coordinate.
- Passed bounded GF(16) checks: 14 determinants, 33,824 local state cases,
  262,144 one-step inverses, 65,536 width-four zero-output checks, and 128
  full-stream inverse fixtures, including an outside-family negative control.
  Root ran them serially with Python -B; the tool call took about 2.7 seconds.
- Whole-code margin: not computed. Native timing: not measured.

### 2026-10-20 Width-two rejection

The first group-union reset argument rules out width two for the initial
outer and route. Restricting to 32 adjacent outer groups leaves dimension
4096. Killing seven post-input gaps per region costs at most 3584 equations.
Uniform pair placement makes the expected deleted fraction exceed 0.71 in
each region. Independence across regions and Hoeffding imply setup failure
at 35% distance with probability above 0.841, regardless of inner coefficients
or outer mixer choices. See [the proof](WIDTH_TWO_OBSTRUCTION.md).

Exact rational checks confirm the two deletion expectations and the Hoeffding
constant. This is an obstruction to the code family, not a failed lower bound.
Do not benchmark or pursue a 35% certificate for this width-two instantiation.
An independent proof review found no quantifier, rank, or probability gap.
The reproducible `theory/check_width_two_obstruction.py` check passed with
Python -B, reporting a failure lower bound above 0.841689975667; the note
uses the conservative value 0.841.

Paired adjacent inputs also expose a cheap singleton-state reset at every
width with natural-cycle wiring. Pair-antipodal cycle orders may remove that
local pattern without an extra global route; their full effect is unresolved.

Next test: extend the group-union reset bound to widths four and eight, and
compare natural versus pair-antipodal wiring. Preserve the local rank lemma,
but design the continuation bound around these multi-group obstructions.
No bulk proof files, raw timing data, or builds were created in this round.

### 2026-10-20 Wider states and a constructive proof path

Classification of the preceding turn: progress. It specified the family and
proved a structural rejection of width two. This continuation adds the
following evidence rather than repeating a stalled parameter search.

- Extended the cluster-gap argument using actual pair-slot placements and
  inward rounding to complete batches. Width four fails 45% distance with
  probability above 0.018, independently of all coefficients. This does not
  reject width four at 35%. The width-eight expression gives no useful
  sub-50% rejection. Both the finite proof and asymptotic interpretation
  received an independent read-only review.
- Exhaustively checked necessary reset masks at widths 2, 4, and 8. A single
  antipodal pair requires state support at least w/2, versus one for natural
  wiring. At width eight, four of the six two-pair unions again permit
  singleton reset masks. These masks are not whole-code witnesses.
- Derived a positive local rank transfer, exact paired-region coefficient
  normalization, and composition in one common quotient with free terminal
  state. Implemented exact rational local kernels with bounded caching.
- Established a support-containment bridge from the actual GF256 outer
  distribution to larger-field shortened dimensions. Identity outer mixers
  admit distance at most 33N/256, so a proof uniform over all outer mixer
  realizations is impossible. This exceptional realization does not refute
  the averaged randomized construction.
- Found and repaired a proof relaxation: repeated singleton localization
  appeared to have zero fresh charge. Fixing the route before inner coins
  permits one-step lookahead. A run of L such steps has L-1 fresh charges,
  without changing the encoder or adding a randomizer. The global numerical
  effect of this improvement has not been evaluated.
  An independent review confirmed the lookahead exposure and composition;
  a separate code review found no error in the quotient tests or charge kernels.

Root ran the following bounded, serial Python -B checks; no files were
generated by them:

| Check | Evidence |
|---|---|
| Reset bounds | 10,816 exact rounding cases, 58 tiny-route expectation profiles, rational K20 failure lower bound >0.018098410105. |
| Support wiring | All masks at widths 2, 4, 8; 28 single-coordinate and 44 physical-pair-union contexts. |
| Rank transfer | 1,250 GF16 quotient fixtures, 34 exact pair profiles covering 6,294 labeled permutation/swap choices, 12 generating rows, 16 two-group shortening fixtures. |
| Lookahead | 818 two-step quotient fixtures, 6,886 retained mask branches including 238 with extra charge; 780 localization runs of lengths one through four. |

All six indexed check commands pass in a fresh serial replay. Actual
linear-form rank fixtures use width two, dimension two, and tau=2 over
GF16; the lookahead fixture list is not exhaustive over all forms. Wider
support and generating-function checks have the narrower scopes stated
above. The analytic lemmas, rather than those fixtures, cover general widths.

These are local checks, not a complete distance certificate. The next step
is a bounded width-four/35% proof screen, with width eight retained for
higher-distance headroom. Start by testing source-support counting and the
lookahead continuation on sparse group unions and a dense class. If a bound
fails, identify the specific lost constraint before expanding the search.
No native implementation, benchmark, production change, or paper edit was
made. The goal remains active.

### 2026-10-20 Dense band, exact support bridge, and first sparse closure

Classification: progress. The code is unchanged. This round produced two
partial distance bounds and identified a concrete accounting loss in the
sparse continuation; it did not expand a blind parameter sweep.

- Closed all source supports with at least 99% nonzero coordinates at 35%
  output distance. The ideal-setup failure bound is below 2^-167780,
  conditional on any injective outer and route. Exact arithmetic checks
  include 640 local rows and the denominator-cleared integer inequality.
- Proved an exact-support bridge from GF256 words to the F-linear shortened
  spaces. Its uniform per-group factor is less than 3.31. Root checked
  277 small GF4 codes, 1091 exact nonzero supports, and their GF16 extensions.
  Keep the containment bridge: its uniform point-cap envelope is tighter
  than the new uniform-factor envelope. This is not an automatic exponent win.
- Evaluated the ordinary lookahead sparse gate with five markers at tau=256.
  All five bounds exceeded one; the smallest diagnostic log2 bound was
  about 2264.59. No construction failure follows from that relaxed screen.
- Identified the lost finite state-span budget. An empty interval retains
  the same span of at most w state forms. A live branch can gain at most
  w-1 independent ranks across the entire interval. Exceptional-root
  cancellations remain possible and retain their probability costs.
- Implemented and independently reviewed a Bellman bound that distinguishes
  those cases. Root checked 280 exact GF16 fixtures over one and two empty
  steps. A float comparison at z=.994, tau256, and 2048 empty steps improves
  the worst rank-four log2 continuation from about -25.13 to -58.49.
  Those values are diagnostics, not a whole-code failure margin.
- Obtained a common positive envelope and verified it exactly. The initial
  4096-iteration proposal cap was insufficient for one vector; the bounded
  checker now allows at most 16384 iterations. At z=.998 and lambda=.993,
  all 32 live-state inequalities pass after rationalizing the proposal.
- Added an ordinary singleton-step interface with outgoing state rank at
  most incoming rank plus one. Fixed-route linear envelopes precede route
  averaging, avoiding an invalid interchange with the nonlinear Bellman
  maximum. The resulting 2-by-2 regional bound retains the free terminal.
- Closed the exact full-support one-group class at 35%, summing all 8192
  group choices: failure below 2^-849. The local constants use exact
  fractions; the final positive calculation rounds upward at 80 and 120
  decimal digits. Both give an upper bound below 1.897e-256.
- The same envelope's 33-full-column bound is about 1.77316e31 and therefore
  uninformative. This is still an open proof gate, not a bad-code witness.

Next: refine the 33-column gate before expanding the support search. It has
only four global message dimensions; retain that global budget across
occupied steps or sharpen their joint coefficient outcomes. Do not reset
global rank at gap boundaries, and do not average a nonlinear Bellman
operator over route positions. Then extend to partial columns and two
adjacent groups using the actual outer distribution.

All numerical runs were bounded and serial. The workstream contains only
source and compact notes; no raw data, caches, builds, native benchmarks,
production changes, or paper edits were created. The goal remains active.

All eleven indexed check commands pass in a fresh serial replay, including
the deliberately still-open 33-column screen reported by the final command.
The workstream occupies approximately 195 KB across 26 source/note files;
its only subdirectory is `theory/`. The existing tracked research-index
diff passes the whitespace check. No commit or push was made.

An independent final read-through found no defect in the state-rank
boundary, common linear route envelope, source multiplicities, or directed
rounding. It specifically checked the strict upward comparison with
2^-849 and the 33-column recurrence's single counting of column choices.

### 2026-10-20 Separate rank costs and close the minimum-column class

Classification: progress. Two further sparse classes close without changing
the encoder. The key change separates an independent equation's rank weight
from an exceptional coefficient's actual probability.

- Proved that the empty-span Bellman operator permits any tau>=1. It keeps
  the factors 1/tau and 1/256 distinct. Ordinary active-step transfers must
  instead use the common charge parameter min(tau,256). The same global
  potential persists across every step; no rank or probability is reset.
- At tau=65536, z=499/500, lambda=993/1000, all 32 exact envelope inequalities
  pass after 16 bounded proposal iterations. The union over all 33..64
  fully populated columns in one group is below 2^-291, including every
  column set and all 8192 group choices. The h=33 term is largest.
  The final positive calculation rounds upward at 80 and 120 decimal digits.
- Derived a shared-outer pair-count bridge at tau=q^2. At exactly 33
  occupied RS columns, a shortened row has one amplitude; a pair's column
  rank is consequently constant. Toeplitz zero constraints have rank at
  least a+t-1 for a zero output coordinates and t-dimensional input span.
  This controls partial scalar masks using the actual mixer distribution.
  It does not assume independent mixers for the two members of a pair.
- The first partial-mask bound at z=.998, lambda=.993 was uninformative,
  about 6.45e230. Five bounded alternative envelopes showed that the output
  marker was too aggressive for this broader class. At z=.999, lambda=.9962,
  tau65536, all 32 inequalities pass exactly after 48 proposal iterations.
  Directed 80/120-digit replay gives failure below 1.084e-66 and explicitly
  verifies 2^-219. This covers every mask on exactly 33 occupied columns,
  all their positions, and every active group, with a free terminal state.
- Extended the empty-gap test to 840 exact GF16 fixtures at tau=1,2,8 with
  q=2. Added 1250 fixed-mask fixtures at tau8>q2, confirming the capped
  common charge. Added 750 exact GF2 Toeplitz subspace fixtures and 15
  ordered coefficient-recurrence fixtures. Root ran all checks serially.
- An independent read-only audit checked the shared-setup pair count,
  extension-field dimensions, Toeplitz conditioning, source multiplicities,
  matrix order, terminal convention, and directed arithmetic. It found no
  defect. The numerical values are root's replay, not independent agent runs.
- Preserved an unimplemented 108-state global-rank refinement separately.
  It is not needed for either new closure. The deterministic information-set
  cap in the partial-column note may help the remaining classes.

The new sparse classes overlap. They do not imply a full certificate.
In particular, pair ranks need not remain constant across 34 or more occupied
columns. Reusing the 33-column product bound there would be unjustified.
Multiple active groups also require the actual paired-route correlations;
they cannot be replaced by independent singleton insertions.

Next: study h=34 partial columns using the shortened degree-one row structure
and information-set dimension cap, while auditing two adjacent active groups
as the first genuinely paired input class. Select a continuation refinement
only after identifying the lost constraint. Keep the dense band and all
proved sparse classes as fixed regression checks. No native timing,
production change, paper edit, raw result archive, commit, or push was made.
The goal remains active.


All fourteen indexed check commands pass in a fresh serial replay. The
older tau=256 command deliberately retains its failed 33-column diagnostic;
the newer commands verify the closures above. The workstream occupies
approximately 245 KB across 32 source/note files, with only `theory/` as
a subdirectory. No caches or generated result files were created. The
research-index diff passes the whitespace check.

### 2026-10-20 Every one-group support and the first paired classes

Classification: progress. The construction is unchanged. This round closes
the entire one-active-group class and a first genuinely paired input class.
It does not establish a whole-code certificate.

- Proved an exact-support pair bridge at tau=q^2: for every realized
  shortened space, its initial dimension factor is at most its exact-joint-
  support pair count divided by q^2-|S|. The pair only bounds the initial
  dimension; it need not be an inner witness.
- Replaced constant column rank by polynomial pair profiles. A nonzero
  minor limits the number of rank-one evaluations. A stronger line-constraint
  lemma retains every such constraint: a dependency in their Vandermonde
  system forces the entire row-pair kernel into one rational line. This
  repairs the information-set count's loss beyond the first b positions.
- Combined the actual Toeplitz column bounds with a positive ordered
  (occupied-column count, rank-one count) recurrence. Four exactly verified
  marker pairs bound the same h-classes. Taking the minimum of complete
  bounds per h and summing gives failure below 2^-147 at 35%, including
  every one-group scalar support and all 8192 group choices. Directed
  80/120-digit replay gives less than 5.511e-45. The worst h is 44.
- Root's exact toy tests cover 2343 GF4 line systems, 910 deficient kernels,
  36 pair-bridge fixtures, and 20 ordered profile coefficients. An independent
  read-only audit found no defect in the proof counts or implementation.
- Derived the actual adjacent-group route: alternate regions contain one
  full pair or two distinct singleton pairs. The latter collide within a
  batch with probability 1/4095. Both pair parities and the free terminal
  remain in the positive region products.
- Retained independent source pivots in the first b=h-32 occupied complete
  columns. This lowers the h=64 bound from about 6.109e843 to 2.121e766,
  but does not close it. The hypothesis is not available for arbitrary
  partial columns or lookahead constraints.
- Identified a second loss: the scalar reset envelope forgets the sparse
  state needed by a cheap paired reset. Applying two exact empty Bellman
  steps to the physical-pair reset terminal before collapsing improves
  its normalized long-gap constant by a factor greater than 21000.
- Integrated four-step tails into the actual route via positive ordered
  matrix-series convolution. All zero/short gaps and collisions remain;
  no nonlinear Bellman operator is averaged over routes. This closes the
  union over two adjacent groups on the same 33..35 fully populated columns
  below 2^-117 at 35%. Both 80/120-digit comparisons pass; the bound is
  below 5.510e-36. The all-h adjacent-group bound remains uninformative,
  about 6.993e689. This is a failed bound, not a low-distance witness.
- Independent proof/code audits found no discrepancy in the endpoint
  refinement or its O(m ell) convolution. Exact small-route tests cover
  six fixtures, 72 full-pair and 784 ordered singleton-pair samples,
  24 insertion coefficients, 72 endpoint convolutions, 24 endpoint regional
  averages, and eight zero-tail recoveries of the old formula.

Next: diagnose the large-h adjacent-group loss before enlarging the class.
The current envelope still spends the full initial dimension factor while
discarding source/state relations after the information-set prefix. Separate
that loss from the remaining reset-boundary loss: inspect the dominant
quiet/live products and test a common global-rank continuation. Choose class-
appropriate rank/output markers only to test a stated hypothesis, not as an
unbounded search. After this gate, extend to different column sets, partial
columns, nonadjacent pairs, and larger group unions.

No native implementation or benchmark, production change, paper edit,
raw-result archive, commit, or push was made. Prior proof checkpoints remain
indexed. The research goal remains active.

All twenty indexed check commands pass in a fresh serial replay. The old
sparse and broad two-group commands deliberately report their failed bounds;
their successful execution is not a certificate for those open classes.
The workstream contains 41 source/note files, approximately 340 KB, with
only `theory/` as a subdirectory. No generated files or caches were created.
The tracked research-index diff passes the whitespace check.

### 2026-10-20 Large paired classes: rank, unused outer randomness, and online routing

Classification: progress in the proof machinery, but no new closed source
class. The construction and all previous partial certificates are unchanged.

- Separated the initial dimension factor, output-threshold factor, and
  regional live/reset costs at four fixed marker choices. Smaller rank
  price helps the large-h class; merely extending the reset tail does not.
- Implemented the paired augmented-span operator. It tracks the old-state
  rank, rank after adjoining current inputs, and which inputs lie outside
  the old span. An independent algebra/code audit found no discrepancy.
  Root's exact tests passed 144 continuation comparisons, including 1297
  retained branches, 1137 old-span pivots, 1244 outside pivots, and 34 resets.
- Directed 80/120-digit checks of that operator with the actual endpoint
  route give a broad shared-full-column union about 9.928e207. It improves
  the preceding bound only modestly and does not close the class.
- Derived a lazy raw-Toeplitz exposure for the final 32 occupied columns
  per group. The prefix remains sampled invertibly for source freshness.
  First-nonzero pre-mixer partitions use the same global quotient; later
  input matches pay unused diagonal coefficients. Conditioning the 64 raw
  matrices back to invertible costs (256/255)^64 once. This changes no
  encoder operation. Directed 80/120-digit evaluation gives about 5.751e162,
  still above one. The removed-reset diagnostic is only about 4.93 bits
  smaller, identifying the live continuation as the main current loss.
- Proved a separate online route recurrence for ordinary/fresh chronological
  equations: enumerate support histories, expose each batch's pair labels
  and swaps, and apply a suffix bound uniform in all remaining labeled forms.
  This avoids the invalid interchange of nonlinear fixed-route bounds and
  route expectations. The lazy-outer extension has an enlarged filtration;
  its obligations are recorded explicitly. An independent proof/code review
  found them satisfied for the shared-full-column class with delayed input
  nonzero checks, genuinely invertible prefix mixers, and no lookahead.
- Implemented the regional recurrence without replacing pairs by independent
  samples or discarding collisions. The first combined float screen gives
  log2 union 286.351, with h33=-43.545, h48=117.969, and h64=286.350.
  That is about 254 bits smaller than the endpoint/lazy calculation but
  remains a failed broad bound. These floating values are diagnostic only.

Next: retain residual ranks 1..4 separately across region boundaries, instead
of collapsing all live states to one value. Test that change first on the
large-h adjacent-pair class, with the same markers and source distribution.
If it helps, audit and verify its positive regional composition before
expanding the source classes. No blind marker sweep or native implementation
is justified by the present failed bound.

No encoder benchmark, production or paper edit, raw-result archive, commit,
or push was made. Compact source and notes preserve all older proof paths.
The goal remains active; no whole-code distance guarantee has been obtained.

Root replayed the twenty earlier indexed commands and the new paired-rank
and online-route tests serially; all 22 passed their stated checks. The
online tests cover 12 exact route-law/linear recurrences over 215 labeled
placement/swap samples, 30 local terminal/boundary comparisons, and 396
directed 80/120-digit entries against exact three-batch nonlinear tables.
These tests do not assert a general interchange of nonlinear averaging.
Older failed-bound diagnostics
remain failed and are not interpreted as certificates. The workstream has
only compact source and notes (about 0.4 MB), with no generated caches or
result archives. The tracked research-index diff passes the whitespace check.

### 2026-10-20 Uniform reference before further hill climbing

Classification: an all-message asymptotic proof, not another partial support
class and not a finite certificate for the optimized construction.

- Specified independent random local [B,RB] outer codes, one uniform global
  permutation, and independent uniform recursive feedback of memory L.
  Setup is fixed for all messages; the terminal state remains free.
- Generalized the exact input/output enumerator. Identified why its naive
  unconditional word count can fail at large alphabets: rare bad routes
  contain entire subspaces of low-weight words.
- Replaced that count by conditioning on two properties of the routed outer:
  small unions of c intervals support dimension at most k(c-1), and larger
  unions have shortened dimension at most their random-code benchmark plus
  o(N). Independent block sampling and concentration establish both properties.
- Counted low-weight outputs through the inverse inner. Its random input
  coordinates form an affine fiber, with one fresh activation per interval.
  Forced zero tails close small unions; the q-ary entropy exponent closes
  larger unions. B,L=Theta(log N) suffice for every fixed strict sub-GV
  distance. These constants have not been made useful at K20.
- Root and two analytic reviewers found no gap in the whole-code argument.
  A bounded exact root check verified activity identities and assignment
  counts for 262140 length-eight nonzero words over four symbols, memories
  one through four, and 322 interval representations. No files were generated.

For context, the rate-half GV boundary is approximately 38.0472% over GF256
and 49.2189% over GF(2^128). These are asymptotic reference distances, not
finite achieved operating points or new performance measurements.

The previous interrupted rank-boundary screen is also preserved. Retaining
five boundary states gave floating log2 upper bound 203.3545 for the broad
adjacent-pair class, still above one. Its optional active-rank-floor
extension has not been evaluated; neither is a new certificate. These
diagnostics do not drive the current uniform-baseline task.

Next: extract finite outer/inner proof interfaces, then replace one reference
component at a time. A block-dispersing route should address terminal
concentration before cheapening the inner randomness. Preserve the 3.798 ms
GF256 baseline and all older partial results. No native benchmark, production
or paper edit, raw-result archive, commit, or push was made.

### 2026-10-20--21 Fixed RS, stratified routing, and Toeplitz batches

Classification: complete finite whole-code certificates for a new reference
ladder, not closure of the old coupled scan and not a new timing result.

- Created the requested unbudgeted hill-climbing goal. Retained the uniform
  theorem and every earlier proof checkpoint.
- Replaced random dense injections by fixed subfield RS blocks. One
  stratified permutation places one coordinate per block in each region.
  Exact MDS shortening, a fixed-size regional moment comparison, and convex
  binomial envelopes bound every interval union.
- Applied a fractional moment to conditional failure, not to the original
  unconditional bad-word count. This avoids overcharging rare bad routes.
  Complementing coordinate sets before the moment comparison preserves
  their fixed total occupancy and removes another substantial loss.
- Closed RS[256,128], scalar memory 256, 45% symbol distance, with at least
  66 bits of setup-failure margin at K=2^20 and rate 1/2 over GF(2^128).
  The complete cover has 1240 rectangles; exact rational route moments
  and directed 80/120-digit checks both pass.
- Proved joint uniformity of fresh Toeplitz feedback on every nonzero
  state. This permits batching without replacing the actual inverse
  distribution by independent marginal approximations. Keeping D(A),
  the shortening on random inverse coordinates, is essential.
- Closed memory 256/batch 64 at the same 45% distance, with at least
  72 bits of margin. Then closed RS[128,64] and RS[64,32] at 35% using
  memory 128/batch 32, each with at least 70 bits.
- Smaller-memory probes exposed slack in the output-volume bound. Added
  the q-ary Chernoff entropy bound, preserving the actual alphabet.
  It does not close every failed probe, but allows RS[64,32], memory 128,
  batch 64 at 35% with at least 40 bits. This final cover has 19935 boxes;
  both directed precisions verify every pair's contribution <=2^-69.
  There are fewer than 2^29 pairs. No source/support class is omitted.
- The 35%/40-bit candidate stores 191/64 full-field diagonals per output
  symbol, instead of 128 independent coefficients for scalar feedback.
  Direct Toeplitz multiplication still costs 128 full-field products per
  output. Actual arithmetic reduction requires a fast middle-product
  kernel; the storage reduction is not a timing claim.
- Independent proof and implementation reviews passed. Root's exact
  regression checks include regional fixed-size distributions, convex
  envelopes, complement identities, scalar/batch interval counts,
  Toeplitz joint uniformity, affine zero-dimension cases, output volumes,
  and small rectangular bounds. All numerical work was bounded and serial.

See [the finite certificate](STRATIFIED_CLOSED_CERTIFICATE.md) for the exact
distributions, table, commands, complete-coverage argument, and failed
probes. An unclosed envelope is not a bad-code witness. In particular,
the first rejected pair is not reported as a global maximum.

Next: implement a precomputed Toeplitz middle product for the certified
RS[64,32], memory-128, batch-64 design. Compare it with batch 32 and the
preserved 3.798 ms GF256 code under matched conditions. Include routing,
full-field products, transformed setup, and scratch traffic. If arithmetic
or route traffic still dominates, change only that component and prove
the new distribution. Do not infer practical speed from the proof alone.

A reviewed exact implementation plan uses evaluation at the 256 GF256
subfield points. Reduction modulo Z^256-Z leaves the needed Toeplitz
middle coefficients unchanged for L=128,b=64. With precomputed coefficient
evaluations, it requires four full-field pointwise products per output,
plus two transforms and any polynomial/field-basis conversions. The
certificate note records the indexing and aliasing argument. This is the
next kernel to test, not a new encoding-time result.

The goal remains active. Production, paper, and older proof files are
unchanged. Only compact research notes and source were added; no generated
data, cache directories, build trees, commits, or pushes were made.

### 2026-10-21 Native gate and conversion-free uniform maps

Classification: progress on exact implementation and a new proved inner;
the practical performance target remains unmet.

- Implemented the RS[64,32], L128/b64 full-field profile. The first exact
  Toeplitz kernel passed direct polynomial and full-encoder oracles, but
  took 201.36 ms at K20. Isolated transforms accounted for about 163 ms;
  256 full-field products per batch accounted for about 26.4 ms.
- Preserved that implementation as a selectable baseline. A derived and
  constexpr-validated Cantor basis makes monomial conversion XOR-only,
  reducing the full time to about 157 ms. This alone is insufficient.
- Proved a different batch map H diag(r) G. The state is directly a
  degree<128 novel-basis polynomial; H extracts the highest64 novel
  coefficients of 256-point interpolation. Any64 H columns are independent.
  Nonzero states have at least129 nonzero evaluations, so the full-field
  diagonal yields exactly joint-uniform64 feedback coordinates.
- Restricted independent random entries to192 points. At least65 active
  evaluations remain for every nonzero state. The same uniformity proof
  applies, so the complete >35%/40-bit certificate transfers unchanged.
  No claim about a particular deterministic benchmark seed is made.
- Implemented the exact three-coset kernel: six64-point transforms,
  192 products, and local XOR/scaling passes. Its output matches the full
  256-point control and an independent polynomial oracle. The outer now
  accepts novel-basis coordinates, which preserves the fixed RS code.
- Serial Peach/GCC15/CPU2 medians give43.407 and43.667 ms for the final
  profile, versus48.295 ms for its unpruned control. GFNI-only field
  multiplication is slower at46.748 ms. No simultaneous benchmarks ran.
- The final profile has96 MiB of keys. Standalone contiguous inner costs
  about29.4 ms and routed inner38.7 ms. Normalized point-product work alone
  is about19.8 ms. These component runs are diagnostic, not additive.
- Direct checks cover monomial and novel transforms, top-coefficient
  extraction, dense feedback, two arithmetic kernels, regional bijections,
  and the complete free-terminal recurrence. The final profile also passes
  the undefined-behavior sanitizer. Analytic reviewers agree on the map,
  pruning, and proof transfer.
- Proved a restricted cost bound: an exact-uniform F-bilinear batch with
  L inputs and b outputs requires at least L+b-1 variable F-products.
  Thus192 is near-minimal for this interface, not for all machine algorithms.

See SPECTRAL_BATCH.md and NATIVE_BASELINE.md for precise claims, commands,
source hashes, and measurements. Current best35% full-field native result
is43--44 ms, not the preserved3.798 ms GF256 result. No near-baseline speed
claim is justified. The45% certificate remains a proof checkpoint without
a native timing in this turn.

Next discriminating step: compare a cheaper realization of full-field
products with a lower-cost inner whose nonuniform inverse law can be
bounded explicitly. Keep the exact-uniform baseline as the control.
Do not return to an unlabelled small-coefficient substitution. Scalar
routing also costs roughly6 ms alone, so future packetization/fusion
must preserve or replace its actual route proof.

The goal remains active. The old candidates, paper, and production library
are unchanged. Local additions are source and compact notes only. The
bounded remote scratch directory totals about2 MiB; no raw data was retained.
No commit, push, promotion, or cleanup of unrelated user work occurred.

### 2026-10-21 Persistent state: smaller exact-uniform construction

Progress: a new whole-code proof, not a micro-optimization of the FFT.
Replace the sliding output window by S'=P S+C Y, with P diagonal over
GF(2^32) and C a fixed GF256 Vandermonde matrix. Scaled evaluation points
are distinct across the whole stream. Every closed state excursion then
has weight at least b+1 and imposes b independent output constraints.
Batch and state width are equal. The feedback remains exactly uniform
on F^b for every nonzero state, using a 2b-point spectral map.

Two independent envelopes are needed. The inverse envelope counts
low-weight outputs satisfying the reset constraints; the forward envelope
charges q^-b for each prescribed reset and counts projective outer inputs.
Taking the better complete envelope per class closes the finite sum.
The actual uniform regional route and fixed RS outer are unchanged.

Complete exact-moment directed replays, at both80 and120 digits:

| Field | B | Batch=state | Distance | Margin at least | Rectangles |
|---|---:|---:|---:|---:|---:|
| GF(2^32) |64|64|>35%|46 bits|2846|
| GF(2^32) |128|32|>35%|44 bits|2939|
| GF(2^128) |256|128|>45%|48 bits|710|

All rows use K=2^20 and rate1/2. Thresholds are734003 and943718.
This is a whole-code ideal-setup claim, not a guarantee about a particular
benchmark seed. Root ran all numerical jobs serially with bounded
coverage/time; no raw output or search state was saved.

Independent analytic review found no defect in the interval counting,
projective correction, complemented route moments, volume bounds, or
directed arithmetic. Exact local tests are recorded in the checkpoint.
The new helper leaves the frozen window-certificate sources untouched.

Unsuccessful probes are retained as failed upper bounds, not bad-code
counterexamples: B64/b16 and B64/b32 at35%, and B256/b64 at45%.
No parameter sweep or native timing was performed in this round.

Remaining obstruction: implementation cost, especially the exact field
products, the fixed state-moment map, and scalar route traffic. The GF32
profile can use four genuine field lanes per128-bit record. That is not
an implementation of a GF128 multiplier, nor a trivial GF256 embedding.
The GF128 row requires correct embedded-GF32 state multiplication in its
chosen field representation. Smaller state than batch remains uncertified:
nonzero vectors in ker C can bypass random feedback in a quiet batch.

Next discriminating test: implement B64/b64 and B128/b32 GF32 controls
with direct moment-map and feedback oracles, then measure outer, route,
inner, and full encoder serially against the preserved native control.
Do not infer runtime from field-operation counts or promote these profiles.
The goal remains active; paper and production code are unchanged.

### 2026-10-21 Native persistent state, additive coordinates, and packet routes

Classification: progress on proofs and implementation; performance target
still unmet. The full GF32 multiplicative control passed direct native
oracles and UBSan, but took34--35 ms for B64/state64 and30--31 ms for
B128/state32. These encode four GF32 payload lanes per128-bit record.

The additive-coset construction changes the state coordinates to reuse
feedback transforms in the update. It removes a separate moment map,
monomial conversion, and fixed GF32 state multiplication. Exact-uniform
feedback and global MDS reset constraints are preserved. The scalar-route
certificates transfer. Native direct Cauchy, full-stream, global-syndrome,
and sanitizer checks pass. B64/state64 now takes23--24 ms; its contiguous
inner falls from23.0 to10.4 ms. B128/state32 takes about25.4 ms.

To address routing, group p independent outer rows into fixed packets.
Each region independently permutes those groups, with fixed lane order.
For whole-packet sets, the shortened dimension is p times the sum of
group excess occupancies. The route moment consequently uses tau^p and
N/(pB) factors; output dimensions and reset probabilities stay scalar.
Complete directed80/120-digit checks give35% distance at K20, rate1/2:

| Packet records | Outer B | State=batch | Margin at least |
|---|---:|---:|---:|
|2|128|64|46 bits|
|4|128|128|48 bits|
|4|64|128|59 bits|

The third row first closed at40 bits using a finer rational marker grid.
Summing each certified rectangle with its own area improves the bound
to59 bits. That sum is an exact dyadic upper bound, independently reviewed;
it is not a floating log-sum approximation. All4096 toy shortened dimensions,
4096 complement identities, and255 exact toy moments pass.

Native packet order, regional bijections, direct stream/global-state
oracles, and UBSan pass. Serial31-sample medians on Peach CPU2 show route
cost4.06 ms for two-record packets and3.15 ms for four-record packets,
versus6.47 ms for scalar routing. Full timings do not improve: packet2
B128/state64 is23.45--23.49 ms; packet4 B64/state128 is23.18--24.05 ms;
packet4 B128/state128 is25.13--25.83 ms. A same-session scalar control
takes23.33 ms. This is not a packet end-to-end win.

Remaining obstruction: the exact-uniform inner still dominates, and the
route improvement costs larger certified state or outer parameters.
Fine-grid and lower per-term-target probes do not close packet2 B64/state64;
the failed bounds are not low-distance counterexamples.

Next discriminating tests: extend the GF256 Cantor chain to GF32, so upper
transform levels use cheap byte coefficients; then examine a one-diagonal
feedback law with a separate proof. The latter is not certified and its
inverse random support can cut packets, so the present route moment must
not be reused without justification. Old proofs, production, and paper
remain unchanged. All jobs were serial; only compact source/notes are kept.

The optional full32 Cantor schedule is implemented and independently
checked. It extends the original GF256 chain by bounded GF2 elimination,
verifies all chain relations and rank32, and makes upper transform levels
use byte coefficients. Distinct global cosets, exact-uniform feedback,
and every complete certificate survive this deterministic schedule change.
Both scalar profiles and both packet-four profiles pass direct native
oracles and UBSan. The packet sanitizer build required `-g0` to stay
inside the2 GiB compiler cap; sanitizer coverage was not removed.

Final same-source raw/full/full/raw comparisons show a modest improvement:
scalar B128/state32 gives21.441 and21.687 ms, against raw22.214 and22.264;
packet4 B64/state128 gives22.426 and22.074 ms, against raw23.072 and23.781.
Those profiles retain35%/44 bits and35%/59 bits, respectively. Earlier
control variability is retained in NATIVE_PERSISTENT.md; do not claim the
larger gain suggested by selecting only the slowest raw run. The other
scalar B64/state64 profile did not improve. No near3.8 ms claim is justified.

The next exact-distribution implementation test is compact raw GF32
feedback keys:16 MiB rather than64 MiB, trading more arithmetic for less
traffic. Existing tower arithmetic supplies a direct control and should
produce identical encoded outputs under the same seeds. Separately,
SINGLE_DIAGONAL.md records a lower-work recurrence and exact local bounds;
its randomized inverse support cuts packets, so the old certificate does
not transfer. That is the next substantial proof gate, not a claimed design.

The older field-kernel checkpoint already measures the raw-key tower
kernel slower in isolation (about1.44 versus1.07 ms per2^20 records).
Limit this diagnostic to one whole-encoder A/B for changed cache pressure;
do not assume compact keys win or restart an unbounded kernel search.

The goal remains active. No paper/production edits, commits, or pushes.
Local files are source and compact notes. Remote scratch totals42 MiB,
with about63 GiB free on its volume; no raw experiment dumps are kept.

### 2026-10-21 One-diagonal implementation and partial-support proof interfaces

Classification: progress on implementation and proof interfaces, not a
new whole-code certificate. The certified profiles remain unchanged.

- Completed the bounded compact-key comparison. Raw GF32 coefficients
  reduce key storage from 64 to 16 MiB and preserve the route, random draws,
  and encoded outputs. They slow scalar B128/state32 from about 21--22 ms
  to 25 ms. Packet B64/state128 results vary in both directions. Keep the
  optional flag, but stop this branch rather than retuning it indefinitely.
- Implemented default-off one-diagonal feedback: V=DZ, Y=X+r.*V,
  Z'=V+Y. It removes one D application and one random diagonal per batch;
  prepared keys fall to 32 MiB. The actual unshuffled packet route stays
  unchanged. The new mode is explicitly marked experimental/uncertified.
- All four exposed scalar/packet-four profiles pass direct scalar field,
  Cauchy-map, routing, full-stream, and global-syndrome oracles. All four
  also pass UBSan with recovery disabled. Native checks are implementation
  evidence, not distance certification.
- Serial 31-sample medians give 18.16--18.53 ms for scalar B128/state32,
  17.15 ms for scalar B64/state64, and 16.01--16.13 ms for packet-four
  B64/state128. The last same-source control/single/single/control replay
  gives 22.250/17.034/16.934/22.536 ms. Report the candidate as about 16--17 ms,
  not only its fastest run. Setup/allocation are excluded. The historical
  3.798 ms baseline was not rerun and remains far faster.
- Proved a partial-packet route moment by adding local lane shuffles in an
  auxiliary distribution. Monotone sampling without replacement, Holder,
  and Bernoulli equalization give the actual p-power envelope. A complete
  invariant support-mask sum has the same moment under the unshuffled
  route. This proof step adds no shuffle to the measured encoder.
- Grouping output masks W at fixed random-coordinate mask H permits the
  W count inside a fractional power. The H multiplicity remains outside.
  Applying that saving to all masks at once would be invalid. Per-batch
  multiplicities distinguish a previously troublesome periodic class;
  that diagnostic is partial and omits its finite boundary count.
- Derived exact linearized output dimensions by propagating the reachable
  syndrome space and its fibers. Retaining a known output-column core
  strengthens the dimension relaxation. Summing the constraints gives an
  amortized excursion-rank bound, with the free terminal charged explicitly.
- Proved that one observation batch sees a column-scaled Cauchy matrix
  from any collection of earlier outputs. Selecting observation times
  therefore charges delayed ranks on disjoint fresh-variable intervals.
  No arbitrary multi-time minor or generic-matrix claim is made.
- Exact serial tests pass 4096 routes across 126 support classes, 516096
  complement identities, 252 rational moments, 4589 monotone comparisons,
  56 generator-rank fixtures, and six mask-orbit checks. Separate rank tests
  pass 1026 dense fixtures, 2428 Cauchy minors, and all 65536 width-two,
  two-batch outputs. Analytic derivations, not these toy instances, support
  the general lemmas.
- An independent final read-only audit found no gap in the known-core cap,
  amortized excursion bound, selected-time Cauchy charges, projective
  counting, or terminal treatment. The partial-packet comparison and mask
  symmetrization also received independent analytic review.

Remaining obstruction: the complete amplitude/support sum. The amortized
rank bound combined only with D(A)<=|A|/2 leaves a positive dense exponent
above roughly one-third distance. Globally replacing per-batch support
counts by binomial(N,H) also discards too much. These are identified proof
losses, not low-distance witnesses. The two-diagonal reset probability
cannot be reused for the one-diagonal recurrence.

Next discriminating test: a bounded support-profile or positive-operator
calculation retaining per-batch counts, the known output-column core,
delayed Cauchy charges, actual partial-packet moments, and D(U)>0 geometry.
An ell-only state drops the fresh-column rank in persistent partial-mask
patterns; do not expand that relaxation blindly. Likewise, counting all
quiet/live paths without outer geometry includes impossible late-only
inputs and cannot close the full bound. First evaluate a small, explicitly
chosen family of difficult profiles and compare the complete valid
envelopes before constructing a large cover.

The goal remains active. No production/paper edits, promotion, commit,
push, raw data, or generated reports. Local additions are compact source
and notes; remote scratch totals 50 MiB with about 63 GiB free. All numerical
runs and benchmarks were serial, with explicit time and memory caps.

### 2026-10-21 Affine compatibility, overlap ranks, and bounded profile gate

Classification: proof progress; no new whole-code certificate and no
reliable full-encoder performance improvement in this round.

- Derived the exact affine membership condition for the inverse
  X=Y+diag(r)LY. Besides paying q^(-rank P_A), a proposed output must
  satisfy PY in image(P_A). Its deterministic compatibility rank is
  rank(P_(A union W))-rank(P_A). This uses the actual routed MDS outer,
  including partial packets; it requires no additional randomization.
- Proved a valid additive joint rank: compatibility equations vanish
  on A intersect W. Inner zero equations restricted to those overlap
  columns therefore contribute independent rank. The selected-time
  Cauchy bound applies with overlap counts as source counts. This is
  not the invalid sum of unrestricted outer and inner ranks.
- Grouped masks by S=A union W before fractionalization, retaining the
  number of S masks outside the fractional power. Both H and S groupings
  remain necessary. The complete mask family has the required packet
  symmetry; the implementation still performs no lane shuffle.
- Root's exact serial replay passes 10240 joint-rank fixtures, including
  overlap ranks, 13225 coefficient assignments, 40 exact joint spaces,
  and 6280 output assignments. Independent read-only audits agree with
  the algebra, grouping, and free-terminal treatment.
- Proved explicit nonzero rank-one 2-by-2 causal impulse minors for both
  b64 and b128 in the actual full-Cantor K20 prefix. The b128 example
  uses the AES reciprocal element 0x8d and a quadratic-extension root.
  Its four coordinates occur at batch/lane pairs (555,82), (556,55),
  (861,81), and (862,52). Exact native-tower tests pass. This rejects
  generic multitime-minor arguments; it does not exhibit a low-distance
  codeword or invalidate the single-observation Cauchy lemma.
- Added a bounded floating diagnostic for 18 named absolute support
  profiles. It retains activation, resets, free terminal, per-batch
  multiplicities, delayed ranks, and outer activity geometry. Its
  selected-time algorithm passes 64 exact comparisons with a quadratic
  recurrence. All eligible tested profiles have strongly negative
  expressions under at least one bound. Late-only activity is excluded
  by D_out(U)=0. This is neither a directed certificate nor coverage of
  all profiles, translations, or placements.
- A targeted nineteenth profile then exposes a real loss in these bounds:
  activation, 10734 batches of (h,w,o)=(65,1,1), then 5649 full batches.
  Its output weight is 733934, below the 734003 threshold. Even the full
  dimension-only relaxation leaves E-H/2 positive; the averaged H/S
  moment expressions also remain uninformative. This is a failed bound,
  not a codeword. The late dense suffix contains 22 complete regions,
  information discarded by those bounds. A convex occupancy argument
  gives D(H)<=533246 for every route, instead of H/2=710391. This makes
  the selected-time exponent E-H+D(H)<=-164398. Counting this profile's
  masks crudely by 2^(135*10734) proves its existence probability at most
  2^-3811646 over inner coefficients, uniformly over the route. This is
  an exact bound for one absolute profile, not a full-code margin.
  The diagnostic reproduces the improvement and passes 48105 pointwise
  and 144315 aggregate integer checks of the regional cap.
- Tested outer store tiling without changing the sampled code. A
  sixteen-row streaming tile improves the isolated outer from 4.140 to
  2.583 ms, but the gain largely disappears in full encoding. A final
  single-diagonal tile/control/tile replay is 16.797/16.489/16.732 ms.
  Keep the default unchanged; no reliable full-encoder win is claimed.
  Direct checks pass for both packet-four profiles. The tiled candidate
  also passed UBSan before the final SFENCE-to-MFENCE read-back barrier
  refinement; final native checks pass after it. The candidate remains
  uncertified and about 16--17 ms, far above the historical 3.798 ms.

Remaining obstruction: a complete support-profile sum. None of the
selected diagnostic patterns now exposes a failed bound, but that does
not make their finite list exhaustive. Arbitrary mixed profiles and
excursion placements remain unaccounted for. A generic-matrix theorem
cannot fill this gap because the explicit causal minors disprove it.

Next discriminating test: carry regional occupancy information into a
complete support count alongside delayed and overlap ranks. The phased
example rejects a universal negative E-H/2 claim from the current inner
rank inequalities alone. Its cure is outer geometry, not more inner
mixing. Treat almost-full regions robustly, for example by virtual
completion or retaining nonidentical regional marginals. Do not substitute
an unbounded pattern sweep for missing coverage. The affine and
causal-minor notes index the new lemmas and exact tests.

The goal remains active. No promotion, production/paper edits, commits,
pushes, or raw dumps. Local free space is about 62.7 GiB; remote scratch
is 56 MiB with about 63 GiB free. All numerical work was serial and
bounded. Final runner SHA-256 is
`d9be438a3bea45b26816028aae38c1ba87ececd54129a0df38d8b516d8491f49`.

### 2026-10-21: complete activity-band closure, not a whole-code certificate

Classification: progress. The encoder and setup distribution are unchanged.
No new benchmark was run. The certified two-diagonal checkpoints, the
uncertified 16--17 ms one-diagonal candidate, and the historical 3.798 ms
reference retain their prior status.

- Proved the stronger deterministic regional cap: outer shortening is at
  most the sum of the smallest half of the regional cardinalities. Nested
  incidence columns attain it in the relaxed incidence model. This handles
  almost-full regions and dominates the earlier full/touched secant; optimal
  virtual completion gives the same cap. Exact tests pass 135424 matrices,
  7467 count classes, 14196 completions, 37800 monotonicity checks, and
  32 large-integer fixtures. Packet-route attainability is not asserted.
- Proved a genuine bad-route obstruction, not a failed rank bound. A valid
  route places the same 65 groups in each region's first 65 packet slots.
  Their 8320 message dimensions exceed the 8192 regional reset equations.
  A nonzero word then has weight at most 24576, or 1.171875%, for every
  inner coefficient setup. That specified alignment has probability
  binom(8192,65)^(-64), not a claimed total failure probability. A whole-code
  proof must retain route randomness or establish a good-route event.
- Added the nested guarded moment E[tau^D(T) 1[D(U)>0], for T subset U
  with whole-packet U. Its bound min(H^G,G J H^(G-1)) retains the chance
  of a nonzero outer word. The auxiliary lane shuffle transfers only
  through the complete invariant mask sum. Exact tests pass 4096 labeled
  toy routes across 8 fixtures, 16 guarded bounds, 32 positive-envelope
  checks, and 16 orbit identities. A scalar-cut U counterexample rejects
  an invalid extension of the shared packet indicator.
- Derived a complete two-index excursion envelope by replacing D(T) with
  the common D(U), then summing all masks and output weights before the
  fractional power. The layout count is
  binom(v-c,c-1)*binom(n-v+c,c). The sharper touched-region cap is
  min(B,2c+floor((v-2c)/L)); it excludes classes before route equalization.
  Exact small tests pass 646 layout/interval checks. A bounded floating
  screen evaluates 101 complete (v,c) classes and excludes another 38.
  Most remaining scores are positive: this homogeneous envelope is not
  sufficient for a whole-code certificate. It is not a codeword search.
- Closed one entire activity band with exact integer arithmetic:
  8096<=v<=8192 and 1<=c<=32, covering all 3104 classes, every placement,
  every support mask, and all nonzero output weights <=734003. The common
  amplitude exponent is -2700496, hence C^nu <= (33/32)^(-84390).
  The actual whole-packet moment is bounded at regional mean 1/2.
  The final rational comparison proves failure below 2^-255 for this
  restricted union. It took 0.295 seconds remotely; this is checker runtime,
  not encoder performance. Independent source and algebra review pass.
  Reproduce with `python -B theory/check_active_union_band.py` from the
  workstream, under a 60-second/256-MiB cap. Other activity ranges remain open.
- Recorded analytic operator-cost reductions: a sparse-activity scalar
  bound leaves two separable products; positive chords can compress the
  8192-term expansion with a bounded loss; tangent bounds can remove
  regional counters. None is claimed to yield a contracting complete
  operator. The unverified small-state transfer proposal was not asserted.

Remaining obstruction: the complete sum outside the closed band. Replacing
the regional occupancy vector by its mean loses too much for fragmented
activity. Replacing D(A) by D(U) also loses the dense-pattern gains.
Neither more named-profile tests nor a deterministic all-route theorem
addresses these gaps.

Next discriminating test: sum regional excursion placements with the guard
retained for a bounded sparse-activity class. Use the two-product envelope
before considering a large general support operator. Keep the regional
information-set cap and delayed/overlap ranks as complementary dense-class
bounds. Reject the homogeneous active-union envelope as a sole whole-code
strategy; do not spend an unbounded sweep trying to repair it.

All numerical work was serial and source-only. The local volume had about
62.7 GiB free and remote /tmp about 62.3 GiB. Only a few small source files
were transferred into the existing scratch directory. No production or
paper changes, promotion, commits, pushes, raw dumps, or new build trees.
The goal remains active.

### 2026-10-21: outer-aware short excursions and structural cost review

Classification: proof progress, with no new full-code certificate or timing.
The performance objective remains approximately the historical 3.798 ms
reference, not merely a certificate for the current 16--17 ms candidate.

- Proved a limitation of the coarse active-union gate. For every route,
  some layout in the complete (v,c)=(66,33) class has positive outer
  shortening. Its coarse conditional amplitude bound is already one.
  Consequently exact regional route probabilities cannot make that
  complete coarse sum smaller than one. This is a loss in the bound,
  not an actual low-distance codeword. Stop the regional-only repair.
- Derived a forward repair that retains the number R of outer rows
  capable of contributing to the activity union. Sparse activation gives
  a state of weight at least b+1-R before a two-batch reset. Its actual
  reset probability is at most q^(-(b+1-R)). A group factorial moment
  bounds the exceptional large-R routes without assuming independent
  actual groups. Projective source counting occurs after fixing the route
  but before exposing inner coins.
- An exact integer checker closes all classes with 1<=v<=512 and
  3c-v-2>=ceil(v/4), subject to 2c-1<=v. There are 10,954 checked classes
  and 54 deterministic exclusions. Each inner and route contribution,
  including every placement, is below 2^-80. The complete restricted
  union is below 21908/2^80<2^-65. This includes (66,33), every mask, and
  every source amplitude. The serial run took 0.373 seconds under a
  60-second/256-MiB cap. Source and algebra received independent review.
  Reproduce with `python -B theory/check_short_excursion_forward.py`.
- Added joint forward support bounds and exact reset/activation cases.
  Nonnegative additive Fourier transforms bound every affine weighted
  sum by the homogeneous MDS split enumerator. The bound explicitly
  retains its impossible-origin overcount: subtracting that term to
  obtain a live-only bound is unjustified. No contraction or complete
  support count is claimed. Exact GF16 regression passes 3072 key
  assignments, 60 shell inequalities, 108 potential inequalities,
  3 split enumerators, 162 affine inequalities, 256 activations, and
  512 two-batch cases covering 131072 keys. The serial run used a
  30-second/256-MiB cap and wrote no files.
- Reviewed the cost structurally. The packet candidate's contiguous
  inner alone takes 9.311 ms. Compact keys and outer tiling did not give
  reliable full-encoder wins. The proposed cheaper interface returns to
  64 inputs and 16 state coordinates with full-GF32 nonzero DHD diagonals,
  retaining the historical GF256-valued outer mixers, paired route,
  and fixed A/C/H maps. It uses 0.5 general field products per output
  and one H16 per 64 outputs, rather than one D128 per 128 outputs.
  These are operation counts, not a timing prediction or a certificate.
  Independent review confirms transfer of H16 total nonsingularity and
  the quiet 40% event bound below 2^-76, including the free terminal.
  The nonquiet proof does not transfer automatically from the old field.
- Closed a first nonquiet gate for that small-state proposal. Full-GF32
  DHD has point density at most (q-1)^-16 for each fixed nonzero input.
  A delayed-exposure argument makes incoming liveness predictable before
  the fresh refresh, avoiding conditioning on future path events. Split
  source words by their exact occupied RS-column union and their number
  of state returns. Many returns cost full-state coincidences; few returns
  force a long live stretch. The complete one-group bound is below
  2^-411 and the two-group bound below 2^-406, hence their union is below
  2^-405. All group choices, source amplitudes, column sets, activity
  masks, and the free terminal are included. The bound is uniform over
  the fixed outer and paired route. The root's integer-only replay passes
  all 64 (group-count,column-count) classes and 8144 interval fixtures,
  serially under 60 seconds/256 MiB. Independent analytic review and root
  source review agree. No native implementation or timing was added.

Remaining obstruction: long excursions and the rest of the support sum
for the one-diagonal checkpoint; three-or-more-group counting for the
cheaper small-state proposal. At three groups with all 64 columns, the
coarse deterministic gate permits only 6118 live batches, too few to force
output zeros at 35%. This is a failed bound, not a low-distance witness.
For five or more groups, nonzero quiet inputs may also lie in ker C;
the all-quiet receipt does not cover partially quiet trajectories.
The GF256 outer density must not be relabeled as uniform-GF32 density.
Neither candidate yet meets the full goal.

Next discriminating test: retain actual paired-route placement in the
three-group small-state gate, or charge source constraints consumed by
returns. Do not spend another round improving its already excessive
singleton margin. Retain the forward kernel as a complementary proof
tool, not an excuse for an unbounded 16-ms proof-only detour. The old
certificates and all partial checkpoints remain indexed and unchanged.

No benchmarks ran in this round. All numerical jobs were serial. Remote
/tmp retained about 62.3 GiB free; only small source files were transferred.
No promotion, production/paper changes, commits, pushes, or raw data.
The full goal remains active.

### 2026-10-21: sparse transfer and the whole-state counting bridge

Classification: new partial proofs and a reviewed complete counting
interface, but no new whole-code certificate or encoding time.

- Transferred the retained direct paired-route bounds to full-GF32 DHD
  diagonals. Independent review checked the activation and projected-rank
  monotonicity, source-rank cutoff, and actual route. Together with the
  one-group gate, every message on one through six actual groups has
  setup-failure probability below 2^-44 at 35% distance. This is an
  analytic domination of existing directed receipts, not a new replay.
- Derived an all-source whole-state bridge. Fix output-zero requests and
  the zero/live path, count the E-linear trajectory kernel minus its
  zero-source fiber, and take the extension-degree root before expanding
  the positive output marker. The argument keeps actual-source support
  separate from auxiliary E support. Empty actual inputs cannot change
  whole-state liveness, and the terminal refresh is never charged.
- Derived two valid source-support corrections: a missing-group cutoff,
  and joint Jensen weighting that cancels exact zero-group stratum costs.
  A positive auxiliary-zero-state penalty is available in either version.
  Independent reviews covered the kernel count, conditional zero-state
  probability, pointwise stratum cancellation, and matrix interface.
- The exact GF4/GF16 checker passes 288 fixed event inequalities,
  288 zero-source fibers, 35 retained-support checks, one empty-input
  path, 135 root-before-marker inequalities, 180 projective-count
  inequalities, 7172 exact source strata, and 35 joint-weighting checks.
  The root reran it serially under a 30-second/256-MiB cap.
- Preserved a partly quiet reduction: deleting quiet constraints costs
  their rank and can be counted by live-run geometry. This is a valid
  unevaluated bound, not a certificate. Auxiliary E source classes cannot
  be deleted merely because actual F sources of those sizes are closed.
- Independent review caught an incorrect local layout in the initial
  floating diagnostic. The encoder stripes each pair across two RS16
  blocks; the diagnostic had placed it within one block. All initial
  scores are withdrawn. The analytic construction and transferred sparse
  proofs were unaffected. Corrected v2-striped code has an exact
  singleton-law regression and a negative control for the wrong layout.
- Ran the corrected cutoff and joint probes serially on the same 50
  selected pair profiles and 600 profile/marker choices each. Both pass
  15 small pair-polynomial checks, 36 matrix checks, 3 profile-mass checks,
  and 289 striped singleton cells. Each probe took about 5.2 seconds and
  50 MiB RSS, writing no files. They are floating diagnostics only.
  Selected g=4096 and g=8192 profiles have negative upper-bound scores;
  selected g=7,16,64,256,1024 still have positive scores. Joint weighting
  yields no qualitative improvement; omega=1 is best on this tested grid.
  See SMALL_STATE_BRIDGE_SCREEN.md for the compact corrected table.

Remaining obstruction: the complete sparse/intermediate source count for
the small-state proposal. The source-density substitution discards the
quiet constraints' RS ranks. The at-mean histogram divisor is too small
to explain the entire observed deficit, although exact coefficient
extraction can remove other losses. Further blind marker sweeps are not
the next step.

Next discriminating test: retain per-column quiet source rank and the
best-32-column nonzero-source cutoff at g=7. Use actual input ranks;
do not extend the g<=4 injectivity shortcut. An exact g=7 route coefficient
is a useful control. Preserve the new whole-state bridge for dense
classes and the separate direct gates for g<=6.

No native kernels, benchmarks, paper edits, promotion, commits, or pushes
were made. Only small source files were transferred to the existing
remote scratch directory. Both local and remote volumes retained more
than 20 GiB free. All jobs completed; no background calculation remains.
The full goal remains active.

### 2026-10-21: exact-route closure through sixteen source groups

Classification: substantial partial-proof progress. The cheaper small-state
construction is unchanged; no whole-code certificate or new runtime.

- Implemented and reviewed a direct GF32 rank-tagged seven-group gate,
  including actual diagonal rank probabilities, one common source quotient,
  the strict best-32-column cutoff, exact striped routing, and a free
  terminal. Its full floating bound improved from +1473.40 to +815.65 in
  log2 units after bounded rational-marker tuning. Raising tau to q-1
  without retuning was worse. Do not continue this seven-group search:
  the separate two-state route gate below closes it more simply.
- Exact small-field regressions pass 1875 source-rank fixtures, 339 full
  cutoffs, 625 threshold identities, 36 diagonal-rank caps over 3510
  assignments, 12 projected caps over 27000 assignments, and 18 MDS
  minors. These check local algebra, not the full finite code.
- Replaced the IID-conditioned regional route majorant by exact striped
  singleton/full-pair coefficient extraction in the two-state bridge.
  The seven-group log2 score changed from the older failed +10118.64
  selected-profile screen to -979.078576 at a fixed rational marker.
  The improvement removes off-histogram paths; it is much larger than
  merely removing the histogram-normalization penalty.
- Directed Arb checks at 128 and 256 bits close every source count
  g=7,...,14 below 2^-900 per count, and g=15,16 below 2^-100 per count.
  They cover all 440 candidate matching profiles at each precision,
  excluding the impossible simultaneous-full endpoint for even g.
  The g16 gate uses alpha=4/5 before route averaging; its binomial source
  count remains outside the fractional power. Its worst profile is near
  -135.28110 in log2 units. No cycle-run census is needed for these
  directed claims.
- Combined with the retained individual g=1,...,6 gates, all messages
  using one through sixteen actual groups have failure probability below
  2^-44 at 35% symbol distance. This does not cover larger source sizes.
  SMALL_STATE_EXACT_ROUTE_GATE.md records the exact rational witnesses,
  source factors, scope, commands, and directed acceptance rules.
- Explored source-cycle run counting as a floating control. It correctly
  weights rare clustered supports by their actual count instead of the
  full binomial. This helps diagnose difficult profiles, but the saved
  directed gates use the simpler full-binomial bound for every profile.
- Preserved two analytic refinements without adding encoder operations:
  average the four-block state-support emission floors in the direct
  operator; and retain quiet source constraints before the whole-state
  relaxation. The latter excludes source-zero paths. Because each outer
  column has rank at most 4g, its full-source cutoff needs only a Boolean
  local rank-deficit flag and a counter of 0..31 fully constrained
  columns, not a large rank-degree or threshold-union operator. Its
  output-zero tags require one common output marker before summing them.
  Neither refinement has yet been implemented or numerically certified.

Independent reviews covered the bridge, actual pair layout, source support
correction, fractional grouping, terminal convention, fixed verifier, and
generalization of its source counts. Root reviewed all executed source and
owned every numerical run. All jobs were serial, bounded by at most
55 seconds/1 GiB, and produced stdout only. The largest tuning run took
45.9 seconds and about44 MB; directed checks stayed below60 MB. Local and
remote storage remain above the20-GiB floor. No raw output files, caches,
new build trees, native kernels, benchmarks, paper edits, commits, pushes,
or promotion were made.

Next: extend the same exact-route gate to source counts above sixteen,
using bounded marker selection and directed fixed-witness replay. First
try this existing two-state interface; the inherited marker's failure
was not a structural obstruction at g16. If a genuine range remains
inconclusive, implement the Boolean full-column source cutoff rather
than another coordinate-support operator. Dense/intermediate coverage and
the final union over all 8192 source counts remain required. Only after
that proof gate should the small-state native implementation be timed
against the preserved 3.798-ms reference. The full goal remains active.

### 2026-10-21: exact-route closure through thirty-two source groups

Classification: further partial-proof progress, with no encoder change.

- Extended the fixed-witness verifier to g=32. Above g16 it requires an
  explicit rational output marker instead of extrapolating the old rule.
  Every g=17,...,32 now passes at 128 and256 bits: g17/18 below2^-300
  each, and g19,...,32 below2^-150 each. The new checks cover2784 profiles
  per precision, giving3224 for g=7,...,32. The combined failure bound
  for all messages on at most32 groups remains below2^-44 at35% distance.
  This still leaves g=33,...,8192 open.
- The change was fractional-marker selection, not added encoder mixing.
  The inherited g16 marker failed at g17. Rational witnesses with smaller
  alpha closed the new classes; one common (z,alpha)=(49/50,1/5) closes
  every g=23,...,32. The exact-route checkpoint records all witnesses.
- Added a bounded FLOAT continuation driver, preserving complete profile
  coverage at each published marker. Its sparse local table uses the
  identical exact rational provider formulas. A selected four-count run
  completed42 markers in17.12 seconds at about51 MB peak RSS. Floating
  proposals were separately replayed with directed arithmetic.
- Directed replays took at most5.03 seconds per source count, below59 MB.
  An initial4-second per-count cap stopped g31 before a result; its later
  complete replay used a10-second cap. The old g7 regression and the new
  explicit-marker guard also pass. No interrupted result is counted.
- Implemented the Boolean full-column source cutoff as an independent
  FLOAT control. The correct local fractional kernel is
  [[C^alpha,D^alpha],[0,(C+D)^alpha]], with entrywise powers before route
  averaging. Exact positive quiet partitions preserve the free-terminal
  branch and do not charge rank twice. The alpha1 regional partition and
  uncut counter controls agree with the old operator to rounding error.
  Selected clustered profiles at g17,24,32 gave no useful cutoff gain.
  This branch is preserved but is not needed for the directed new gates.
- A fresh independent audit found no gap in the projective fourth-root
  bridge, GF256-only source density, missing-group corrections, full-field
  DHD atom bound, or route fractionalization. Extension guards remain:
  h=0 requires g<256, and the full-source g8192 case must include its
  fully paired profile rather than apply the proper-cycle exclusion.

The next computational step is shared positive route-coefficient extraction
instead of repeating a polynomial for each matching count. Preserve exact
striped geometry and chronological multiplication. Larger classes still
need their full profile cover; isolated successful counts do not imply
intervals. Do not spend more time on the tested clean-column cutoff unless
a new profile shows that discarded source-zero paths actually dominate.
The construction, native library, paper, and preserved timings are unchanged.
No promotion, commit, or push is authorized. The full goal remains active.

### 2026-10-21: shared coefficients close source counts through sixty-four

Classification: further partial-proof progress in the same continuation.

- Added a shared triangular route engine: retain every singleton/full-pair
  coefficient with a+2b<=G, and square seven times. Initial coefficients
  also obey a+b<=32. The FLOAT version normalizes positive coefficients by
  their exact combinatorial masses; the directed version uses arbitrary-
  exponent Arb coefficients and one exact final normalization. Neither
  changes the route distribution or introduces an IID approximation.
- Nine floating regional comparisons at g7,17,32 agree with the old engine
  within 7.2e-15 in log2 entries. Directed g7,23,32 regressions also match.
  An initial directed smoke test exposed early dependency hashing before
  lazy provider imports; root moved hashing after evaluation and reran it.
- Every source count g=33,...,64 now passes below 2^-60 at both 128 and 256
  bits. The five rational witness windows are in the exact-route note.
  All 20672 new profiles pass per precision; total coverage g7,...,64 is
  23896 profiles per precision. The g64 worst logarithm is near -67.78357.
  The combined bound for every message on at most 64 groups remains below
  2^-44 at 35% symbol distance. No claim covers g65,...,8192 yet.
- Shared FLOAT G64 screens took about 1.0--1.3 seconds at 51 MB. Directed
  witness windows took at most 1.91 seconds at 60 MB. These are verifier
  costs, not encoder timings. Root owned all serial numerical jobs; agents
  only derived, reviewed, or edited disjoint source. All computations used
  bounded time/memory and stdout, with no caches or raw result files.

Next: extend the same shared proof engine to a bounded 65..128 continuation,
then replay any successful witnesses with directed arithmetic. Preserve
the source-correction and full-cycle guards before further enlargements.
If a class fails after reasonable marker selection, inspect its dominant
path/profile before adding encoder work. Full coverage and implementation
timing remain outstanding; the goal remains active.

A separate source audit of both shared engines found no normalization,
matrix-order, profile-coverage, terminal, or acceptance-reporting issue.
Final storage checks show 67.37 GB free locally and 66.89 GB in remote
scratch, both above the 20-GiB floor. The final process check shows no
remaining numerical job. No raw data, bytecode caches, build trees, native
encoder changes, paper edits, commits, pushes, or promotion were produced.

### 2026-10-21: run census closes source counts through 128

Classification: further partial-proof progress; the construction is unchanged.

- Extended the shared exact-route coefficient domain from G64 to G128.
  Normalizer conversion is safe in this bounded FLOAT domain; the largest
  observed normalizer has 945 bits. Initial local geometry is unchanged.
  Independent review retained the source-correction, proper-cycle, and
  positive-underflow guards. Larger domains are not silently supported.
- Full-binomial profile maxima stopped closing near g65. Their worst
  profiles were strongly clustered. Bridge Section 11 now partitions
  actual source sets by occupied cyclic runs, using exact integer counts
  outside alpha and a worst ordered profile within each class. This is
  valid with actual-empty pruning and adds no encoder work. Exact toy
  checks covered 5446 source sets and 90 classes on five small cycles.
- Directed replays at both 128 and 256 bits close every g65,...,128 below
  2^-150 per g. Witness windows are 65..77 at z191/200,alpha1/8; 78..100
  at z47/50,alpha1/10; and 101..128 at z23/25,alpha1/10. The verifier
  compares exact rational sums of certified Arb upper endpoints against
  the exact dyadic threshold. All 159104 new profiles are included per
  precision; cumulative coverage g7,...,128 is 183000 profiles.
- Combined with the direct gates, all messages on at most 128 actual
  outer groups have setup-failure probability below 2^-44 at 35% symbol
  distance. Actual sizes 129..8192 remain open. The new windows took at
  most 16.39 seconds and under 60 MB peak RSS. These are proof-verifier
  costs, not encoding times. Old-mode g7 regression still gives the prior
  -979.078576 bound; its optional census gives about -987.813450.
- Added a bounded FLOAT dense-anchor helper covering every matching count
  through IID histogram cells. It uses raw empty-input matrices so a
  successful anchor may extend to a range with actual-size source factors.
  Bridge Section 9 now also proves averaged-anchor monotonicity by double
  counting nested source sets. Neither fact is a directed dense receipt.
  The full-support FLOAT bound is about -526601 after controlling rare
  off-histogram empty batches. A single shared IID tilt fails at several
  intermediate anchors. Quiet and live entries favor different tilts;
  a positive entrywise envelope is the next discriminating test.

The proof, census implementation, and dense-helper interface received
independent read-only reviews. Root owned all serial numerical runs.
Only source and compact notes were retained. No native encoder changes,
paper edits, commits, pushes, or promotion occurred. The goal remains
active: finish the all-source proof before claiming this cheap proposal
reaches the preserved GF256 reference's encoding time.

The final dense control adds at most five fixed IID proposals and takes
their positive regional entrywise envelope. Complete g4096 coverage still
fails near +14748.49 at z370391/5732639,alpha2939512/14963605,tilts1,5/3.
The worst matching cell is 2040..2048 in both parities. A direct comparison
with the older selected-profile helper finds +29273.96 at the realizable
profile (2048,2047), versus -9543.45 at its previously selected
(2048,1024) profile. Raw-versus-pruned costs only 0.012 bits there.
The remaining dense obstruction is therefore not explained by raw
empty-state births or interpolation alone. The selected old grid missed
the strongly clustered two-parity class. Bridge Section 12 records the
complete cell bound; the diagnostic note records the failed controls.

Next discriminating test: isolate this clustered class and separate
off-histogram coefficient slack from local-kernel slack. Retain the
proved census and g<=128 checkpoint. Do not turn the old selected dense
successes or the full-support floating score into an interval claim.

Final checks: 67.36 GB free locally and 66.88 GB in remote scratch; no
numerical job remains running. The new helper and verifier were exercised
serially from their retained source. No raw data or bytecode caches were
created, and no benchmark ran during the proof work.

### 2026-10-21: first directed dense interval for the small-state proposal

Classification: verified partial-proof progress, with no encoder change.

- Isolated the previously failing g4096 one-run profile. Optimized IID
  tilts reduce its log2 bound from about +29274 to +4146.05 after exact
  run counting. A new positive-log coefficient control reduces it further
  to +1611.75. Coefficient conditioning matters, but does not close the
  inherited marker. These positive bounds do not exhibit bad codewords.
- The fractional parameter had not been bracketed. Keeping the same z
  and lowering alpha to 3/20 closes the one-run control strongly. A
  complete 121-cell anchor, with a single proposal tilt 3/2, then gives
  a log2 bound near -3247.12 for every source set at g4096. The earlier
  apparent dense obstruction therefore did not require more encoder work.
- Added a directed dense-anchor verifier. It checks exact integer striped
  occupancy masses, uses certified lower histogram probabilities, and
  takes exact rational maxima of certified regional upper endpoints.
  Bridge Eq. (36) bounds every region by one common matrix U. Its power
  U^256 covers all matching profiles without an extra union or any
  commutation assumption. Raw-operator monotonicity then extends the
  anchor while retaining each actual source size's density and count.
- Serial replays at 128 and 256 bits certify all messages on 4096..4108
  actual groups below 2^-100 combined, at 35% symbol distance. Both
  cover all 121 cells and all 13 source sizes. The log2 interval is near
  -1396.44545192524. Runtime is at most 1.81 seconds with 61 MB RSS;
  these are verifier costs, not encoding timings. The earlier g<=128
  checkpoint remains unchanged. Counts 129..4095 and 4109..8192 are open.
- Independent static audits passed both new helpers and the common
  regional-envelope argument. The exact-route control agrees with the
  shared g8/g16 engine within 9.4e-13 log2 units and with a constant-local
  fixture within 1.8e-15. The new integer occupancy weights agree with an
  independent floating implementation within 1.2e-16 for three complete
  distributions. A time-cap test fails closed. A complete g2048 control
  correctly rejects its positive bound near +25532.07.
- A conservative grouped-Jensen control for two, four, or eight batches
  gives no gain at the inherited dense witness. Preserve it as a failed
  diagnostic, not a reason to add inner operations. At g=2048, z=1/5,
  alpha=1/10 and 1/20 still fail; further marker or coefficient analysis
  is needed there.

The directed source hash is
`b0f0f3e01d6fe7cadb9a2b482791d4b7f9f35adf2c47b74f03d6be2618c968ec`.
The targeted FLOAT coefficient source hash is
`0a97176434433d9411fda9fc1b9adfb3d17f5e73b62d41bc9ba729238e4e33fb`.
Reproduction commands and fixed witnesses are in the exact-route gate.

Next: tune complete intermediate anchors, starting at g2048, and extend
successful raw anchors into certified intervals. Retain exact coefficients
for cells whose IID bound loses too much. The goal stays active: all-source
coverage and a competitive measured implementation are still outstanding.
Root owned every serial numerical run; agents only derived, edited disjoint
source, and reviewed. Only source and compact notes were retained. No native
encoder changes, paper edits, raw dumps, builds, commits, pushes, or promotion.
Final storage checks show 67.36 GB free locally and 66.88 GB remotely.
No numerical job remains running. The final documentation review retained
the partial-proof scope and clarified the regional averaging chronology.

### 2026-10-21: two-batch grouping closes an intermediate interval

Classification: verified partial-proof progress, with the encoder unchanged.

- Diagnosed g2048 with exact tail coefficients and a complete cell census.
  A retuned one-run class passes near -118.90, but the complete old census
  still fails near +3386.50. Passing a clustered example is not enough.
- Grouped two ordered batches before taking the fractional power. This
  retains the intermediate-state sum instead of paying for both branches
  separately. A positive ratio-bin calculation covers all 4845 occupancies
  without quadratic enumeration. The one-run control improves to -3478.61;
  the complete grouped census, not that example, is used for acceptance.
- Added the grouped FLOAT control and an independent directed verifier.
  Exact rational bin assignment, positive Arb aggregation, W^64, and one
  full regional histogram divisor implement bridge Section 14. Section 13
  records the complete cell census and averaged-anchor transport. It also
  proves mixed-marker transport, but the accepted witness uses common markers.
- Directed replays at 128 and 256 bits certify every message on actual
  group counts 2048..2051 below 2^-100 combined at 35% symbol distance.
  The fixed witness is z=23/100, alpha=19/200, tilt=33/25, pseudocount
  1/65536, bin width 1/4, maximum128 bulk cells, and32 singleton tail cells.
  Both replays cover157 cells,24648 compatible ordered cell pairs,2048
  run classes, and all four actual source sizes. Their log2 bound is near
  -125.5887508163504. The maximum verifier cost is4.43 seconds and79 MB;
  these are not encoder timings. Source hash:
  `bf9f75a0f8f3984cd53edafbc84095b5af8d27ef6e27e8b551f40909ef5bc83e`.
- Cross-checked the FLOAT grouped kernel against all23474025 ordered
  profile pairs for one proposal. Every binned entry is an upper bound;
  alpha1 recovers the one-batch square within1.43e-13. All1809 exact bin
  boundary tests pass. A directed noncommuting g8 census fixture agrees
  with direct ordered-profile enumeration within2.6e-38 relative error.
  A complete positive-bound control rejects, and an immediate timeout
  produces no bound or passing claim. Independent static audits passed.

Current directed coverage is g1..128,2048..2051,4096..4108. Counts
129..2047,2052..4095,4109..8192 remain open. The goal stays active:
the complete all-source proof and competitive native timing are not done.
Next evaluate a complete g1024 grouped anchor, starting with retained
z=2825730/5340629, alpha=1/10, tilt=33/25, and the same bin/cell settings.
That is an untested starting witness, not an extrapolated certificate.

Root owned every serial numerical run. Agents derived, edited disjoint
source, and reviewed without executing. Only source and compact notes were
retained; no encoder changes, raw dumps, builds, paper edits, commits,
pushes, or promotion. Fresh checks show more than60 GiB free locally and
remotely. No numerical job remains running.

### 2026-10-21: wider directed interval and a longer-window witness

Classification: verified partial-proof progress, plus a new complete FLOAT
witness. The construction and implementation are unchanged.

- Added fixed-parity reuse of unconditioned regional numerators. Recondition
  the complete matching-count partition and count source runs at each actual
  size. This avoids the larger source factor/smaller route mismatch in raw
  subset transport. Bridge Section 15 records the proof.
- The new directed verifier closes 2048..2073 at 35% distance, below 2^-100
  combined. Both 128- and 256-bit replays pass; their log2 union bound is
  approximately -276.66625679167865. Each covers 26 sizes and 53,573 run
  classes with complete 157-cell parity families. Maximum verification cost:
  10.06 seconds and 86 MB. Source hash:
  `79644db7e509a017340918e847b0685f51fc8910c72edb84bcbbf7fd385711ff`.
  Extending the same families through 2074 produces a complete positive bound
  near 43.0572 and correctly rejects. An immediate timeout makes no claim.
- At g1024, a marker bracket lowers the complete two-batch bound to +2505.27,
  but does not close it. A new projective-row calculation keeps 4..128
  consecutive batches together. It propagates one rounded ratio per
  history and retains all 4845 occupancy profiles. At z=1/2, alpha=3/40,
  tilt=6/5, pseudocount=1/65536, tail count 32, and grid ratio 9/8, its
  complete 153-cell scores are +812.67 for 4 steps, -35.57 for 8,
  -450.52 for 16, and -758.88 for 64.
  The 16-step witness is the next directed-replay target; none of these
  FLOAT values is added to certified coverage.
- Rank-one and exhaustive noncommuting short-window fixtures check the
  implementation and the analytic quantization-loss bound. A proposed
  common weighted functional remains too loose, near +12582 bits. Its
  supremum forgets the evolving state ratio. Preserve this failed control
  rather than treating it as a need to increase encoder work.

Current directed coverage is g1..128, 2048..2073, and 4096..4108. The open
ranges are 129..2047, 2074..4095, and 4109..8192. No all-source certificate or competitive
timing exists for this small-state proposal. The goal remains active.
Next: a bounded directed 16-step replay at g1024, followed by direct interval
extension if it passes. Preserve the 3.798 ms GF256 performance reference.

Root owned all serial numerical jobs. Agents only derived, edited disjoint
source, and reviewed. Only source and compact notes were retained; there
were no raw dumps, builds, encoder changes, paper edits, commits, pushes,
or promotion. Fresh checks show 67.04 GB free locally and 66.88 GB remotely.
No numerical job remains running.

### 2026-10-21: park small-state work and close the window design at 49.1%

Classification: complete finite-proof progress for the user-selected
high-distance construction. The 35% objective is superseded.

- Saved the small-state specification, partial coverage, source hashes,
  FLOAT-only window witness, and proof-method ceiling in
  [the parked checkpoint](PARKED_SMALL_STATE.md). No sources were deleted
  and no partial result was promoted to a whole-code certificate.
- Returned to fixed RS[256,128] over GF(2^128), one stratified regional
  permutation, and uniform full-field Toeplitz window feedback. Preserved
  the old >45%/72-bit certificate and its sources unchanged.
- Identified an exact Jensen obstruction to certifying 49% with the
  old memory 256/batch 64 envelope. At c128,v16512 the forced-tail gain
  is smaller than the mean shortened-dimension penalty. This is a bound
  obstruction, not a low-distance codeword. An exact integer check passed.
- Increasing memory and batch to 1024/256 retains the four-batch window
  architecture. Added a bounded wrapper using the existing entropy volume
  and finer rational markers. All retained numerical helpers are unchanged.
- Closed >49%/40 bits with 851 boxes, then >49.1%/40 bits with 926 boxes,
  without another encoder change. Each cover includes all 6714983 feasible
  activity-size/count pairs. Separate 80- and 120-digit directed replays
  pass; an independent integer adjacency check establishes complete,
  nonoverlapping coverage. Removed, duplicated, and boundary-corrupted
  covers reject. An immediate timeout produces no certificate.
- [The near-half note](NEAR_HALF_WINDOW.md) records exact thresholds,
  source and partition hashes, reproduction commands, first rejected
  probes, and the cost boundary. Maximum final-source replay time was
  20.05 seconds, not encoder latency. Sampled coefficients per output
  change only from 319/64 to 1279/256, but state and batch are four times
  larger and no runtime equivalence is inferred.

Independent agents reviewed the larger-batch proof interface and wrapper
without executing. Root owned every serial numerical run. Complete-cover
and replay jobs were limited to 50 seconds internally, 55 seconds externally,
and 1 GiB virtual memory; selected-point diagnostics used external caps.
Only source and compact notes were saved. No raw dumps, new native builds,
paper/library changes, commits, pushes, or promotion. More than 60 GiB
remains free locally and remotely. No numerical job remains running.

Next: measure the actual high-distance Toeplitz profile and optimize its
middle products without changing its certified distribution. Near-50%
distance is the objective; do not resume the parked 35% coverage search.

### 2026-10-21 — Native near-half implementation and coordinate variant

- Implemented the exact 256-by-1024 Toeplitz map using twelve small middle
  products, then an exact nine-product block factorization. Cached state
  transforms and combined output spectra need only two forward/two inverse
  transforms per batch. The prepared tape falls from 384 to 288 MiB.
- The first exact control measures about 154 ms at K20. Deferred arithmetic
  and the nine-product decomposition reduce that to about 130 ms. This is
  much slower than the retained RS reference; the proof alone did not imply
  an efficient implementation.
- Introduced a separately labeled coordinate variant, T'=diag(U,U) T diag(Vx8).
  Its fixed invertible maps remove runtime monomial conversions. Independent
  review checked that the original fixed-output inverse distribution, and
  hence the >49.1%/40-bit finite certificate, applies unchanged. This is not
  claimed to be identical encoding or a weight-preserving transformation.
- The fastest tested coordinate variant uses deferred reduction and exact
  raw GF256 scaling. Repeated 31-sample K20 medians are 70.90 and 74.30 ms.
  Its final phase run has 6.89-ms outer and 66.58-ms routed inner. The
  standalone 5.82-ms routing time must not be added to that fused inner.
- Checked all 1024 state-basis inputs, dense independent field products,
  diagonal boundaries, coordinate conversion by independent root products,
  full recursion through ring wraps, both routing branches, and arithmetic
  variant equality. UBSan passes for the exact control and fastest variant.
- [The native note](NEAR_HALF_NATIVE.md) retains construction details,
  certificate transfer, source interfaces, checks, sample counts, timings,
  and build commands. No proof helper or existing encoder was changed.

Root ran all builds, checks, and timings serially with explicit time/memory
limits. Agents derived, implemented disjoint headers, and reviewed without
running numerical work. No raw data, random tapes, generated archives,
paper/library changes, commits, or pushes. The unresolved performance cost
is now primarily 18.9 million full-field products per codeword, not basis
conversion. Next: seek fewer or cheaper products while preserving the
implemented inverse-law proof; do not infer a near-5-ms encoder from this result.

### 2026-10-21 — Phase diagnosis against the 3--5 ms target

- Added two measurement-only native files, leaving the retained encoder,
  field headers, and certificate sources unchanged. Split outputs agree
  byte-for-byte with the frozen encoder at K12 and K20.
- K20 medians: original 74.319 ms, untimed split 74.332 ms, instrumented
  split 75.255 ms. Instrumented phases are outer 7.788 ms, pointwise
  field work/keys 55.280 ms, inverse transforms 3.650 ms, routed output
  4.835 ms, state append/transforms 3.470 ms, residual 0.203 ms.
- Cache-hot reuse of one prepared key still takes 49.474 ms for the
  isolated pointwise stage, versus 55.765 ms streaming keys. This is not
  a replacement construction. Arithmetic remains a major obstruction.
- Identical RS arithmetic with contiguous tiled output takes 3.547 ms
  versus 6.984 ms with regional stores. Values agree after reindexing;
  the resulting routing-locality tradeoff has not been measured end-to-end.
- The prepared 288-MiB tape alone reads in 6.870 ms. Both setup-material
  traffic and arithmetic must fall substantially; microtuning one does
  not meet the whole-encoder target.
- The native note records measurement overhead, nonadditive diagnostics,
  the fixed-family bilinear-rank bound, and an explicitly aspirational
  phase budget. The proof-side next test is whether outer constraints
  permit a much cheaper, rank-dependent feedback law, rather than joint
  uniformity for every arbitrary nonzero state. No new certificate claimed.

Root owned the serial bounded runs; agents derived and reviewed without
executing. No library/paper changes, raw dumps, commits, pushes, or promotion.
The target remains near-half symbol distance, one structured route, and
GF(2^128)-linearity. Next: jointly choose a cheap feedback law and a proof
that accounts for its actual deficient-rank trajectories, while testing
exact outer-layout improvements against the full encoder.

### 2026-10-21 — Optimize the exact small-RS path before changing the code

- User constrained this pass to implementation changes, benchmarking the
  outer alone and outer plus routing. No outer, route, inner distribution,
  field representation, or certificate change was made.
- Added `native/RsOuterOptimized.h` and `native/rs_outer_bench.cpp`, retaining
  the frozen implementation. The selected path fuses source loads and zero
  handling, batches 128 blocks, and fuses register-resident final FFT layers
  with coordinate-major output stores. Routing uses batched loads after a
  sequential cache-fill pass. No hot-path allocation or indirect dispatch.
- Two matched 31-sample huge-page runs reduce the outer 8.297/8.258 to
  3.303/3.245 ms, and outer plus routing 12.372/12.234 to 5.923/5.930 ms.
  A matched ordinary-page run gives 7.951 to 3.554 ms and 12.778 to 6.519 ms.
  The driver verifies actual huge-page backing and alternates case order.
- Rejected tiled output and non-temporal stores despite favorable outer-only
  timings: both lose after routing. Register leaves without fixing the
  output-store schedule are insufficient.
- Exact output comparisons, position/field-basis tests, an independent RS
  generator oracle, non-power-of-two block tails, input preservation, guards,
  and K12/K20 UBSan pass. The diagnostic includes an exact full-encoder path
  with only the outer replaced; standalone routing improvements are not
  assumed to transfer to the fused inner.
- Final-source K20 confirmation gives 3.329 ms outer and 5.952 ms outer plus
  routing against 8.312/12.366 ms controls. The integrated full encoder gives
  69.588 ms against 74.406 ms, with every output matching. It uses the old
  fused routing and inner. Final-source ordinary-page results are 3.585 and
  6.517 ms against 8.184 and 13.602 ms; the native note retains exact scope.

All builds and numerical jobs were serial, with explicit time/memory limits
and more than 20 GiB free. No raw dumps, proof changes, library/paper edits,
commits, or pushes. The earlier outer timing was implementation/layout cost,
not evidence that the small-RS construction is intrinsically expensive.
Next: audit the inner's implementation with the same standard before making
construction changes or claiming that its measured runtime is unavoidable.

### 2026-10-21 — New goal: purpose-built coded-state inner

The user explicitly requested a new inner design, with outer and permutation
fixed, rather than further tuning the certified inner. The active goal is
>=49% all-message symbol distance with >=40-bit setup-failure margin at K20.
The [design note](coded_state/DESIGN.md) records the actual setup distributions,
local proofs, failure modes, countermeasures, and remaining global problem.

- Constructed fixed GF256 maps A: F^64 -> F^128 and B: F^128 -> F^64
  from interpolation on disjoint Cantor cosets. BA=I; image and kernel
  have distance65. Stronger combined input/output bounds follow from
  [192,64,129] and [192,128,65] RS codes. Small additive transforms implement
  the maps; dense matrix multiplication is unnecessary.
- Identified a one-batch reset in Q_t=B Y_t. Replacing it with
  Q_t=Q_{t-1}+B Y_t retains state across zero output at one extra state XOR.
  Excluding diagonal entries0,1 also preserves support on empty input.
- The one-tap inverse distribution has only r random coordinates for a
  state of weight r. A naive intersection/support bound loses too much
  information; this is not a code counterexample or a closed certificate.
- Zoomed out to a two-tap alternative with complementary A,A2 and
  BA=BA2=I. It uses one full-field product/output, versus nine in the
  certified reference, and shares its 128-point transform. A full-support
  state gives uniform feedback conditioned on64 affine inequalities.
  Partial-support trajectories remain the unresolved global proof case.
- Added a compact exact checker, without touching frozen helpers.
  Two serial bounded remote runs passed. The final run (1.712s) checks
  65,535 nonzero toy states/kernel words, all494 toy B minors, actual
  GF256 maps against independent Lagrange interpolation, both recurrence
  inverses, complementary-map rank, and two-tap support preservation.
  Limits:60s/1GiB, >20GiB free. No benchmark, native timing claim, dumps,
  library/paper edits, commit, push, or promotion.

Next: derive a joint outer/trajectory enumerator for the two-tap recurrence,
especially how partial state support is created and maintained. Retain
one tap as the lower-cost comparison. Do not confuse dense-state uniformity
with a whole-code proof, or multiplication counts with end-to-end timing.

### 2026-10-21 — Two-tap partial-state structure and a closed global class

The previous goal turn was progress: it created proved local maps, exact
checks, and two candidate recurrences. This continuation preserved the
outer and route and advanced the actual two-tap proof.

- Proved that a state of weight r gives two independent linear RS[64,r]
  feedback spaces before the affine sampling exclusions. Prescribing z_L
  and z_R output coordinates has exact rank min(r,z_L)+min(r,z_R).
- Added a lower bound for joint output and next-state observations, without
  assuming those observations independent. Exhaustive GF16 checking covers
  20480 cases; the bound is strict in1548, leaving a concrete refinement.
- Reparameterized the SAME distribution into independent nonzero state
  multipliers H and unrestricted kernel multipliers E. State evolution is
  Q'=BX+H Q. Conditional output randomness is an affine coset of an
  [128,r,129-r] MDS code. This separates state thinning from output counting.
- Closed an all-message, arbitrary-route trajectory class using an exact
  integer bound: at49%, total state deficit M<=2293 gives failure <2^-255;
  at49.1%, M<=1245 gives <2^-127. The mandatory initial empty state counts
  toward these budgets. These are NOT whole-code certificates.
- Removed whole inactive regions using exact shortened outer dimensions.
  Taking all32896 aligned intervals gives <2^-239 and <2^-111 for the same
  internal deficit budgets. Terminal state is free when the interval reaches
  the end. Arbitrary partial-region endpoints and larger deficits remain open.
- Added compact source checkers and PARTIAL_PROOF.md. Remote runs were
  serial with30/60s time and1GiB memory caps, >20GiB free; no large output,
  benchmarks, changes to retained controls, or performance claims.

Next: combine the separated state transition and MDS coset enumerator with
outer activity/shortening. A single input activates all state coordinates,
but batches with >=65 active inputs may lie in ker(B). Bound long partial
state episodes and nonaligned intervals jointly with the outer; do not
silently extrapolate the bounded-deficit class to all messages. The goal
remains active, and neither candidate has been promoted.

### 2026-10-21 — Structural countermeasure: randomize after expansion

The prior turn was progress: joint rank lemmas and finite all-message
trajectory-class bounds were proved and checked. This round reassessed
the non-coordinate intersection obstruction instead of only tuning it.

- Added the post-expansion diagonal candidate Y=X+D*A*Q, Q'=Q+B*Y, keeping
  the outer, route, and existing small-field maps unchanged. It has one
  full-field product/output, like two taps, but one fewer64-point transform.
  No native timing is assigned.
- Proved exact uniform next state in F^64 for EVERY nonzero incoming Q:
  A*Q has >=65 nonzero coordinates and B restricted there has rank64.
  The tradeoff is an empty-input reset probability q^-64, which must be
  counted globally. Zero output still preserves the state.
- Proved an exact coordinate inverse distribution. The existing outer
  shortened-dimension formula now applies directly, removing the two-tap
  non-coordinate intersection issue. Enumeration of masks remains open.
- Closed bounded-hole outputs at49% with H<=4587 below2^-128. Separately,
  full-scan live trajectories with <=31 nonzero outputs per batch have
  failure <=2^-1044480, even with large persistent hole budgets. Neither
  result covers the full bad-output set.
- Found a precise proof pitfall: AB has8064 zero entries, and single-coordinate
  output trajectories can sustain expanded states of weight65. Repeated holes
  do not automatically impose fresh equations. A fixed GF256 row scaling
  of B preserves the local theorems and reduced the diagnostic zero count
  to65, but was not selected; fewer zero entries do not prove better distance.
- Added POST_DIAGONAL.md and one bounded source checker. Exhausted65535
  toy states/all93 expanded masks,48 exact refresh fixtures, and checked
  actual GF256 minimum-weight expansion, full refresh rank, and inverse.
  All numerical runs serial,60s/1GiB caps, >20GiB free. No benchmark,
  raw data, production edit, commit, push, or promotion.

Next: prioritize the post-diagonal candidate's mask/interval enumerator
against actual regional moments, comparing to two taps as a retained control.
Treat empty-state intervals and correlated holes explicitly. Add map scaling
only if it removes a demonstrated global loss. The >=49%/40-bit full-code
goal remains active and unachieved.

### 2026-10-21 — Complete sparse-source certificate through31 blocks

The previous turn was progress: the post-diagonal countermeasure and its
exact uniform-refresh/coordinate-inverse lemmas changed the proof strategy.
This turn kept that construction unchanged and closed an actual source class.

- Proved a joint-in-time tail bound for TOTAL holes of nonzero refreshed
  states, N^m*q^-m, rather than charging63 holes for each exceptional batch.
  Reset and random-output-zero tails are separate union events; no false
  independence assumption is used.
- For p<=64, every nonzero input batch activates an empty state, because
  its weight is <=p and ker(B) has distance65. This gives an explicit episode
  cover of all nonzero inputs, with free terminal state.
- Counted projective messages by active outer blocks p and occupied regions
  h. The shortened source dimension is p(h-128). A deterministic regional
  span bound then converts few resets/defects into a weight lower bound.
- Exact integer replay covers all3968 cells p=1..31, h=129..256. Every
  nonexceptional output has weight >=1039505, exceeding floor(.491N)=1029701.
  The union of all three exceptional-event families is <2^-50. This holds
  for every regional route realization and covers setup-adaptive messages.
- Added SPARSE_PROOF.md, sparse_certificate.py, and an independent small-
  timeline/coverage test. All9580 toy timeline summaries and all3968 integer
  cells pass. The final state may be nonzero. Runs were serial,30s/1GiB caps,
  >20GiB free, with no data files or benchmarks.
- First uncovered p=32 has a failed envelope at h=256: r=64, m=4906,
  L>=8063, weight lower bound1018159. This is NOT a counterexample. The span
  bound allows resets to skip large gaps without using regional randomness.

Next: use actual regional occupation/shortening to limit the inactive gaps
for larger supports, rather than change the construction prematurely. Beyond
p=64, nonzero input in ker(B) is an additional branch requiring proof. The
full >=49%/40-bit objective remains open; no promotion, code-library change,
paper edit, commit, push, or native-performance claim was made.

### 2026-10-21 — Regional episode certificate through64 active blocks

The previous turn was progress: every p<=31 source class closed with exact
coverage. This turn used regional sampling to resolve the next range without
changing the construction, maps, or setup distribution.

- Derived a joint prescribed-reset/defect bound. Refreshed nonzero-state
  holes cannot occur at reset steps; output-zero witnesses exclude those
  batches. This permits chronological multiplication of reset and defect
  costs without assuming the two defect families independent.
- Counted projective words in the outer shortened to an episode union U,
  before averaging the route. Positive per-block occupancy moments retain
  the requirement that every selected source block has >=129 positions in U.
- The first probe improved p=48 substantially but did not close all selected
  p=64 points. Splitting terminal and nonterminal episodes fixed that loss:
  terminal intervals have one fewer partial regional boundary; nonterminal
  intervals require an additional reset. No code changes were needed.
- Complete coverage for p=32..64 and all feasible episode lengths/counts in
  both terminal classes:15709 boxes,31416 checks,4429185024 discrete cells.
  Exact rational route factors and outward replays at80/120 digits pass.
  Worst log2 upper bound is -80.173269845194262..., giving <=2^-47 after
  the full cell union. Combined with p<=31, failure is <2^-46 at >49.1%.
- Independent checks:128 enumerated interval-count cases,2500 exact positive
  route-moment comparisons, and disjointness on66 cover slices plus a separate
  count of the full domain. All pass. A large-marker floating overflow was
  fixed in a separate evaluator; retained proof helpers were not edited.
- Sources:ROUTE_EPISODES.md, route_episode_probe.py,
  route_episode_certificate.py, test_route_episodes.py. Runs serial, capped
  at60/200s and1GiB, >20GiB free, no persisted cover/report or raw data.

Next obstacle is genuinely different: p>=65 permits nonzero input batches
in ker(B). They emit >=65 nonzero symbols but do not activate the state,
so episodes no longer cover all input. Audit that quiet branch and its
outer compatibility before selecting a new cheap countermeasure. The full
all-message >=49%/40-bit goal remains active and unachieved. No benchmarks,
library/paper changes, commits, pushes, or promotion.

### 2026-10-21 — Small-field column scaling closes the quiet branch

The previous turn closed all messages on <=64 outer blocks. This turn
attacked nonzero inputs that fail to activate an empty state, rather than
extending that proof by a loose source-count sweep.

- Added a separate inner-only variant: Y=X+D A Q, Q'=Q+B R Y, with each
  R entry sampled independently from GF256 minus zero. The fixed outer,
  regional route, A, and B are unchanged. Cost is one additional GF256
  scaling per output and 2MiB raw coefficients; no extra full-F products.
  Runtime and prepared storage are not measured.
- The <=64-block certificates transfer for every fixed R, since kernel
  distance65, column independence, and exact uniform state refresh persist.
- For trajectories with zero incoming state on every batch, all constraints
  are over GF256. Scalar extension preserves this subcode's minimum symbol
  distance exactly. The last state is free, including final-batch activation.
- Each nonzero quiet GF256 batch costs <=255^-64 under column scaling.
  Combined actual-route shortening moments and a fractional union bound
  cover all15714 possible nonzero quiet-batch counts, across all messages.
  Four boxes suffice; exact rational moments and directed80/120-digit
  replays pass. Worst upper log2 summand=-11327.1404373998...; conservatively
  retain2^-80 per cell and <=2^-66 for the entire branch at >49.1% distance.
- Independent checks:15 small-field quiet fixtures,432 refresh fixtures,
  a nonvacuous exhaustive GF4/GF16 scalar-extension test with free terminal
  state, the integer coefficient comparison, and complete count coverage.
- Added coded_state/QUIET_BRANCH.md, quiet_certificate.py, and
  test_quiet_branch.py; indexed the distinct variant. Runs serial, capped
  at160s/60s and1GiB, >20GiB free. No data files, benchmarks, library/paper
  changes, commits, pushes, or promotion.

Remaining obstruction: large-source messages can mix live episodes with
quiet nonzero kernel batches. The all-quiet descent does not apply once
full-F D coefficients participate. Next retain the common message-space
dimension after quiet constraints and combine it with live-episode costs;
do not multiply independently averaged partial certificates. The full
near-half/40-bit goal remains active, with this turn making genuine progress.

### 2026-10-21 — Common quiet rank and a dense mixed-class obstruction

Previous turn: progress, with the entirely quiet branch closed. This turn
kept the encoder unchanged and derived a composable common-space interface.

- For fixed route and allowed input coordinates, quiet constraints have
  matrix L diag(R) C over GF256. A common independent coordinate set of
  size r reduces a minor to random nonzero diagonal plus a fixed matrix.
  Sequential pivots bound rank deficit by Bin(r,1/255), giving an explicit
  moment for the remaining message dimension over both GF256 and GF(2^128).
- MDS shortening makes the common independent set a capacitated bipartite
  matching between outer blocks and quiet batches. Source capacities are
  shortened dimensions, batch capacities64, edges actual routed coordinates.
  This is an exact structural-rank interface, not a generic-route assumption.
- Conditioning on R before counting F-projective messages and exposing D
  gives the combined bound in coded_state/MIXED_RANK.md. Quiet batches pay
  their65-symbol output contribution before setting the live defect threshold.
- Important negative result: v=j=8192, one terminal episode, full allowed
  support defeats this forward envelope EVEN with perfect quiet rank and no
  route/pattern penalties. Residual dimension>=524288, T=14743637/64, leaving
  positive q-exponent4702683/16. This is a bound failure, not a codeword.
  It changes the next action: do not spend another round only improving rank.
- Exact tests pass:704 small-field diagonal tail comparisons,5184 rectangular
  product ranks with rational moments,4096 max-flow/brute-force comparisons,
  and the exact K20 obstruction arithmetic. Serial run,60s/1GiB caps,53GiB
  local free and63GiB remote free; no output files, benchmarks, or commits.

Next use output-side enumeration for dense mixed classes, retaining quiet
kernel ranks and the dependence of feedback holes across zero-output runs.
The rank interface remains useful for small allowed supports. No full
near-half certificate or native timing was obtained; goal remains active.

### 2026-10-21 — Output count closes the prior dense mixed obstruction

Previous turn made progress by isolating a failure of the forward count.
This turn changed the counting direction, not the encoder or its setup.

- For fixed output and compression scalings, inverse randomness has support
  S contained in the live intervals A_live, with |S|>=65(v-c). Bound its
  affine outer intersection using d_C(A_live), avoiding a separate union
  over internal feedback-hole masks. Quiet output batches have kernel
  rank64 on every possible nonzero support, giving2^u*q^(W-64j) outputs.
- Derived the conditional output bound and its route-fractional version
  in coded_state/OUTPUT_INTERVAL.md. It counts setup-adaptive messages via
  fixed output candidates, without conditioning on future diagonal events.
- Full integer coverage:252978 feasible pairs for one episode of at most
  8192 batches and at most32 outside zero batches. Every specified pattern
  contributes <=2^-44672; all placements cost <482 bits, giving <=2^-44190
  at >49.1%. The bound holds for every route and every nonzero R. The
  formerly failed aligned dense example now has log2 bound -1359104.
- Added a deterministic output-rank interface for state-hole constraints:
  new variables in batch t contribute rank >=max(0,h+min(64,w_t)-64),
  or min(64,w_t) when the next state is zero. Fresh columns make these
  ranks additive; no constraint is placed on the final state. Zero-output
  runs preserve state and therefore preserve the hole pattern exactly.
- Exact checks:512 interval-region cases,65535 toy outputs/129 trajectory
  groups,23808 local submatrix ranks,512 joint state-rank fixtures, and
  all integer cover cells. Local bound strict7604 times; h=w=3 toy minors
  include552 rank2 and2584 rank3, so neither generic-full-rank nor uniform
  worst-rank treatment is sharp. Sources:output_interval.py and its test.
- Serial checks,60s/1GiB caps,>20GiB free, no data files or benchmarks.
  No library/paper changes, commits, pushes, or promotion.

Next: retain AB's actual half-block structure and exceptional-minor counts
in the output/state-hole enumeration for longer and fragmented trajectories.
Combine this with the interval envelope; do not independently charge holes
that persist on zero-output runs. Full near-half/40-bit goal remains open.

### 2026-10-21 — Alternating evaluations and live-to-live rank checks

The user resumed after the prior pending local test. Replayed it from
source after confirming no numerical process remained. Classification:
new local proof and structural candidate, not a whole-code certificate.

- Alternate A between the two disjoint128-point evaluation halves of
  GF256. Each evaluates the same degree-<64 state polynomial; B, the
  outer, regional route, and D/R distributions stay fixed. This adds no
  full-field products or randomness. Native cost is not measured.
- Proved the pair distance193 for unchanged nonzero state, and the run
  hole cap63*ceil(length/2). The previous fixed-A bound allowed63 holes
  per batch. All established sparse, quiet, and output-interval partial
  certificates transfer analytically through unchanged local hypotheses.
- Derived a joint two-observation rank lower bound. Eliminate the old
  state before charging constraints to new outputs; do not charge the
  state rank twice. The exact projected rank is joint rank minus the
  rank on old-state columns. Recorded proof in ALTERNATING_MAPS.md.
- The existing half-count diagnostic has449 profiles and two ambiguous
  ranks. The new map also has deficient minors: GF16 A1B rows(0,4),
  columns(2,6) have rank1. Generic-full-rank reasoning remains invalid.
- A more informative exact test separates zero-state fibers. Among
  50176 maximal-hole/one-output fixtures,896 are rank deficient, but
  every such fixture forces either activation or reset. None has two
  live endpoints. Exhausting1048544 projectively normalized live pairs
  gives minimum combined feedback weight11 in both phases, versus13
  for unchanged state. This is GF16 evidence, not a state64 theorem.
- Other exact checks pass:65535 toy states with run parities,864 refresh
  fixtures, actual GF256 FFT/Lagrange and inverse checks,17768 square
  AB minors,4096 joint-rank fixtures, and two eliminated-state fiber
  counts. The first two-observation bound is strict in239 fixtures.
- Sources are alternating_maps.py, test_alternating_maps.py,
  test_alternating_rank.py, and half_rank_probe.py. Serial runs stayed
  under60s/1GiB with >20GiB free; no output files or benchmarks. No
  production/paper change, promotion, commit, or push.

Next derive a sparse-update split enumerator in the actual maps, retaining
zero-state fibers before applying a worst rank. Then combine it with the
inverse route/shortening bound over complete mixed classes. Alternation
resolves persistent zero-output holes, not arbitrary changing-state paths.
The full >=49%/40-bit target remains unachieved. The app goal reached a
usage limit before this explicit user-requested continuation; no completion
or new performance claim is made.

### 2026-10-21 — Actual-size sparse-update theorem and restricted extensions

The user requested further exploration and internet research. The outer,
route, and alternating-map candidate are unchanged in this iteration.

- Proved the exact compression image for output supported in one interpolation
  half: d=g_(C_e minus W)*v, deg(v)<|W|. Root counting on the two other
  cosets yields a live-state feedback-pair floor of 150/256 for weight-one
  output updates, and 149,146,141,135 for weights 2,4,8,12. This is a proof
  over GF(2^128) for the actual state64 maps, not toy extrapolation.
- Derived exact zero-state/zero-update fiber subtraction and a dimension
  cap. Allowed support is distinguished from exact support; shared state
  boundaries cannot be independently resampled or charged twice.
- Combined these bounds with output enumeration. Three complete restricted
  one-episode classes now close at >49.1%, always with at most32 outside
  zero batches: v<=8700 with >=90% weight-one canonical updates; v<=8800
  with all weight-one canonical updates; v<=10000 with >=75% zero canonical
  updates. Integer coverage includes 269742,273042,312642 cells respectively.
  Their full placement bounds are 2^-706206,2^-770206,2^-72734. These very
  large margins apply ONLY to the stated overlapping classes, not the code.
- Final-source checks pass: 524272 exhaustive GF16 live pairs, six exact
  zero-fiber counts and dimension caps, 304 GF256 structured-hole fixtures
  across state sizes4..64, and all restricted-class integer cells. An
  independent JavaScript calculation agrees with the cell counts/maxima.
  Serial checker runs took about15.6s and0.15s, not encoding benchmarks.
- Zoom-out obstruction: even all zero updates in the canonical pairing
  do not close the current bound for a full-length episode. Its feedback
  floor1580928 and containing-interval dimension1048576 leave positive
  field exponent497349 before support counting. This diagnoses the coarse
  outer-shortening relaxation, not a bad codeword or failed construction.
- Consulted primary MDP and weighted-RS convolutional-code literature.
  Their multitime-minor criteria support checking more than block MDS, but
  their existence/field-extension hypotheses do not certify our maps.
  Links and scope are recorded in coded_state/HALF_UPDATES.md.
- No native benchmarks, generated reports, production/paper changes,
  commits, pushes, or promotion. Tests were bounded to60s/1GiB, with the
  last storage check showing34GiB local and63GiB remote free.

Next retain actual feedback support S, instead of its whole interval,
in the regional outer-shortening moment. Compare this saving against a
joint output/state-hole rank count on the full-length stress class before
expanding parameter sweeps. Dense cross-half updates and multiple episodes
remain required. The >=49%/40-bit whole-code goal is unachieved. The app
goal remains usage-limited; the requested replacement was rejected because
that unfinished goal exists. This authorized turn continued regardless.

### 2026-10-22 — Actual feedback support closes a full-length class

The user requested sharper routing with the construction fixed. No outer,
route, inner map, setup distribution, or encoder implementation changed.

- Retained S(y,R), the actual feedback support, before route averaging.
  For fixed output and R it is independent of route/D. Exact rate-half
  complementation changes the inverse exponent to -K+d_C(S^c), instead
  of using the entire live interval. Derived its nonidentical regional
  Poisson-binomial moment and a convex equalized binomial majorant.
- Kept mask multiplicities OUTSIDE the fractional power. The new argument
  conditions on a fixed actual mask, counts its outputs, fractionalizes
  that event, and then sums masks. It does not average one fixed mask's
  route moment over a setup-adaptively selected support without a union.
- Closed the full-length live class Q_1,...,Q_(n-1)!=0 with
  y_2=y_4=...=y_(n-2)=0; the terminal state and last output remain free.
  Paired states have at most63 roots across256 points. There are8191 such
  pairs and one128-point final observation, hence516097 possible hole totals.
- Actual routing alone was insufficient in the first broad diagnostic.
  Added fresh-output constraints on the holes and the exact local tradeoff
  w-r(h,w)<=128-h-(64/127)(128-w). Triangular fresh columns justify adding
  these ranks; incoming state dimensions are not counted twice. Explicit
  truncated mask polynomials and a weight-aware output-volume bound remove
  the remaining loss. The initial FLOAT highest-hole score was positive;
  the revised bound passes with no new encoder operation.
- Directed80/120-digit final-source replays cover every H=0..516096 in37
  intervals. Every summand is <=2^-80, giving full restricted-class failure
  <2^-61 at >49.1% symbol distance. Weakest box log2 is about -232.10.
  These are complete mask/output/route counts for this class, NOT a
  whole-code certificate or a claim about arbitrary intervening outputs.
- Exact tests pass:4096 incidence matrices,768 route-moment comparisons,
  512 multitime fresh-column rank fixtures, all scalar tradeoffs through
  state64, small output-volume domains, directed rational operations and
  independent rational-series checks of logarithm enclosures. Immediate
  timeout rejects without producing a certificate. Final directed runs
  take about2.5--2.8s, not encoder latency.
- Sources and proof: coded_state/FEEDBACK_SUPPORT.md, feedback_support.py,
  feedback_support_certificate.py, test_feedback_support.py. The bounded
  FLOAT probe is separate and cannot accept a certificate. Existing proof
  helpers remain unchanged. No data files, caches, native benchmarks,
  production/paper changes, commits, pushes, or promotion. Local storage
  exceeds20GiB; jobs are serial and time-bounded, with small bounded arrays.

Next allow nonzero intervening outputs while retaining actual-support
route moments and fresh-output constraints. Replace the single unchanged
polynomial root count by the actual two-state live-to-live rank interface;
do not simply reuse the63-hole cap. Multiple episodes remain required too.
The near-half/40-bit whole-code goal is still open. This turn makes proof
progress without changing the construction or claiming a new runtime.

### 2026-10-22 — Nonzero intervening outputs and two-state constraints

Kept the construction and production code fixed. Extended the actual-mask
route count to the full-length live class with at most32 nonzero even
outputs among y_2,...,y_(n-2). Their locations, supports, and amplitudes are
arbitrary. Both 80/120-digit directed replays cover all518113 hole totals
in74 intervals; every term is <=2^-80 and the class union is <2^-61 at
>49.1% distance. Weakest box H=515584..515710 has log2 upper -81.6787.
This remains a restricted trajectory class, not a whole-code certificate.

Derived the fresh two-state rank on preceding odd and intervening even
output variables, retaining incoming-state dependence. The certified
extension uses its simpler triangular part. Tested the stronger joint
rank on2048 three-pair GF16 trajectories; it improves1125 local pairs.
Also proved that each64-row coset of A_1 B_W is GRS for any |W|<=64,
even with cross-half input support. All988 corresponding small-field
square minors pass. No generic full-rank claim across both cosets is made.

The bounded weighted-pair FLOAT diagnostic preserves the stronger rank
and its actual-map refinement. Selected d=64 and larger points still
fail; the map refinement does not improve the sampled optimized scores.
The d=34 directed replay fails the chosen -80-per-term rule, not a proved
40-bit impossibility. Kept the conservative d<=32 result.

Identified an explicit relaxation artifact: at maximal holes per state,
the coarse pair envelope permits (u,v)=(1,128). Repeating it violates the
already proved150/256 feedback bound across the OTHER pair boundaries.
Thus the next test should carry boundary hole/update profiles and impose
the sparse-update constraint on both parities, before local maximization.
Do not redesign the code to accommodate this impossible relaxed profile.

Added coded_state/CHANGING_PAIRS.md, changing_pair_probe.py,
changing_pair_certificate.py, test_changing_pairs.py, joint_pair_probe.py;
indexed the note in README/DESIGN. Tests include mask coefficients,
dimension bounds, exact extremizers, d=0 regression, and timeout rejection.
No raw data, benchmarks, caches, commits, pushes, promotion, or paper edits.
Bounded serial runs use small arrays; local free space exceeds20GiB.
Unrestricted changes, later activations, and multiple episodes remain open.

### 2026-10-22 — Random-feedback bridge finds a real compact-state obstruction

The user asked for incremental construction/proof development rather than
remaining attached to the structured inner. Preserved all prior work and
replaced feedback by independent uniform128-by64 matrices as a proof model.
Its inverse law has full randomness on every live batch. This immediately
closes all low-weight outputs whose zero-incoming-state batches touch at
most128 regions, including unrestricted full-length live episodes, with
failure <=2^-318848 at49.1%. It is still a restricted class.

The decisive next check was a reset-supported dimension argument. Choose
128 intervals, each joining a suffix/prefix in a pair of adjacent regions.
Select32 batches in108 regions and31 in148, for1029632 symbols, below the
1029701 bad-weight threshold. The shortened outer has exact expected
dimension17728.325879..., while state64 permits at most8192 reset equations.
Whenever the shortened dimension exceeds that cap, a nonzero message has
output entirely inside the intervals, for EVERY inner coefficient choice.
The exact expectation plus a pointwise dimension cap gives setup-failure
lower bounds >1.88% for state64 and >0.269% for state128. States<=138 are
excluded at40-bit margin by the same exact check. No block-independence
assumption, simulation, or weak numerical upper bound is used.

This applies to the structured candidate too. Its partial proofs remain
valid, but the ensemble cannot meet the desired whole-code guarantee.
It does not say every setup is bad or certify state139. Do not continue
the former cross-pair hill climb toward an impossible ensemble-level target.

The next simple ideal bridge uses batch=state=256, Q'=Q+Y, and uniform
feedback; same outer and regional route. This avoids the fixed-B reachable
space cap and removes nonzero quiet batches. It is not certified or timed.
Bounded FLOAT forward/inverse screens retain partial-region geometry;
some selected cases pass and others remain positive. Next close or reject
this more-memory bridge before thinning/structuring its randomness.

Added coded_state/RANDOM_FEEDBACK_BRIDGE.md, random_feedback_bridge.py,
compact_reset_obstruction.py, uniform_feedback_screen.py, and
test_random_feedback_bridge.py. Exact tests pass:48 uniform feedback laws,
4096 toy outputs/672 count comparisons,1296 route restrictions,128 explicit
reset-witness fixtures, and exact K20 threshold decisions. CWC kept the
new lower-bound obstruction separate from failed upper-bound diagnostics.
No production/paper changes, native benchmarks, generated data, commits,
pushes, or promotions. Local free space remains above20GiB; runs are serial
and bounded. The north star remains a fast near-half code with one route.

### 2026-10-22 — Test state256 and isolate the outer/state/reset tradeoff

The user requested structural checks first, then proof iteration and a
comparison with the already-proved uniform inner. The256-state running
sum remains open, not rejected. The new ideal MDS-reset controls below
are complete proofs, not implementations or timing claims.

- Extended the exact expectation-based reset obstruction screen to every
  interval count for batch256. Even allowing a free terminal reset, its
  maximum expected shortened dimension per reset is184.223582 at c26,
  below256. This does not control route tails or setup-adaptive interval
  choices; it only prevents extrapolating the old64/128 exclusion.
- Derived an exact parity-check output enumerator for the running sum,
  including sparse two-output cancellations. Its improved Chernoff bound
  still leaves positive complete-class upper bounds. No bad codeword for
  the256-state ensemble was found or asserted.
- Compared with ideal globally distinct Vandermonde reset columns H_t,
  while retaining independent uniform full-field feedback. The explicit
  update Q'=Q+H_tY has minimum closed-segment weight b+1 and rank b on
  every valid support. This meets the frozen persistent-MDS counting
  interfaces without extending its old b<=255 multiplicative schedule.
- Closed RS[256,128], batch=state512, at >49.1%/40 bits with the original
  route. Separate80/120-digit replays cover4194712 pairs in535 rectangles.
  Every contribution is <=2^-63; the conservative union has <2^23 cells.
  Weakest directed log2 bound is below-63.33389521022477.
- Tested the different-component hypothesis: RS[512,256], batch=state256,
  using GF(2^16) outer coefficients and512 regions of4096 blocks. This
  also closes >49.1%/40 bits. Both80/120-digit replays cover16773302 pairs
  in517 rectangles. Every contribution is <=2^-65, with <2^25 cells in
  the conservative union. Weakest directed bound is below-75.60640837817.
  This is a different outer/route profile, not a proof for the running sum.
- Exact independent adjacency checks and removed/duplicated/moved-boundary
  controls pass. Added40 parity-enumerator checks,1883 reset-column rank
  fixtures,65536 MDS output assignments,15 full-feedback laws,64 inverse
  streams,16 small partitions and48 corruption controls. The older bridge
  regression passes. Immediate timeout rejects without claiming a bound.

Root cause: uniform feedback on a nonzero state is not global uniform
encoding. The1024-output window forces long zero tails; the running sum
permits sparse cancellation. The outer and route also determine how much
shortened source dimension can hide in fragmented intervals. Stronger
reset algebra plus either more state or a larger outer gives a certified
intermediate design. Original state256/RS256 remains a proof question.

Saved IDEAL_RESET_LADDER.md, uniform256_probe.py, ideal_reset_certificate.py,
and test_uniform256.py; updated research indexes. Frozen proof helpers,
production code, paper, old proofs, and old timings remain untouched.
Runs were serial and source-only. Some proposal/replay caps expired and
were explicitly rejected. The larger-outer120-digit complete replay
used separate120-second phase caps and finished in113.15 seconds total;
the corresponding80-digit run took65.45 seconds. These are verifier times.
Only bounded coefficient/cover arrays are retained in memory; no raw data,
generated reports, caches, benchmarks, commits, pushes, or promotions.

Next: retain both certified controls and cheapen the reset maps first,
holding the feedback distribution fixed. Compare the same weakened reset
interface under both outer/state profiles before reducing feedback
randomness too. Preserve the original256-state running-sum path; its
remaining losses are sparse cancellation and separately maximizing route
moments before counting all interval placements. Do not present either
failed bound as a structural rejection or infer encoding time from size.

### 2026-10-22 — Cheaper exact reset maps before weakening their proof

User requested cheaper resets, unchanged uniform feedback, and comparison
of both certified profiles. Retained the exact MDS reset interface instead
of accepting a proof loss and retuning a bound.

- Extended additive-coset coordinates to state512/original RS256 outer and
  state256/RS512 outer. Use degree-graded novel moments and Gray-code cosets.
  The recurrence is Y=X+UZ, Z'=L(delta)Z+C(0)Y, with independent full uniform
  U over GF(2^128). C(0) is a fixed transpose additive FFT; L is an h-layer
  triangular butterfly. The conceptual syndrome still has globally distinct
  MDS columns. Both existing >49.1%/40-bit certificates transfer exactly.
- Counted all reset products. Gray coordinates use16781312 and14688256
  nontrivial subfield-product sites, respectively. A direct shifted FFT has
  fewer sites but about4.23 million GF32 products. Gray has only16128/8064
  GF32 sites, occurring in32 batches per profile. No native timing claimed.
- Iterated with a deterministic diagonal state normalization. Half the
  Gray translations become XOR-only. Including the added reset scaling,
  total products become14147584/12566528. However GF65536 sites increase
  by524288 in each profile. Partial normalization avoids that increase and
  is retained as a native comparison candidate. Fewer scalar products do
  not establish less CPU time or an end-to-end winner.
- Screened deletion of the fixed transform: the highest-index state unit
  vector is invariant under every L(delta), producing a weight-two closed
  output segment. This invalidates the retained proof interface, not a
  claim of a bad outer-compatible codeword. More generally identity
  injection requires an MDS branch map; a sparse transition cannot supply
  the same closed-segment minimum weight. Do not tune a bound around a
  silently weakened reset or hide injection work in a state basis change.

Added CHEAP_RESET_COORDINATES.md and cheap_reset_maps.py, and indexed them.
Exact checks pass14 transform/translation fixtures (including256/512),
58 state trajectories,1883 small-field ranks,32 coupled encoder/inverses,
32 shortcut negative controls, and4 normalized-coordinate identities.
The proof transfer is algebraic, using the prior complete directed receipts;
these toy and actual-dimension algebra checks are not new numerical certs.
Runs are serial, source-only, and bounded at45 seconds with small vectors.
No native benchmark, raw output, cache, production/paper change, commit,
push, or promotion. Uniform feedback remains a proof control, not a fast
implementation claim. CWC kept local obstructions distinct from code failure.

Next: measure direct FFT, Gray, and normalized Gray reset kernels on both
profiles, then include their respective outer and route. Choose by native
cost before weakening reset algebra or reducing feedback randomness.

### 2026-10-22 — Time both resets; exploit activation holes for small groups

Completed the user's timing pass, then investigated a weaker reset while
retaining uniform full-field feedback and both outer/route profiles.

- Added native/reset_native.cpp using existing GFNI arithmetic, batched
  fixed-size butterflies, aligned scratch, and no hot allocations. Exact
  scalar field, polynomial-column, coordinate-stream, and grouped-stream
  checks pass. Optimized and UBSan builds pass on Peach.
- Serial Ryzen7950X/GCC15.2/CPU2 measurements process N=2^21 128-bit
  elements, excluding setup, outer, route, and feedback. Gray reset alone
  takes about6.0ms for state256 and6.5ms for state512. Direct shifted FFT
  is about11.5ms with large prepared-table traffic. Neither normalization
  variant improves the measured Gray timing.31-sample replays support
  the conclusion; do not quote these as full-encoder latencies.
- Tested independent reset groups while preserving total state b and full
  uniform feedback across groups. Groups32 take about2.80ms; groups8 about
  2.0ms, on either profile. These timings have no new distance certificate.
- Derived the actual weaker interface: a closed episode with a active
  groups has rank da and minimum weight(d+1)a. Its first input batch is
  supported only in those groups because state zero implies X=Y. For c
  episodes, v batches, and A total active closed groups, shortened outer
  support has size b(v-c+1)+dA. Its partial regions are still <=2c.
  Inactive-group output zeros also reduce capacity to b(v-2(c-1))+2dA.
  Combined these facts with the existing forward/inverse route moments,
  charging active-group choices. No independence shortcut is used.
- Coarse rank-floor covers fail and were not tuned indefinitely. The new
  activation-aware bound resolves selected classes, e.g larger-outer d32
  at v4096,c64,A252 has FLOAT log2 upper<-249.29. A broader selected screen
  still fails: worst+851.80 over105 larger-outer classes and+1443.94 over
  100 original-outer classes. These are failed bounds, not codewords or
  directed certificates. Intermediate activity is the remaining gap.
- Exact Python checks pass4096 grouped ranks,24 pattern counts,755 support
  geometries; a d2 weight-three reset is retained as a negative control
  against claiming full-state MDS. Added striped_reset_probe.py,
  test_striped_reset.py, and RESET_TIMING_AND_GROUPS.md; indexed results.

Runs were serial and bounded, with >=20GiB free. Builds stayed remote in
the existing scratch directory. The old WinSCP executable could not read
the configured key; source transfer used SCP over the existing SHA-pinned
Plink connection, without changing keys or relaxing host verification.
No generated data, local binaries, library/paper changes, commits, pushes,
or promotion. Prior certified controls and frozen proof helpers remain intact.

Next: inspect a failing intermediate-activity class using its actual
episode lengths, activation holes, and reset/output-zero constraints.
Test structural feasibility before refining joint route/output counting.
If that stalls, add cheap inter-group mixing rather than demanding full
MDS again or shrinking total memory below the earlier obstruction.

### 2026-10-22 — Joint group32 count; screen intergroup countermeasures

User requested continued group32 closure with uniform feedback, then cheap
intergroup mixing if progress stalled. Kept both outer/route profiles and
the full uniform feedback distribution fixed throughout.

- Derived the missing inactive-interior cost H=sum_i(b-da_i)(l_i-2).
  Forward fixed-input probability is at most q^[-b(c-1)-H]; inverse output
  capacity is E-H with E=b(v-2(c-1))+2dA. Terminal inactive coordinates
  are not charged twice. A positive generating function sums active-group
  masks and excursion lengths jointly. Entropy tangents let the inverse
  argument retain the same H rather than maximizing capacity separately.
- On the retained selected screens, the larger-outer group32 worst bound
  improves from+851.80 to+27.01, and the original-outer worst from+1443.94
  to+772.22. Finer markers give+5.81 for the selected larger-outer critical
  class. These are failed FLOAT upper bounds, not distances or failure
  lower bounds. The full group32 certificate remains open.
- Derived a routing variance cap using actual short-excursion geometry.
  Its Bennett envelope did not improve the tested critical classes. Also
  proved a zero-coset weight-enumerator envelope that sums output group
  choices before fractional averaging. Its length-uniform relaxation is
  still too loose (+731.09 at larger-outer v4224,c128); retained separately.
- Tested sparse intergroup companions over GF256. Found irreducible
  x^8+x^4+x^3+x+2 and x^16+x^4+x^3+x+15, with no cyclic-vector return
  through8192 steps and favorable consecutive-batch ranks. However a
  two-batch output can follow an unmixed shift edge. An exact weight-(d+1)
  reset of rank d survives. Rejected this as a short-excursion repair;
  neither a long period nor favorable one-group examples is sufficient.
- The direct fallback pairs32-symbol groups into64-symbol MDS resets,
  preserving total state and uniform feedback. Required transforms and
  translations must actually change; this is not a free XOR relabeling.
  The new joint screen passes84 larger-outer classes (worst<-521.84) and100
  original-outer classes (worst<-185.54). Additional intermediate-activity
  checks pass too. No full certificate or new native timing is claimed.
- Added a three-dimensional FLOAT coverage proposal. Bounded45/150/180s
  runs expired with unchecked frontier regions; e.g150s accepted2234
  rectangles after4476 checks. Feasible-polytope linear maxima and bounded
  proposal caches were tested, but no complete cover or directed replay
  was obtained. The open work is not represented by a fabricated margin.
- Tests pass12 exact weighted-coefficient comparisons with FLOAT envelopes,
  470 activation/variance geometries,1212 tangent comparisons,1820 exact
  paired-reset ranks,15 exact uniform-feedback laws,200 exact polytope
  maxima, and64 exact coset comparisons. Sparse-mixer tests include the
  negative control, not just favorable rank fixtures.

Added GROUPED_RESET_JOINT.md and bounded source helpers; updated research
indexes. CWC kept local algebra, numerical proposals, failed bounds, and
whole-code claims separate. Prior proofs and frozen helpers are unchanged.
No production/paper edits, native benchmarks, raw dumps, caches on disk,
commits, pushes, or promotion. All numerical runs were serial and bounded;
free local space remains above20GiB.

Next: complete the all-class argument for the paired64 reset, then obtain
an independent coverage check and directed replay. The current3D proposal
needs fewer overly conservative subdivisions; selected points are not a
substitute. If the proof closes, measure the actual added transform cost
before calling it a performance winner. Keep group32 and its stronger
joint count as the cheaper open candidate. Feedback cost remains separate.
