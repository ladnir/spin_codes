# SpinCodes — Lean formalizations

The public native theorem below verifies the earlier Structured SPIN
construction. It does not verify the revised Scalable SPIN or finite Turbo
SPIN results in the current paper.

The new [uniform-mixer kernel](SpinCodes/UniformMixer/Kernel.lean) checks a
first part of the Scalable SPIN proof: the exact state transition under a
uniform invertible linear map, the full-state formula for an empty gap,
and the one-bit impulse quotient. The [regional algebra](SpinCodes/UniformMixer/Regional.lean)
checks full-state kernel composition, the zero/live lift and projection,
the weighted impulse bound, and both gap comparisons conditional on a
scalar live-weight estimate. The [mean](SpinCodes/UniformMixer/Mean.lean),
[production map](SpinCodes/UniformMixer/Map22.lean), and
[tilt bounds](SpinCodes/UniformMixer/Tilt.lean) prove the exact live-state
mean for the paper's 22-by-64 matrix and a quadratic error bound for its
empty-gap tilt. They also give a full-state empty-gap sandwich for an
arbitrary entering state, with an explicit `1±C/L` factor at tilt `θ/L`
when the gap fits in one region. The
[production impulse](SpinCodes/UniformMixer/Regional22.lean) identifies
the actual one-bit output and feedback maps, proves each feedback column
nonzero, and checks its weighted bounds and zero/live quotient. The
[regional product](SpinCodes/UniformMixer/RegionalProduct.lean) composes
these bounds for any fixed placement with separated active steps. It
retains arbitrary full-state boundaries and an explicit finite-length
factor. The
[matrix bridge](SpinCodes/UniformMixer/RegionalMatrix.lean) now proves
the reduction to the two-state product for every fixed separated
placement, with the full-state boundary lift and projection intact.
The [placement bound](SpinCodes/UniformMixer/RegionalPlacement.lean)
uses the verified uniform shuffle law to bound boundary and close-pair
events and proves that their complement has separated, interior steps.
Its exclusion event is stronger than the one in the paper and still has
probability `O(a²/L)` for fixed occupation `a`. The
[support conversion](SpinCodes/UniformMixer/RegionalSupport.lean) sorts
each good support into local coordinates and positive gaps, then applies
the two-state bound to its marked kernel. The
[stepwise encoder](SpinCodes/UniformMixer/RegionalEncoder.lean) defines
`Y=Aq+x` and feedback `Aᵀx`; its empty and occupied steps equal the
kernels used by that bound. The
[chronological trace](SpinCodes/UniformMixer/RegionalTrace.lean) proves
that a list expanded from marked gaps has exactly the marked-region
kernel and that a good support expands to the correct number of steps.
It verifies the empty-support case and transfers the regional upper bound
to the actual step law. For every good support, it proves the expanded
list equals the encoder's chronological input list, so the regional bound
now applies directly to the encoder. The
[regional average](SpinCodes/UniformMixer/RegionalAverage.lean) proves that
the actual shuffled kernel is nonnegative and substochastic, charges the
discarded placements by their `O(a²/L)` probability, and averages the
good-placement two-state envelope under the verified uniform shuffle law.
The [site-sum bridge](SpinCodes/UniformMixer/RegionalSiteSum.lean) rewrites
that average as an exact fixed-cardinality support sum, ready for the
discrete-to-continuum comparison. The
[grid-limit theorem](SpinCodes/UniformMixer/RegionalGridLimit.lean) proves
convergence of that uniform support sum for any continuous function of
the normalized sorted bit positions. The
[continuum product](SpinCodes/UniformMixer/RegionalContinuum.lean) defines
the paper's zero/live gap-and-impulse integrand, proves its continuity,
and verifies convergence of its ideal support sum. The
[product stability bound](SpinCodes/UniformMixer/RegionalContinuumError.lean)
shows that its two-state rows change by at most `γ` times the `ℓ₁` change
in gap lengths. The [rounding bound](SpinCodes/UniformMixer/RegionalRounding.lean)
shows that every good support with `a` active bits has total normalized
gap error at most `3(a+1)/R`, where the region has `R` encoder steps.
The [recursive product comparison](SpinCodes/UniformMixer/RegionalProductError.lean)
applies this error to the actual rounded and ideal two-state products.
The [product bridge](SpinCodes/UniformMixer/RegionalProductBridge.lean)
identifies the rounded product exactly with the marked regional envelope
and gives a uniform good-support row error of `3γ(a+1)/R`.
The [site-product bridge](SpinCodes/UniformMixer/RegionalSiteProduct.lean)
identifies the ideal recursive product with the continuum integrand at
the normalized sorted bit positions. The
[boundary error](SpinCodes/UniformMixer/RegionalBoundaryError.lean)
transfers the good-support row bound to every full-state endpoint entry.
The [production envelope comparison](SpinCodes/UniformMixer/RegionalEnvelopeComparison.lean)
inserts this bound into the actual finite regional envelope, with its
proved multiplicative error factor. The
[averaged continuum comparison](SpinCodes/UniformMixer/RegionalAverageContinuum.lean)
then sums over the exact uniform shuffle law and gives a finite bound
for every fixed active-support size. The
[regional limit](SpinCodes/UniformMixer/RegionalLimit.lean) proves the
uniform asymptotic upper bound for the actual shuffled setup law.
The [Scalable setup](SpinCodes/UniformMixer/ScalableSetup.lean) packages
the already verified shared Golay–BAA outer and one-bit route with
independent uniform GL(22,2) state maps. It proves the outer and route
marginals, the selected-outer failure limit, and an expectation identity
for routing outer words under the joint law. The
[one-step law](SpinCodes/UniformMixer/ScalableStepLaw.lean) proves that a
sampled state map has exactly the transition probability used by the
regional kernel. The [realized encoder moment](SpinCodes/UniformMixer/ScalableEncoderMoment.lean)
defines the sampled 64-bit-step encoder and proves its moment recurrence.
The [encoder–kernel identity](SpinCodes/UniformMixer/ScalableEncoderKernel.lean)
identifies that moment with the chronological regional kernel for every
fixed routed support. The
[full realization](SpinCodes/UniformMixer/ScalableRealization.lean) defines
the routed encoder weight under the joint setup, with state carried across
region boundaries. The [flattened-map law](SpinCodes/UniformMixer/ScalableMapLaw.lean)
shows that the region-by-step map array becomes independent uniform maps
in encoder order. The [full moment identity](SpinCodes/UniformMixer/ScalableFullMoment.lean)
integrates the realized output weight over the actual outer, route, and
inner setup and identifies its moment with the chronological kernel.
The [region composition](SpinCodes/UniformMixer/ScalableRegionComposition.lean)
identifies that kernel with the ordered product of regional full-state
kernels, preserving the state across boundaries. The
[independent shuffle product](SpinCodes/UniformMixer/ScalableRouteAverage.lean)
and [casted shuffle law](SpinCodes/UniformMixer/ScalableCastShuffle.lean)
show that the actual region shuffles average to the previously defined
regional kernels. The [routed moment](SpinCodes/UniformMixer/ScalableRoutedMoment.lean)
expresses the full setup moment as an expectation over row shuffles of
their ordered product. The
[route domination bridge](SpinCodes/UniformMixer/ScalableRouteDomination.lean)
applies the checked uniform-permutation comparison to that exact
nonnegative moment. Under the selected shared Golay–BAA outer event,
every nonzero message has the weight range needed by this comparison.
The [iid Bernoulli step](SpinCodes/UniformMixer/ScalableBernoulliStep.lean)
defines the exact 64-bit reference input law and proves that its independent
step kernel is a matrix power. The [iid wiring](SpinCodes/UniformMixer/ScalableIidWiring.lean)
shows that the route's iid region bits become these independent step inputs.
The [selected-outer moment bound](SpinCodes/UniformMixer/ScalableFirstMomentBound.lean)
then compares the actual route-and-inner moment of every nonzero message to
the row sum of that matrix power. The quantitative matrix estimate,
outer-word summation, numerical first-moment gates, and finite Turbo
certificates remain open.
The [weighted-row formulas](SpinCodes/UniformMixer/ScalableWeightedRows.lean)
identify the zero-state event and bound a live state's individual target
atoms and return probability under the actual production maps. The
[three-coordinate envelope](SpinCodes/UniformMixer/ScalableEnvelope.lean)
proves one-step domination for any nonnegative kernel satisfying the
listed zero, atom, and averaged-row bounds. The
[row instantiation](SpinCodes/UniformMixer/ScalableEnvelopeRows.lean)
verifies those conditions for the actual Bernoulli step kernel, including
both minimum bounds on zero returns. The
[iteration](SpinCodes/UniformMixer/ScalableEnvelopeIterate.lean)
proves the full-state moment bound by the corresponding three-by-three
matrix power. The [transfer first moment](SpinCodes/UniformMixer/ScalableTransferFirstMoment.lean)
applies this bound to every nonzero message under a selected shared outer.
The [positive-vector lemma](SpinCodes/UniformMixer/ScalableTransferCollatz.lean)
turns any checked Collatz witness for that matrix into a geometric bound
on the full-state moment and applies it to an actual selected message.
The [Bernoulli row identity](SpinCodes/UniformMixer/ScalableMomentSpectrum.lean)
expresses each live row exactly in terms of the production expansion
weight and expresses zero return as the actual feedback-kernel weight
polynomial. The production map's spectrum and kernel distance, numerical
witnesses, the outer-word sum, and distance gates remain open.
The [row-mask ANF check](SpinCodes/UniformMixer/Map22Quadratic.lean)
kernel-checks that every production row has Boolean degree at most two.
The [pivot inverse certificate](SpinCodes/UniformMixer/Map22Rank.lean)
proves the production expansion has full rank 22.
The [Boolean-lattice parity lemma](SpinCodes/UniformMixer/Map22SetParity.lean)
and [six-bit subset equivalence](SpinCodes/UniformMixer/Map22SubsetEquiv.lean)
give the [submask parity](SpinCodes/UniformMixer/Map22SubmaskParity.lean)
used in [ANF reconstruction](SpinCodes/UniformMixer/Map22ANFReconstruct.lean).
The [production polynomial theorem](SpinCodes/UniformMixer/Map22Polynomial.lean)
then proves that each actual output coordinate evaluates an ANF of degree
at most two. The [quadratic code-space equivalence](SpinCodes/UniformMixer/Map22QuadraticSpace.lean)
counts the 22 allowed monomials and proves that every such polynomial is
realized by the production expansion. Its image spectrum, minimum live
weight, and the dual feedback-kernel distance remain open.
The [packed spectrum pilot](SpinCodes/UniformMixer/Map22Packed.lean) pins the
old blockwise checker to the production row masks, proves agreement on each
input basis vector, and checks two 256-state blocks at opposite ends of the
22-bit input range. It does not yet connect arbitrary packed inputs to
`coordMap` or certify the complete image spectrum.
The final Scalable and Turbo distance theorems remain open. The older
native theorem below does not verify those revised constructions.

The native Structured SPIN distance and rate result is complete. Under the
actual setup distribution, the probability that minimum distance is at most
`floor (0.11 * N)` tends to zero along the native lengths. Every realization
has exact rate `1/2`.

## Use the result

```lean
import SpinCodes.Native

open Spin.Structured.ConcreteNativeFamily

#check native_minimum_distance_failure_tendsto
#check realizedCode_rate
#check relative_distance_success_tendsto
#check eventually_exists_rate_half_distance_gt_eleven_percent
```

`SpinCodes.Native` is the documented public import. It exposes the original
checked declarations without aliases or replacement definitions.
[SpinCodes/Native.lean](SpinCodes/Native.lean) explains their objects and scope.
[SpinCodes/NativePin.lean](SpinCodes/NativePin.lean) prints all four full
statements and their recursive axiom dependencies.

The law samples one shared Golay-BAA outer constituent, independently of the
actual routing permutations and IMT transvections. The final probability is
unconditioned. Relative distance strictly above `11/100` holds with probability
tending to one, and a positive-probability good setup exists at every
sufficiently large native index.

## Build and verification

Use the pinned toolchain and dependency manifest; run commands from this directory.
With elan and the pinned dependencies installed:

```sh
lake exe cache get
lake build SpinCodes.Native
lake env lean -j1 -M 20000 SpinCodes/NativePin.lean
```

The build command creates missing or stale project dependencies. On a fresh
checkout, this includes substantial numerical certificate checking. It is
not a quick smoke check. The final command imports the built public result
and prints its statements and axiom dependencies.

The new `UniformMixer` modules were checked with the pinned toolchain
and cached dependency objects. Their final declarations report only
`propext`, `Classical.choice`, and `Quot.sound`. A clean replay of the
revised proof is still required after its remaining modules are complete.

Bare `lake build` still selects the historical `SpinCodes` default target.
That target and `scripts/check.sh` check the earlier abstract framework;
they do not by themselves establish the concrete native theorem.

The completed closure checks report only `propext`, `Classical.choice`, and
`Quot.sound` in the final theorem dependencies. Numerical certificates are
checked by Lean's kernel. The closing run reused previously built objects;
it was not a clean replay of every transitive project source.
[FINAL_REPRODUCTION.md](FINAL_REPRODUCTION.md) gives the exact workflows,
resource limits, platform prerequisites, and evidence distinctions.

## Read and maintain the proof

| Document | Purpose |
| --- | --- |
| [NATIVE_RESULT.md](NATIVE_RESULT.md) | Exact construction, theorem declarations, and result scope |
| [PROOF_GUIDE.md](PROOF_GUIDE.md) | Reading order and dependencies from the encoder to the final theorem |
| [STATUS.md](STATUS.md) | Current completion status and recommended follow-up |
| [CLOSURE_AUDIT.md](CLOSURE_AUDIT.md) | Requirement-by-requirement closure evidence |
| [FINAL_EVIDENCE_REVIEW.md](FINAL_EVIDENCE_REVIEW.md) | Source/object provenance and archival coverage limits |
| [FINAL_REPRODUCTION.md](FINAL_REPRODUCTION.md) | Reproduction commands and full source replay |
| [CONTRIBUTING.md](CONTRIBUTING.md) | Editing proofs, generated certificates, and verification records |
| [scripts/README.md](scripts/README.md) | Verification drivers, generators, and report layout |

`ROADMAP.md`, `TASKS.md`, checkpoint notes, and `LOOP_LOG.md` retain the
history of the proof's development. Their earlier open obligations are
superseded by the current status and final checked declarations.

The verified result covers native-length distance and exact rate. It does
not include encoder/transposed-encoder operation counts or the complete
arbitrary requested-length wrapper. Some fixed-occupation intermediate
bounds are more conservative than the paper's; see
[FIXED_NATIVE_SEMANTIC_REVIEW.md](FIXED_NATIVE_SEMANTIC_REVIEW.md).
