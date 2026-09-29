# Dense semantic review

> **Historical snapshot — archived 2026-09-28.** Status statements and task recommendations below describe an earlier stage.
> See [STATUS.md](STATUS.md) and [NATIVE_RESULT.md](NATIVE_RESULT.md) for the closed native distance/rate result and verification scope.

Reviewed 2026-09-28 by the route-law worker. The reviewed source chain connects the mixed dense certificates to the actual native code's minimum-distance failure event. I found no reversed inequality, missing counting factor, density mismatch, or boundary gap in this chain. The review is a source-level semantic check, not a substitute for the final aggregate Lean compilation. At the structural-audit snapshot, 209 of 307 Fourier boxes had passed; the remaining boxes were still being generated and checked. The final mixed certificate and native theorem must therefore be reported as complete only after their own successful checks.

## Experiment and tilted probability

`ConcreteOuter.failureProbability` is defined in `ConcreteOuterCountingProbability.lean:9` as the probability of `ConcreteRoutedEncoder.weight e rows ≤ d` under the actual `experimentLaw L b R`. Its route-expectation identity, `failureProbability_eq_route`, uses the proved permutation-law equality. The dense interface does not replace this experiment with an assumed row law.

The three numerical alternatives reach this same probability. `ConcreteFourier.routed_probability` (`ConcreteFourierRouted.lean:49`) uses the actual five-class Fourier transfer; `DenseOccupationFixed.routed_counted_rate` (`DenseOccupationFixedRate.lean:29`) uses the actual occupation transfer; `DenseScalarExact.routed_counted_rate` uses `ConcreteScalar.routed_scalar_probability`. The scalar transfer has prefactor one. The occupation and Fourier transfers retain their positive witness-coordinate ratios. The route domination contributes exactly `(b+1)^Q (L+1)^b`, with `Q = activeRows rows`.

All alternatives divide a moment bound by `z^d`, with `0 < z ≤ 1`. This is the correct direction: output weight at most `d` implies `z^weight ≥ z^d`. The tilted reference input has density `p*y`; it need not equal the actual density. The likelihood cost is `exp(L*b*binKL (density rows) (p*y))`.

`finite_exponent_le` (`DenseOccupationFixedRate.lean:9`) bounds the exact finite expression by the exponential of

```text
alpha*(m*x+c) + KL(alpha,p) + alpha*KL(x,y)
  + log(radius)/128 - (11/100)*log(z).
```

The KL-chain inequality is used as an upper bound on `KL(alpha*x,p*y)`. The identity `128*R=L*b` gives the `1/128` coefficient. Since `log(z) ≤ 0`, the assumption `d ≤ (11/100)*L*b` gives the displayed upper bound on `-d*log(z)`. The sign at this threshold conversion is correct.

## Common certificate interface and coverage

`PointRate` (`DenseOccupationFamilyRateDefs.lean:10`) quantifies over arbitrary positive `L,b`, compatible round count `R`, threshold `d`, wiring equivalence, and concrete rows. Its conclusion bounds `K` times the actual experiment probability whenever `0 ≤ K ≤ exp(L*b*alpha*(m*x+c))`. It explicitly requires `density rows = alpha*x`, interior total input weight, and the exact dimension and threshold relations. No route-identification or one-step transfer assumption survives into a completed `PointRate` certificate.

`CertifiedBox` (`DenseGeometryRate.lean:9`) refers to the original global rectangle by index and returns a support belonging to the actual `Majorant.refined.supports`. Each box proves a uniform exponent bound on its whole rectangle. Separate convexity of `boxExponent` in `alpha` and `x` justifies the four-corner check; the argument uses the original `(alpha,x)` rectangle. It does not presume that a nonlinear coordinate change preserves boxes.

`Rect.Contains` uses non-strict inequalities on all four sides. `DenseGeometry.check_sound` handles each subdivision by `≤` versus `>`, so the dividing boundary is included. `global_cover` (`DenseOccupationGeometryCover.lean:7`) covers the entire closed rectangle

```text
alpha in [1/10000, 1],   x in [13/125, 112/125].
```

`cover_indices` and `denseRates_of_families` (`DenseOccupationGlobalIndices.lean:10,12`) cover all 1023 global indices using the 433 occupation, 283 scalar, and 307 Fourier entries. Repeated boundary coverage is harmless. Missing coverage would fail this explicit finite index theorem.

The prepared `DenseGeometry.certified_denseRates` uses `eta=4/10000000` and the fixed common prefactor `C=24000000000000`. `PointRate.mono_constant` and `CertifiedBox.mono_constant` only increase the nonnegative right-hand side. Thus scalar prefactor 1 and occupation prefactor 7000000 can be raised to the common `C`; each Fourier box separately checks its witness prefactor against `C`. The numerical size of this fixed constant affects the finite bound but does not prevent its normalized logarithm from tending to zero.

## Parameters, shared seed, and counting

`row_parameters` (`DenseOccupationRowParameters.lean:8`) and `selected_message_parameters` derive the parameters from the actual rows:

```text
Q = number of active rows,
H = total row weight,
alpha = Q/L,
x = H/(b*Q),
density rows = H/(L*b) = alpha*x.
```

Selected active rows lie in the closed density window, and inactive rows are zero. Consequently `x` lies in the covered interval. Even at `alpha=1`, the input density is strictly below one because `x ≤ 112/125 < 1`; the interior-density hypotheses are therefore supplied, not omitted at the upper occupation boundary.

`selected_spectrum_support` and `selected_profile_support` (`DenseOccupationOuterSupport.lean:8,27`) apply the proved finite spectrum envelope and the selected-seed `b^2` penalty. A support chosen at mean density `x` is valid at every individual row density because it belongs to the global refined support family. Affine averaging then gives the exponent `b*Q*(m*x+c)`.

The outer seed remains fixed and shared throughout this argument. `profileMessages_card` factors the number of message tuples for a fixed seed; this is a Cartesian-product counting identity, not an independence assertion about fresh outer seeds. `profile_failure_eq` (`DenseOccupationProfileFailure.lean:26`) uses equality of row cardinalities to identify the routed failure probability for every message in that profile.

`DenseRates.profile_bound` (`DenseOccupationFamilyRateProfiles.lean:35`) handles empty profiles explicitly. For a nonempty profile it selects an actual message, derives its parameters, bounds its profile multiplicity, and invokes `PointRate.scaled`. The positive scale is `outerCost(b)^Q`, where `outerCost(b)=b^2*exp(16*log(b+1)+13)`.

`DenseRates.layer_bound` then sums profiles and active positions. Its complete finite prefactor is

```text
choose(L,Q) * (b+1)^(2*Q) * outerCost(b)^Q * (L+1)^b * C.
```

The two powers `(b+1)^Q` have distinct sources: profile enumeration and row-route domination. The factor `choose(L,Q)` counts active position sets. The region-route factor is `(L+1)^b`. None of these factors is discarded before the asymptotic estimate.

## Native schedule, selection, and final event

`native_dense_density` (`ConcreteNativeDenseRates.lean:9`) maps every occupation in `Ico (nativeCut m+1) (Lsched m+1)` into the covered alpha interval. The integer-floor cutoff actually gives a strict lower bound before weakening it to `1/10000`. The upper endpoint `Q=Lsched m` is included.

`dense_layers_of_rates` (`ConcreteNativeDenseRates.lean:24`) supplies `b=nativeBlocks m*24=bsched m`, `128*rounds m=Lsched m*bsched m`, and the actual `tupleWiring m`. Its threshold premise follows from `threshold_le`; the final event uses `threshold_floor` to state exactly `floor((0.11:Real)*Nsched m)`. It does not replace that integer threshold by a rounded asymptotic surrogate.

`EZ_le_tuple_layer_bound` (`ConcreteNativeDenseIntegration.lean:12`) averages the pointwise selected-seed bound under the actual conditional outer law. It therefore needs no additional reciprocal selection probability. Eventual positivity of `nativeGood` supplies equality with the total fallback selection definition. The eventual-regime assembly separately pays the actual unconditioned selection-failure probability, proved to tend to zero by the outer-tail closure.

`densePrefactor_le_exp` and `denseRemainder_tendsto` (`DenseNativeRemainder.lean:24,77`) retain all finite costs, with

```text
denseRemainder(C,m) = log(C)/N
  + (log(2)+20*log(b+1)+13)/b + log(L+1)/L.
```

Every term tends to zero along the native schedule. `native_dense_sum` pays at most `L` occupation layers. The eventual dense bound is therefore `L*exp(-eta*N+o(N))`, which tends to zero for the positive certified `eta`. `minimum_distance_of_fixed_and_dense_rates` combines this with the checked fixed and growing-sparse regimes and returns the original seed-law probability of the actual realized code's minimum-distance failure event.

## Incremental Fourier identity audit and verification at this checkpoint

`scripts/route_fourier_index_audit.py` independently reads the source Fourier list, the four replay manifests, emitted module text, and the prepared aggregate. At the review snapshot it passed for 209 checked boxes and 213 emitted boxes, with no duplicate ownership. It checked the source SHA-256 hashes of checked box modules; each record contains four successful axiom audits. It also checked every emitted filename/namespace, its data and witness imports, its global geometry index, and its `certified_uniform` index and constants.

The aggregate imports and proof references are exactly `B000` through `B306` in order. The global Fourier index list agrees with all 307 source entries, and the deduplicated witness list has 105 entries. All 105 witness replays have completed. Missing future box files during an import-graph dry run are expected until their owning lane emits them; for example, B094 was subsequently emitted and passed.

Completion update: all four lanes subsequently finished, and `scripts/map_data/dense_fourier_complete_verification.json` reports PASS for all 105 witnesses and all 307 boxes. The final verifier checked all 824 current source/object pairs, all 1543 standard-only axiom audits, exact lane ownership, and the complete structural audit. These were individual kernel checks with previously checked dependencies reused; the record does not claim a fresh full dependency-closure replay.

At this checkpoint, the remaining action was the separately owned mixed aggregate plus final native theorem/pin checks. These checks are now complete; see [STATUS.md](STATUS.md). No source-level semantic repair was identified by this review.
