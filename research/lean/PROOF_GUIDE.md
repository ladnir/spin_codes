# Reading the native distance proof

Start with the public import:

```lean
import SpinCodes.Native
```

The principal namespace is `Spin.Structured.ConcreteNativeFamily`. The checked
result concerns the actual finite encoder and its sampled setup. Its final
statements have no remaining selection, numerical, or asymptotic premises.

## Construction and probability law

For native index `m`, write `L = Lsched m`, `b = bsched m`, and `N = Nsched m`.
[Schedule.lean](SpinCodes/Structured/Schedule.lean) defines `L = 128*(m+1)`,
chooses the least positive multiple of 24 satisfying
`b ≥ (39/4)*log₂(L*b)`, and sets `N = L*b`.
Messages have `N/2` binary coordinates; encoded words have `N` coordinates.

The setup has three sources of randomness:

1. [ConcreteOuter.lean](SpinCodes/Structured/ConcreteOuter.lean) samples two
   independent uniform permutations for the Golay direct sum followed by two
   permute–accumulate stages. **The same outer seed is shared by every row.**
   [ConcreteOuterNative.lean](SpinCodes/Structured/ConcreteOuterNative.lean)
   supplies `nativeSeedLaw` and the exact native dimensions.
2. [ConcreteRoutePermutation.lean](SpinCodes/Structured/ConcreteRoutePermutation.lean)
   samples independent uniform permutations within each row and within each
   transposed region. `permutation_law_eq` identifies this experiment with the
   support-distribution law used in the estimates.
3. [ConcreteEncoder.lean](SpinCodes/Structured/ConcreteEncoder.lean) samples one
   independent uniform element of `Spin.transPairs 19` per 128-bit round.
   A pair `(u,v)` satisfies `u ≠ 0` and `⟨u,v⟩ = 0`.

The inner state is a 19-bit vector, initially zero. For input block `X_i`,
the round emits `Y_i = X_i + A q_i` and updates
`q_(i+1) = (I + u_i v_iᵀ) q_i + C X_i`, over `𝔽₂`.
The state persists across regions, and there is no terminal flush.
[ConcreteMaps.lean](SpinCodes/Structured/ConcreteMaps.lean) defines the table
maps; [ConcreteSerialization.lean](SpinCodes/Structured/ConcreteSerialization.lean)
fixes the region-major stream order.

[ConcreteRoutedMoment.lean](SpinCodes/Structured/ConcreteRoutedMoment.lean)
combines the route and transvection randomness in `experimentLaw`.
[ConcreteNativeFamily.lean](SpinCodes/Structured/ConcreteNativeFamily.lean)
uses the independent product of `nativeSeedLaw m` and this inner law in
`nativeSetup`. That original product law is the probability space in the
final theorem; conditioning appears only inside the proof.

## Statements to inspect first

The declarations below are in `Spin.Structured.ConcreteNativeFamily`.

| Declaration | Checked conclusion |
| --- | --- |
| `native_minimum_distance_failure_tendsto` | `Pr[minimumDistance ≤ floor(0.11*N)] → 0` as `m → ∞`. |
| `realizedCode_rate` | Every realization has binary dimension divided by length exactly `1/2`. |
| `relative_distance_success_tendsto` | The probability of strict relative minimum distance above `11/100` tends to one. |
| `eventually_exists_rate_half_distance_gt_eleven_percent` | Every sufficiently large native index admits such a realization with positive setup probability. |

Read [ConcreteNativeTheorem.lean](SpinCodes/Structured/ConcreteNativeTheorem.lean)
for the short final assembly and
[ConcreteNativeDistanceConsequencesFinal.lean](SpinCodes/Structured/ConcreteNativeDistanceConsequencesFinal.lean)
for the last two statements.
[ConcreteNativeLinearCode.lean](SpinCodes/Structured/ConcreteNativeLinearCode.lean)
and [ConcreteNativeLinearCodeDistance.lean](SpinCodes/Structured/ConcreteNativeLinearCodeDistance.lean)
connect the encoder to the realized binary linear code and its minimum distance.

## How the three ranges close

Let `Q` be the number of active outer rows. The first-moment argument counts
nonzero messages whose emitted weight is at most the distance threshold.
[Framework.lean](SpinCodes/Framework.lean), particularly
`Spin.Setup.prob_bad_le_cond_sum`, bounds failure probability by the outer
selection failure plus the sum of conditional first moments over `Q`.

`concreteFamily.EZ m Q` is that conditional first moment.
[ConcreteNativeTotal.lean](SpinCodes/Structured/ConcreteNativeTotal.lean)
makes conditioning total using `selectedGood`; its fallback does not change
the encoder or the unconditioned failure event.
[ConcreteOuterTailClosure.lean](SpinCodes/Structured/ConcreteOuterTailClosure.lean)
proves that selection fails with probability tending to zero.

For sufficiently large `m`, the remaining sum splits into these ranges:

| Range | Main reading path |
| --- | --- |
| Fixed: `1 ≤ Q < 4096` | [ConcreteNativeFixedAll.lean](SpinCodes/Structured/ConcreteNativeFixedAll.lean) combines separate `Q=1`, `Q=2`, and `Q≥3` arguments. The finite experiment reaches its continuum kernel through [ConcretePlacementLimitActual.lean](SpinCodes/Structured/ConcretePlacementLimitActual.lean); [ConcreteFixedInsertionNorm.lean](SpinCodes/Structured/ConcreteFixedInsertionNorm.lean) supplies the matrix-norm interface. Only a finite sum of fixed-`Q` limits is used. |
| Growing sparse: `4096 ≤ Q ≤ floor(L/10000)` | [ConcreteNativeSparseSelected.lean](SpinCodes/Structured/ConcreteNativeSparseSelected.lean) combines the selected outer spectrum with the sparse IMT bound. [ConcreteOuterTailClosure.lean](SpinCodes/Structured/ConcreteOuterTailClosure.lean) exposes the completed `native_sparse_eventually` result. |
| Dense: `floor(L/10000) < Q ≤ L` | [DenseOccupationMixedCertified.lean](SpinCodes/Structured/DenseOccupationMixedCertified.lean) combines scalar, occupation, and Fourier certificates. [ConcreteNativeDenseRates.lean](SpinCodes/Structured/ConcreteNativeDenseRates.lean) transports those rates to actual message counts and the distance event. |

[EventualRegimes.lean](SpinCodes/Structured/EventualRegimes.lean) supplies the
eventual range assembly. Earlier bridge files intentionally retain explicit
premises; the final theorem supplies them with the checked results above.

## Source layout and numerical checks

| Location | Role |
| --- | --- |
| `SpinCodes/Prob.lean`, `FiniteLaw.lean`, `Framework.lean`, `Distance.lean`, `Selection.lean` | Finite probability, first moments, distance, and selection. |
| `SpinCodes/Structured/` | Construction and analytic bridges. Most `Concrete*` modules connect an actual finite experiment to a bound. Table-data modules and generated numerical families also live here. |
| `SpinCodes/Numeric/` | Integer interval evaluators and proofs that their successful checks imply real inequalities. |
| `SpinCodes/Cover/`, `SpinCodes/Majorant/` | Outer exponent covers, support data, and their assemblies. These directories mix generated certificates with mathematical interfaces. |
| Data subdirectories of `SpinCodes/Structured/` | `MapSpectrumData/`, `FiberNumericsData/`, `LowCancellationData/`, `SparseData/`, and `SparseBridge/` hold generated tables, finite checks, and polynomial certificates. |
| Dense families in `SpinCodes/Structured/` | `DenseOccupationFixed*`, `DenseScalarExact*`, and `DenseFourierExact*` mix generated data/check wrappers with handwritten soundness and transfer modules. |
| `scripts/` and its `map_data/`, `majorant_data/`, `sparse_data/` | Generators, check orchestration, witnesses, logs, and verification records. |

Read an assembly theorem and its soundness interface before opening thousands
of generated entries. Inspect file headers and generator scripts before editing
these families; regenerate their generated files through the corresponding
scripts. Python and remote execution arrange checks; Lean checks the proof
terms and numerical reductions. The final axiom audits contain only
`propext`, `Classical.choice`, and `Quot.sound`.

[FINAL_REPRODUCTION.md](FINAL_REPRODUCTION.md) gives the check commands.
[FINAL_EVIDENCE_REVIEW.md](FINAL_EVIDENCE_REVIEW.md) distinguishes individual
source/object compiler records from dependency snapshots and historical
source-only records. Final assembly reused checked imports; it was not a
fresh replay of every transitive project source.

## Scope relative to the paper

The final distance and rate statements are checked. The dense assembly uses
`DenseRates (4/10000000) 24000000000000`: a uniform exponent margin of
`4×10⁻⁷` and a sufficient finite prefactor. The
[paper](../paper/structured_proof.tex) reports the sharper numerical margin
`4.10335×10⁻⁷`; that stronger value is not the margin asserted by this Lean
assembly. The checked margin suffices for the same asymptotic 11% distance
conclusion. Encoder bit-operation bounds remain outside the checked result.
