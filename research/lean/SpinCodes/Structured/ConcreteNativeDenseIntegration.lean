import SpinCodes.Structured.ConcreteNativeTwoRegimes
import SpinCodes.Structured.ConcreteNativeLinearCodeDistance
import SpinCodes.Structured.DenseNativeSum

noncomputable section
attribute [local instance] Classical.propDecidable
namespace Spin.Structured.ConcreteNativeFamily
open Finset Filter ConcreteOuter ConcreteRoute DenseOccupationFixed

/-- A pointwise bound for each selected shared seed bounds the conditional first moment
without any extra reciprocal selection probability. -/
theorem EZ_le_tuple_layer_bound (m Q : ℕ)
    (hpos : 0 < (nativeSeedLaw m).prob (nativeGood m)) {B : ℝ}
    (h : ∀ seed, nativeGood m seed →
      (∑ x ∈ univ.filter (fun x : NativeMessage m => occupation x = Q),
        failureProbability (tupleWiring m) (rowSupports seed x) (threshold m)) ≤ B) :
    concreteFamily.EZ m Q ≤ B := by
  change (((nativeSeedLaw m).condition (selectedGood m) (selectedGood_pos m)).prod (nativeSetup m).Pin).expect
    (fun ω => ((nativeSetup m).ZQ (threshold m) (occ m) Q ω : ℝ)) ≤ _
  rw [Setup.expect_ZQ_eq_outerLaw]
  calc
    _ ≤ ((nativeSeedLaw m).condition (selectedGood m) (selectedGood_pos m)).expect (fun _ => B) := by
      apply FinPMF.expect_condition_mono
      intro seed hg
      rw [selectedGood_eq m hpos] at hg
      rw [qd_sum_eq_tuple]
      apply le_trans _ (h seed hg)
      apply sum_le_sum_of_subset_of_nonneg
      · intro x hx
        exact mem_filter.mpr ⟨mem_univ _, (mem_filter.mp hx).2⟩
      · intro x _ _
        exact FinPMF.prob_nonneg _ _
    _ = _ := FinPMF.expect_const _ _

/-- Uniform dense tuple-layer bounds supply the final dense asymptotic regime. -/
theorem dense_eventually_of_layer_bound {C η : ℝ} (hC : 0<C)
    (h : ∀ᶠ m in atTop, ∀ Q ∈ Ico (nativeCut m+1) (Lsched m+1),
      ∀ seed, nativeGood m seed →
        (∑ x ∈ univ.filter (fun x : NativeMessage m => occupation x = Q),
          failureProbability (tupleWiring m) (rowSupports seed x) (threshold m)) ≤
        densePrefactor (Lsched m) (bsched m) Q C * Real.exp (-η*(Nsched m : ℝ))) :
    ∀ᶠ m in atTop, ∑ Q ∈ Ico (nativeCut m+1) (Lsched m+1), concreteFamily.EZ m Q ≤
      (Lsched m : ℝ)*Real.exp (-η*(Nsched m : ℝ)+denseRemainder C m*(Nsched m : ℝ)) := by
  filter_upwards [native_good_positive_eventually, h] with m hm hh
  apply native_dense_sum m (nativeCut m) (concreteFamily.EZ m) hC η
  intro Q hQ
  exact EZ_le_tuple_layer_bound m Q hm (hh Q hQ)

/-- Exact paper minimum-distance event; fixed moments and uniform dense layer estimates
are the remaining hypotheses. All construction, selection, sparse, and remainder steps are proved. -/
theorem minimum_distance_of_fixed_and_dense_layers
    (hfixed : ∀ Q ∈ Ico 1 4096, Tendsto (fun m => concreteFamily.EZ m Q) atTop (nhds 0))
    {C η : ℝ} (hC : 0<C) (hη : 0<η)
    (hdense : ∀ᶠ m in atTop, ∀ Q ∈ Ico (nativeCut m+1) (Lsched m+1),
      ∀ seed, nativeGood m seed →
        (∑ x ∈ univ.filter (fun x : NativeMessage m => occupation x = Q),
          failureProbability (tupleWiring m) (rowSupports seed x) (threshold m)) ≤
        densePrefactor (Lsched m) (bsched m) Q C * Real.exp (-η*(Nsched m : ℝ))) :
    Tendsto (fun m => ((nativeSeedLaw m).prod
      (ConcreteRoutedEncoder.experimentLaw (Lsched m) (bsched m) (rounds m))).prob
        (fun ω => minimumDistance m ω.1 ω.2 ≤ ⌊(0.11 : ℝ)*Nsched m⌋₊)) atTop (nhds 0) := by
  have h := distance_of_fixed_dense hfixed η hη (denseRemainder C) (denseRemainder_tendsto C)
    (dense_eventually_of_layer_bound hC hdense)
  apply h.congr
  intro m
  rw [concrete_probBad_minimumDistance, threshold_floor]

#print axioms EZ_le_tuple_layer_bound
#print axioms dense_eventually_of_layer_bound
#print axioms minimum_distance_of_fixed_and_dense_layers

end Spin.Structured.ConcreteNativeFamily
