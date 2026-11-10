import SpinCodes.Turbo.SparseExactCell22
import SpinCodes.Turbo.RS8ReferenceStepLaw22

/-! Audit: `RS8StepValueFactor22` and `RS8ReferenceStepLaw22` first retain
complete values at a fixed selected set, then average the joint ordered
count profile. Replace their zero-inclusive selected values by exact
nonzero values; the existing shared-support hypergeometric law is reused. -/
noncomputable section
namespace Spin.Turbo
open Finset

theorem regionalNonzeroStepSupport22 {R : ℕ}
    (word : Fin (R * 16) → Packet4) (step : Fin R) :
    regionalPacketSupport22 ((regionalStepWordEquiv22 R) word step) =
      regionalStepSelectedSet22 (regionalPacketSupport22 word) step := by
  ext slot
  simp only [regionalPacketSupport22, regionalStepSelectedSet22,
    Finset.mem_filter, Finset.mem_univ, true_and]
  rfl

theorem regionalNonzeroStepSupport22_eq_iff {R : ℕ}
    (T : Finset (Fin (R * 16))) (word : Fin (R * 16) → Packet4) :
    regionalPacketSupport22 word = T ↔ ∀ step : Fin R,
      regionalPacketSupport22 ((regionalStepWordEquiv22 R) word step) =
        regionalStepSelectedSet22 T step := by
  constructor
  · intro h step
    rw [regionalNonzeroStepSupport22, h]
  · intro h
    apply (regionalStepSetsEquiv22 R).injective
    funext step
    exact (regionalNonzeroStepSupport22 word step).symm.trans (h step)

theorem regionalNonzeroGiven22_step_product {R : ℕ}
    (T : Finset (Fin (R * 16))) :
    (regionalNonzeroGiven22 T).map (regionalStepWordEquiv22 R) =
      Spin.piPMF (fun step : Fin R => regionalNonzeroGiven22 (regionalStepSelectedSet22 T step)) := by
  apply FinPMF.ext
  funext words
  rw [FinPMF.map_equiv_apply, regionalNonzeroGiven22_apply, Spin.piPMF_apply]
  simp_rw [regionalNonzeroGiven22_apply]
  have hsupport : regionalPacketSupport22 ((regionalStepWordEquiv22 R).symm words) = T ↔
      ∀ step : Fin R, regionalPacketSupport22 (words step) = regionalStepSelectedSet22 T step := by
    simpa only [Equiv.apply_symm_apply] using
      regionalNonzeroStepSupport22_eq_iff T ((regionalStepWordEquiv22 R).symm words)
  by_cases h : ∀ step : Fin R,
      regionalPacketSupport22 (words step) = regionalStepSelectedSet22 T step
  · rw [if_pos (hsupport.mpr h)]
    simp_rw [if_pos (h _)]
    rw [regionalStepSelectedSet22_card_sum]
    have hpow : (15 : ℝ) ^ (∑ step : Fin R, (regionalStepSelectedSet22 T step).card) =
        ∏ step : Fin R, (15 : ℝ) ^ (regionalStepSelectedSet22 T step).card :=
      (Finset.prod_pow_eq_pow_sum Finset.univ
        (fun step : Fin R => (regionalStepSelectedSet22 T step).card) (15 : ℝ)).symm
    rw [hpow]
    simp only [one_div]
    rw [Finset.prod_inv_distrib]
  · rw [if_neg (fun hglobal => h (hsupport.mp hglobal))]
    obtain ⟨step, hs⟩ := not_forall.mp h
    symm
    apply Finset.prod_eq_zero (Finset.mem_univ step)
    exact if_neg hs

theorem uniformStepSelectedSet22_bind_nonzero (j : Fin 17) :
    (uniformStepSelectedSet22 j).bind regionalNonzeroGiven22 = packetStepLaw22 j := by
  apply FinPMF.ext
  funext word
  rw [FinPMF.bind_p]
  simp_rw [regionalNonzeroGiven22_apply]
  simp only [mul_ite, mul_zero, Finset.sum_ite_eq, Finset.mem_univ, ite_true]
  rw [uniformStepSelectedSet22_apply, packetStepLaw22_apply]
  change (if (packetStepSupport22 word).card = j.val then
      1 / (Nat.choose 16 j.val : ℝ) else 0) *
      (1 / (15 : ℝ) ^ (packetStepSupport22 word).card) = _
  by_cases h : (packetStepSupport22 word).card = j.val
  · simp only [h, ite_true, div_eq_mul_inv, mul_inv_rev, one_mul]
    ring
  · simp only [h, ite_false, zero_mul]

theorem sparseExactRegionLaw22_step_mixture (R : ℕ)
    (T : Finset (Fin (R * 16))) :
    (sparseExactRegionLaw22 R T).map (regionalStepWordEquiv22 R) =
      (regionalStepProfileLaw22 R T).bind (fun counts =>
        Spin.piPMF (fun step : Fin R => packetStepLaw22 (counts step))) := by
  rw [sparseExactRegionLaw22_eq_support_mixture, FinPMF.map_bind]
  simp_rw [regionalNonzeroGiven22_step_product]
  have hsets : (Spin.Structured.Routing.shuffleLaw T).bind
      (fun U => Spin.piPMF (fun step : Fin R => regionalNonzeroGiven22 (regionalStepSelectedSet22 U step))) =
      (regionalSelectedStepSetsLaw22 R T).bind
        (fun sets => Spin.piPMF (fun step : Fin R => regionalNonzeroGiven22 (sets step))) := by
    rw [regionalSelectedStepSetsLaw22, FinPMF.bind_map]
    rfl
  rw [hsets, regionalSelectedStepSetsLaw22_eq_profile_mixture, FinPMF.bind_assoc]
  congr 1
  funext counts
  change (Spin.piPMF (fun step : Fin R => uniformStepSelectedSet22 (counts step))).bind
    (fun sets => Spin.piPMF (fun step : Fin R => regionalNonzeroGiven22 (sets step))) = _
  rw [Spin.piPMF_bind]
  congr 1
  funext step
  exact uniformStepSelectedSet22_bind_nonzero (counts step)

#print axioms uniformStepSelectedSet22_bind_nonzero
#print axioms sparseExactRegionLaw22_step_mixture
end Spin.Turbo
