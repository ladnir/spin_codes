import SpinCodes.Turbo.RS8StepGeometry22

/-!
# Count profiles and uniform selected sets within physical steps

This is the conditional target for the regional placement law. Each
step has a uniform selected subset at its specified count. The remaining
gate is to prove that the actual joint law disintegrates into this
product law after observing the entire count profile.
-/

noncomputable section

namespace Spin.Turbo

open Finset

def regionalStepCounts22 {R : ℕ}
    (sets : Fin R → Finset (Fin 16)) : Fin R → Fin 17 :=
  fun step =>
    ⟨(sets step).card, by
      have h := Finset.card_le_univ (sets step)
      have h16 : (sets step).card ≤ 16 := by simpa using h
      omega⟩

def uniformStepSelectedSet22 (j : Fin 17) :
    FinPMF (Finset (Fin 16)) :=
  (FinPMF.uniform (SelectedSlotSet22 j)).map Subtype.val

theorem uniformStepSelectedSet22_apply (j : Fin 17)
    (T : Finset (Fin 16)) :
    (uniformStepSelectedSet22 j).p T =
      if T.card = j.val then
        1 / (Nat.choose 16 j.val : ℝ) else 0 := by
  classical
  rw [uniformStepSelectedSet22, FinPMF.map_p]
  by_cases h : T.card = j.val
  · rw [if_pos h]
    let tt : SelectedSlotSet22 j := ⟨T, h⟩
    have he (U : SelectedSlotSet22 j) :
        U.val = T ↔ U = tt := by
      constructor
      · intro hU
        exact Subtype.ext hU
      · rintro rfl
        rfl
    simp only [he, Finset.sum_ite_eq', Finset.mem_univ, if_pos]
    change 1 / (Fintype.card (SelectedSlotSet22 j) : ℝ) = _
    rw [selectedSlotSet22_card]
  · rw [if_neg h]
    apply Finset.sum_eq_zero
    intro U _
    have hU : U.val ≠ T := by
      intro hEq
      apply h
      rw [← hEq]
      exact U.property
    simp [hU]

def uniformStepSelectedSets22 {R : ℕ}
    (counts : Fin R → Fin 17) :
    FinPMF (Fin R → Finset (Fin 16)) :=
  Spin.piPMF (fun step => uniformStepSelectedSet22 (counts step))

theorem uniformStepSelectedSets22_apply {R : ℕ}
    (counts : Fin R → Fin 17)
    (sets : Fin R → Finset (Fin 16)) :
    (uniformStepSelectedSets22 counts).p sets =
      ∏ step : Fin R,
        if (sets step).card = (counts step).val then
          1 / (Nat.choose 16 (counts step).val : ℝ) else 0 := by
  simp only [uniformStepSelectedSets22, Spin.piPMF_apply,
    uniformStepSelectedSet22_apply]

theorem uniformStepSelectedSets22_apply_profile {R : ℕ}
    (counts : Fin R → Fin 17)
    (sets : Fin R → Finset (Fin 16)) :
    (uniformStepSelectedSets22 counts).p sets =
      if ∀ step : Fin R, (sets step).card = (counts step).val then
        ∏ step : Fin R,
          1 / (Nat.choose 16 (counts step).val : ℝ)
      else 0 := by
  rw [uniformStepSelectedSets22_apply]
  by_cases h : ∀ step : Fin R,
      (sets step).card = (counts step).val
  · rw [if_pos h]
    apply Finset.prod_congr rfl
    intro step _
    rw [if_pos (h step)]
  · rw [if_neg h]
    obtain ⟨step, hs⟩ := not_forall.mp h
    apply Finset.prod_eq_zero (Finset.mem_univ step)
    rw [if_neg hs]

#print axioms uniformStepSelectedSet22_apply
#print axioms uniformStepSelectedSets22_apply
#print axioms uniformStepSelectedSets22_apply_profile

end Spin.Turbo
