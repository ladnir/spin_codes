import SpinCodes.Turbo.RS8StepProfileMixture22

/-!
# Factor selected four-bit values across physical steps

Once a region's latent selected set is fixed, the independent uniform
four-bit values split across its ordered 16-slot physical steps. This
identity uses complete values, including selected positions whose value
is zero, and precedes averaging over the count profile.
-/

noncomputable section

namespace Spin.Turbo

open Finset

def regionalStepWordEquiv22 (R : ℕ) :
    (Fin (R * 16) → Packet4) ≃ (Fin R → PacketStep16) where
  toFun word step slot :=
    word ((finProdFinEquiv : Fin R × Fin 16 ≃ Fin (R * 16))
      (step, slot))
  invFun words index :=
    let pair := (finProdFinEquiv : Fin R × Fin 16 ≃ Fin (R * 16)).symm index
    words pair.1 pair.2
  left_inv := by
    intro word
    funext index
    change word ((finProdFinEquiv : Fin R × Fin 16 ≃ Fin (R * 16))
      (((finProdFinEquiv : Fin R × Fin 16 ≃ Fin (R * 16)).symm index).1,
        ((finProdFinEquiv : Fin R × Fin 16 ≃ Fin (R * 16)).symm index).2)) =
      word index
    simp only [Prod.mk.eta, Equiv.apply_symm_apply]
  right_inv := by
    intro words
    funext step slot
    change words
      ((finProdFinEquiv : Fin R × Fin 16 ≃ Fin (R * 16)).symm
        ((finProdFinEquiv : Fin R × Fin 16 ≃ Fin (R * 16))
          (step, slot))).1
      ((finProdFinEquiv : Fin R × Fin 16 ≃ Fin (R * 16)).symm
        ((finProdFinEquiv : Fin R × Fin 16 ≃ Fin (R * 16))
          (step, slot))).2 = words step slot
    simp only [Equiv.symm_apply_apply]

theorem regionalStepWord_support_subset_iff {R : ℕ}
    (T : Finset (Fin (R * 16)))
    (word : Fin (R * 16) → Packet4) :
    regionalPacketSupport22 word ⊆ T ↔
      ∀ step : Fin R,
        regionalPacketSupport22 ((regionalStepWordEquiv22 R) word step) ⊆
          regionalStepSelectedSet22 T step := by
  let e : Fin R × Fin 16 ≃ Fin (R * 16) := finProdFinEquiv
  constructor
  · intro h step slot hslot
    apply Finset.mem_filter.mpr
    refine ⟨Finset.mem_univ _, ?_⟩
    apply h
    apply Finset.mem_filter.mpr
    exact ⟨Finset.mem_univ _, (Finset.mem_filter.mp hslot).2⟩
  · intro h index hindex
    let pair := e.symm index
    have hslot : pair.2 ∈
        regionalPacketSupport22 ((regionalStepWordEquiv22 R) word pair.1) := by
      apply Finset.mem_filter.mpr
      refine ⟨Finset.mem_univ _, ?_⟩
      change word (e (pair.1, pair.2)) ≠ 0
      have he : e (pair.1, pair.2) = index := by
        change e ((e.symm index).1, (e.symm index).2) = index
        simpa only [Prod.mk.eta] using e.apply_symm_apply index
      rw [he]
      exact (Finset.mem_filter.mp hindex).2
    have hs := h pair.1 hslot
    have he : e (pair.1, pair.2) = index := by
      change e ((e.symm index).1, (e.symm index).2) = index
      simpa only [Prod.mk.eta] using e.apply_symm_apply index
    have hsT : e (pair.1, pair.2) ∈ T :=
      (Finset.mem_filter.mp hs).2
    rw [he] at hsT
    exact hsT

theorem regionalSelectedGiven22_step_product {R : ℕ}
    (T : Finset (Fin (R * 16))) :
    (regionalSelectedGiven22 T).map (regionalStepWordEquiv22 R) =
      Spin.piPMF (fun step : Fin R =>
        regionalSelectedGiven22 (regionalStepSelectedSet22 T step)) := by
  apply FinPMF.ext
  funext words
  rw [FinPMF.map_equiv_apply, regionalSelectedGiven22_apply,
    Spin.piPMF_apply]
  simp_rw [regionalSelectedGiven22_apply]
  have hsub : regionalPacketSupport22
      ((regionalStepWordEquiv22 R).symm words) ⊆ T ↔
      ∀ step : Fin R,
        regionalPacketSupport22 (words step) ⊆
          regionalStepSelectedSet22 T step := by
    simpa only [Equiv.apply_symm_apply] using
      regionalStepWord_support_subset_iff T
        ((regionalStepWordEquiv22 R).symm words)
  by_cases h : ∀ step : Fin R,
      regionalPacketSupport22 (words step) ⊆
        regionalStepSelectedSet22 T step
  · rw [if_pos (hsub.mpr h)]
    simp_rw [if_pos (h _)]
    rw [regionalStepSelectedSet22_card_sum]
    have hpow :
        (16 : ℝ) ^ (∑ step : Fin R, (regionalStepSelectedSet22 T step).card) =
          ∏ step : Fin R,
            (16 : ℝ) ^ (regionalStepSelectedSet22 T step).card :=
      (Finset.prod_pow_eq_pow_sum Finset.univ
        (fun step : Fin R => (regionalStepSelectedSet22 T step).card)
        (16 : ℝ)).symm
    rw [hpow]
    simp only [one_div]
    rw [Finset.prod_inv_distrib]
  · rw [if_neg (fun hglobal => h (hsub.mp hglobal))]
    obtain ⟨step, hs⟩ := not_forall.mp h
    symm
    apply Finset.prod_eq_zero (Finset.mem_univ step)
    exact if_neg hs

#print axioms regionalSelectedGiven22_step_product

end Spin.Turbo
