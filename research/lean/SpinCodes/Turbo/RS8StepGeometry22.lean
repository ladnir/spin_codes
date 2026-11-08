import SpinCodes.Turbo.RS8ReferencePhysical22

/-!
# Ordered 16-slot views of a regional selected set

One region contains `R` consecutive physical steps, each with 16 slots.
This module only changes coordinates and counts selected positions. It
does not assume that the step subsets are independent unconditionally.
-/

noncomputable section

namespace Spin.Turbo

open Finset

def regionalStepSelectedSet22 {R : ℕ}
    (T : Finset (Fin (R * 16))) (step : Fin R) :
    Finset (Fin 16) :=
  Finset.univ.filter (fun slot =>
    (finProdFinEquiv : Fin R × Fin 16 ≃ Fin (R * 16))
      (step, slot) ∈ T)

def regionalStepSetsEquiv22 (R : ℕ) :
    Finset (Fin (R * 16)) ≃
      (Fin R → Finset (Fin 16)) where
  toFun T := regionalStepSelectedSet22 T
  invFun sets := Finset.univ.filter (fun index =>
    let pair := (finProdFinEquiv : Fin R × Fin 16 ≃ Fin (R * 16)).symm index
    pair.2 ∈ sets pair.1)
  left_inv := by
    intro T
    ext index
    simp only [Finset.mem_filter, Finset.mem_univ, true_and,
      regionalStepSelectedSet22]
    have he : (finProdFinEquiv : Fin R × Fin 16 ≃ Fin (R * 16))
        (((finProdFinEquiv : Fin R × Fin 16 ≃ Fin (R * 16)).symm index).1,
          ((finProdFinEquiv : Fin R × Fin 16 ≃ Fin (R * 16)).symm index).2) =
          index := by
      simpa only [Prod.mk.eta] using
        (Equiv.apply_symm_apply
          (finProdFinEquiv : Fin R × Fin 16 ≃ Fin (R * 16)) index)
    rw [he]
  right_inv := by
    intro sets
    funext step
    ext slot
    simp only [Finset.mem_filter, Finset.mem_univ, true_and,
      regionalStepSelectedSet22]
    simp only [Equiv.symm_apply_apply]

theorem regionalStepSelectedSet22_card_sum {R : ℕ}
    (T : Finset (Fin (R * 16))) :
    T.card = ∑ step : Fin R,
      (regionalStepSelectedSet22 T step).card := by
  let e : (Fin R × Fin 16) ≃ Fin (R * 16) := finProdFinEquiv
  have hfiber (step : Fin R) :
      (T.filter (fun index => (e.symm index).1 = step)).card =
        (regionalStepSelectedSet22 T step).card := by
    let emb : Fin 16 ↪ Fin (R * 16) :=
      ⟨fun slot => e (step, slot), by
        intro a b hab
        have h := e.injective hab
        exact congrArg Prod.snd h⟩
    have hs : T.filter (fun index => (e.symm index).1 = step) =
        (regionalStepSelectedSet22 T step).map emb := by
      ext index
      constructor
      · intro hi
        have hiT : index ∈ T := (Finset.mem_filter.mp hi).1
        have hp : (e.symm index).1 = step :=
          (Finset.mem_filter.mp hi).2
        refine Finset.mem_map.mpr ?_
        let slot := (e.symm index).2
        have he : e (step, slot) = index := by
          rw [← hp]
          simpa only [Prod.mk.eta] using e.apply_symm_apply index
        have hslot : slot ∈ regionalStepSelectedSet22 T step := by
          apply Finset.mem_filter.mpr
          refine ⟨Finset.mem_univ _, ?_⟩
          change e (step, slot) ∈ T
          rw [he]
          exact hiT
        refine ⟨slot, hslot, ?_⟩
        exact he
      · intro hi
        obtain ⟨slot, hslot, hindex⟩ := Finset.mem_map.mp hi
        subst index
        apply Finset.mem_filter.mpr
        constructor
        · exact (Finset.mem_filter.mp hslot).2
        · change (e.symm (e (step, slot))).1 = step
          simp
    rw [hs]
    exact Finset.card_map emb
  have hcount : T.card =
      ∑ step : Fin R,
        (T.filter (fun index => (e.symm index).1 = step)).card := by
    simpa using (Finset.card_eq_sum_card_fiberwise
      (s := T) (t := (Finset.univ : Finset (Fin R)))
      (f := fun index => (e.symm index).1)
      (by intro index _; exact Finset.mem_univ _))
  calc
    T.card = ∑ step : Fin R,
        (T.filter (fun index => (e.symm index).1 = step)).card := hcount
    _ = _ := by simp only [hfiber]

def regionalSelectedStepSetsLaw22 (R : ℕ)
    (S : Finset (Fin (R * 16))) :
    FinPMF (Fin R → Finset (Fin 16)) :=
  (Spin.Structured.Routing.shuffleLaw S).map
    (regionalStepSetsEquiv22 R)

/-- Under the actual uniform regional shuffle, every ordered array of
16-slot selected sets with total size `|S|` has the same probability.
The local sets are not independent before their count vector is fixed. -/
theorem regionalSelectedStepSetsLaw22_apply (R : ℕ)
    (S : Finset (Fin (R * 16)))
    (sets : Fin R → Finset (Fin 16)) :
    (regionalSelectedStepSetsLaw22 R S).p sets =
      if (∑ step : Fin R, (sets step).card) = S.card then
        1 / (Nat.choose (R * 16) S.card : ℝ) else 0 := by
  rw [regionalSelectedStepSetsLaw22,
    FinPMF.map_equiv_apply,
    Spin.Structured.Routing.shuffleLaw_apply]
  have hc : ((regionalStepSetsEquiv22 R).symm sets).card =
      ∑ step : Fin R, (sets step).card := by
    rw [regionalStepSelectedSet22_card_sum]
    apply Finset.sum_congr rfl
    intro step _
    exact congrArg Finset.card
      (congrFun ((regionalStepSetsEquiv22 R).right_inv sets) step)
  rw [hc]

#print axioms regionalStepSetsEquiv22
#print axioms regionalStepSelectedSet22_card_sum
#print axioms regionalSelectedStepSetsLaw22_apply

end Spin.Turbo
