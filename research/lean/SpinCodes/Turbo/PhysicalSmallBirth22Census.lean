import SpinCodes.Turbo.PhysicalSmallBirth22Checks
import SpinCodes.Turbo.PhysicalSmallBirth22CensusDefs

noncomputable section
namespace Spin.Turbo
open Finset Spin.UniformMixer Spin.Structured

private theorem smallBirthCountOfFn22 {n : Nat} (f : Fin n → Nat) (b : Nat) :
    (List.ofFn f).countP (fun x => x == b) = ∑ i, if f i = b then 1 else 0 := by
  induction n with
  | zero => simp
  | succ n ih =>
    simp only [List.ofFn_succ, List.countP_cons, Fin.sum_univ_succ]
    rw [ih]
    simp only [beq_iff_eq]
    omega

private theorem smallBirthGetDOfFn22 {n : Nat} (f : Fin n → Nat) (b : Nat) (hb : b < n) :
    (List.ofFn f).getD b 0 = f ⟨b,hb⟩ := by
  rw [List.getD_eq_getElem (List.ofFn f) 0 (n := b) (by simpa only [List.length_ofFn] using hb),
    List.getElem_ofFn]
private theorem smallBirthHistogram22_count (xs : List Nat) (b : Nat) (hb : b < 99) :
    (smallBirthCodeHistogram22 xs).getD b 0 = xs.countP (fun x => x == b) := by
  rw [smallBirthCodeHistogram22, MapSpectrum.foldl_bump_count _ _ _ _ (by simpa using hb)]
  simp only [List.getD_replicate 0 hb, Nat.zero_add]

theorem smallBirthOneHistogram22_count (b : Nat) (hb : b < 99) :
    smallBirthOneHistogram22.getD b 0 =
      ∑ a : SmallBirthOne22, if smallBirthOneCode22 a = b then 1 else 0 := by
  rw [← smallBirthOneHistogram22_checked, smallBirthHistogram22_count _ _ hb]
  simp only [smallBirthOneCodes22, List.countP_flatten, List.map_ofFn, List.sum_ofFn,
    Function.comp_def, smallBirthCountOfFn22, Fintype.sum_prod_type]

private theorem smallBirthSupportTwo22_sum (f : SmallBirthSupportTwo22 → Nat) :
    (∑ ab : Fin 16 × Fin 16, if h : ab.1 < ab.2 then f ⟨ab,h⟩ else 0) = ∑ ab, f ab := by
  let g : Fin 16 × Fin 16 → Nat := fun ab => if h : ab.1 < ab.2 then f ⟨ab,h⟩ else 0
  have h := Fintype.sum_subtype_add_sum_subtype (fun ab : Fin 16 × Fin 16 => ab.1 < ab.2) g
  have hy : (∑ ab : SmallBirthSupportTwo22, g ab.val) = ∑ ab, f ab := by
    apply sum_congr rfl
    intro ab _
    simp only [g, dite_eq_left ab.property]
  have hn : (∑ ab : {ab : Fin 16 × Fin 16 // ¬ ab.1 < ab.2}, g ab.val) = 0 := by
    apply sum_eq_zero
    intro ab _
    simp only [g, dite_eq_right ab.property]
  rw [hy, hn, add_zero] at h
  exact h.symm

theorem smallBirthTwoShardHistogram22_count (a : Fin 16) (b : Nat) (hb : b < 99) :
    (smallBirthTwoShardHistograms22 a).getD b 0 =
      ∑ t : Fin 16, if a < t then
        ∑ u : Fin 15, ∑ v : Fin 15, if smallBirthTwoCode22 a t u v = b then 1 else 0 else 0 := by
  rw [smallBirthTwoShardHistograms22_sum_checked a ⟨b,hb⟩]
  apply sum_congr rfl
  intro t _
  split_ifs with h
  · rw [← smallBirthPairHistograms22_checked a t h, smallBirthHistogram22_count _ _ hb]
    simp only [smallBirthPairCodes22, List.countP_flatten, List.map_ofFn, List.sum_ofFn,
      Function.comp_def, smallBirthCountOfFn22]
  · rfl
theorem smallBirthTwoHistogram22_count (b : Nat) (hb : b < 99) :
    smallBirthTwoHistogram22.getD b 0 =
      ∑ a : SmallBirthTwo22,
        if smallBirthTwoCode22 a.1.val.1 a.1.val.2 a.2.1 a.2.2 = b then 1 else 0 := by
  have he := congrArg (fun xs : List Nat => xs.getD b 0) smallBirthTwoHistogram22_sum_checked
  have he' : (∑ a : Fin 16, (smallBirthTwoShardHistograms22 a).getD b 0) =
      smallBirthTwoHistogram22.getD b 0 := by
    have hget := smallBirthGetDOfFn22
      (fun k : Fin 99 => ∑ a : Fin 16, (smallBirthTwoShardHistograms22 a).getD k.val 0) b hb
    exact hget.symm.trans he
  rw [← he']
  simp_rw [smallBirthTwoShardHistogram22_count _ _ hb]
  rw [Fintype.sum_prod_type]
  simp_rw [Fintype.sum_prod_type]
  rw [← smallBirthSupportTwo22_sum]
  rw [Fintype.sum_prod_type]
  apply sum_congr rfl
  intro x _
  apply sum_congr rfl
  intro y _
  by_cases hxy : x < y
  · simp only [hxy, ite_true, dite_true]
  · simp only [hxy, ite_false, dite_false]

theorem smallBirthOneLayer22_count (b : Nat) (hb : b < 99) :
    (univ.filter (fun x : PacketStepLayer22 1 =>
      9 * (physicalBirthIndex22 (stepFeedback (packetStepInput22 x.val))).val +
        stepWeight (packetStepInput22 x.val) 0 = b)).card = smallBirthOneHistogram22.getD b 0 := by
  rw [card_filter]
  rw [← Equiv.sum_comp smallBirthOneEquiv22]
  change (∑ a : SmallBirthOne22, if
    9 * (physicalBirthIndex22 (stepFeedback (packetStepInput22 (smallBirthOneWord22 a)))).val +
      stepWeight (packetStepInput22 (smallBirthOneWord22 a)) 0 = b then 1 else 0) = _
  simp_rw [← smallBirthOneCode22_actual]
  exact (smallBirthOneHistogram22_count b hb).symm

theorem smallBirthTwoLayer22_count (b : Nat) (hb : b < 99) :
    (univ.filter (fun x : PacketStepLayer22 2 =>
      9 * (physicalBirthIndex22 (stepFeedback (packetStepInput22 x.val))).val +
        stepWeight (packetStepInput22 x.val) 0 = b)).card = smallBirthTwoHistogram22.getD b 0 := by
  rw [card_filter]
  rw [← Equiv.sum_comp smallBirthTwoEquiv22]
  change (∑ a : SmallBirthTwo22, if
    9 * (physicalBirthIndex22 (stepFeedback (packetStepInput22 (smallBirthTwoWord22 a)))).val +
      stepWeight (packetStepInput22 (smallBirthTwoWord22 a)) 0 = b then 1 else 0) = _
  simp_rw [← smallBirthTwoCode22_actual]
  exact (smallBirthTwoHistogram22_count b hb).symm
#print axioms smallBirthOneHistogram22_count
#print axioms smallBirthTwoHistogram22_count
end Spin.Turbo





