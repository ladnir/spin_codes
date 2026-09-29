import SpinCodes.Structured.ConcretePlacementEncoderPositions

/-! The gap-expanded experiment is exactly the serialized singleton-support experiment. -/
noncomputable section
namespace Spin.Structured.Placement
open Finset ConcreteEncoder
attribute [local instance] Classical.propDecidable

theorem mem_impulseInputs_getD {a : Nat} (gaps : Fin a → Nat) (coords : Fin a → Fin 128)
    (last n : Nat) (c : Fin 128) :
    c ∈ (impulseInputs gaps coords last).getD n ∅ ↔
      ∃ i, n = impulsePosition gaps i ∧ c = coords i := by
  induction a generalizing n with
  | zero => simp [impulseInputs, List.getD_eq_getElem?_getD]
  | succ a ih =>
    rw [Fin.exists_fin_succ]
    simp only [impulsePosition_zero, impulsePosition_succ]
    rw [impulseInputs]
    by_cases hn : n < gaps 0
    · rw [List.getD_append _ _ _ _ (by simpa using hn), List.getD_replicate (∅ : Input) hn]
      simp only [Finset.notMem_empty, false_iff, not_or, not_exists, not_and]
      constructor
      · omega
      · intro i hi; omega
    · have hgn : gaps 0 ≤ n := by omega
      rw [List.getD_append_right _ _ _ _ (by simpa using hgn), List.length_replicate]
      by_cases he : n = gaps 0
      · subst n
        simp only [Nat.sub_self, List.getD_cons_zero, Finset.mem_singleton, true_and]
        constructor
        · exact Or.inl
        · rintro (h | ⟨i, hi, _⟩)
          · exact h
          · omega
      · have hp : 0 < n - gaps 0 := by omega
        obtain ⟨k, hk⟩ := Nat.exists_eq_succ_of_ne_zero (Nat.ne_of_gt hp)
        rw [hk, List.getD_cons_succ, ih]
        constructor
        · rintro ⟨i, hi, hc⟩
          exact Or.inr ⟨i, by omega, hc⟩
        · rintro (⟨h, _⟩ | ⟨i, hi, hc⟩)
          · exact (he h).elim
          · exact ⟨i, by omega, hc⟩

def serializedInput {R : Nat} (T : Finset (Fin (128 * R))) (b : Fin R) : Input :=
  univ.filter (fun c => blockBit b c ∈ T)

theorem serializedInput_singletonSupport {R a : Nat} (blocks : Fin a ↪ Fin R)
    (coords : Fin a → Fin 128) (b : Fin R) (c : Fin 128) :
    c ∈ serializedInput (singletonSupport blocks coords) b ↔
      ∃ i, b.val = (blocks i).val ∧ c = coords i := by
  simp only [serializedInput, mem_filter, mem_univ, true_and, singletonSupport, mem_image]
  constructor
  · rintro ⟨i, hi⟩
    have hb := congrArg (fun x : Fin (128 * R) => x.val / 128) hi
    have hc := congrArg (fun x : Fin (128 * R) => x.val % 128) hi
    simp only [blockBit_div] at hb
    simp only [blockBit_mod] at hc
    exact ⟨i, hb.symm, Fin.ext hc.symm⟩
  · rintro ⟨i, hb, rfl⟩
    exact ⟨i, by rw [Fin.ext hb]⟩

theorem serializedInput_eq_impulseInputs_getD {R a : Nat} (B : BlockSubset R a)
    (coords : Fin a → Fin 128) (b : Fin R) :
    serializedInput (singletonSupport (orderedBlocks B) coords) b =
      (impulseInputs (fun i => emptyGaps B i.castSucc) coords (emptyGaps B (Fin.last a))).getD b.val ∅ := by
  ext c
  rw [serializedInput_singletonSupport, mem_impulseInputs_getD]
  simp only [impulsePosition_emptyGaps]

theorem serializedInput_list_eq_impulseInputs {R a : Nat} (B : BlockSubset R a)
    (coords : Fin a → Fin 128) :
    List.ofFn (serializedInput (singletonSupport (orderedBlocks B) coords)) =
      impulseInputs (fun i => emptyGaps B i.castSucc) coords (emptyGaps B (Fin.last a)) := by
  apply List.ext_getElem
  · simp [placement_impulseInputs_length]
  · intro n hn hl
    have hnR : n < R := by simpa using hn
    rw [List.getElem_ofFn]
    rw [serializedInput_eq_impulseInputs_getD]
    exact List.getD_eq_getElem _ _ hl

end Spin.Structured.Placement

