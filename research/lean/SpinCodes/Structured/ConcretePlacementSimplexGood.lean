import SpinCodes.Structured.ConcretePlacementSimplexGaps

noncomputable section
namespace Spin.Structured.Placement
open Finset Routing
attribute [local instance] Classical.propDecidable

theorem good_emptyGaps_lower {R a : Nat} (B : BlockSubset R a) (H : Nat) (hHR : H ≤ R)
    (hg : ¬ badBlockIndices (orderedBlocks B) H) (i : Fin (a + 1)) : H ≤ emptyGaps B i := by
  have hb (j : Fin a) : H ≤ (orderedBlocks B j).val ∧ (orderedBlocks B j).val + H < R := by
    constructor
    · by_contra h; exact hg (Or.inl ⟨j, Or.inl (by omega)⟩)
    · by_contra h; exact hg (Or.inl ⟨j, Or.inr (by omega)⟩)
  have hp (j k : Fin a) (hjk : j < k) : (orderedBlocks B j).val + H < (orderedBlocks B k).val := by
    have hm : (orderedBlocks B j).val < (orderedBlocks B k).val := orderedBlocks_strictMono B hjk
    have hd : H < Nat.dist (orderedBlocks B j).val (orderedBlocks B k).val := by
      by_contra h
      exact hg (Or.inr ⟨j, k, ne_of_lt hjk, by omega⟩)
    unfold Nat.dist at hd
    omega
  have hh : gapLeft B i + H ≤ gapRight B i := by
    refine Fin.cases ?_ (fun j => ?_) i
    · cases a with
      | zero => simpa [gapLeft, gapRight, Fin.snoc] using hHR
      | succ a => simpa [gapLeft, gapRight, Fin.snoc_apply_zero] using (hb 0).1
    · simp only [gapLeft, Fin.cons_succ]
      by_cases hj : j.val + 1 = a
      · have he : j.succ = Fin.last a := Fin.ext (by simpa using hj)
        rw [he, gapRight, Fin.snoc_last]
        have := (hb j).2
        omega
      · have hjlt : j.val + 1 < a := by have := j.isLt; omega
        let k : Fin a := ⟨j.val + 1, hjlt⟩
        have he : j.succ = k.castSucc := Fin.ext rfl
        rw [he, gapRight, Fin.snoc_castSucc]
        have := hp j k (by change j.val < j.val + 1; omega)
        omega
  unfold emptyGaps
  omega

end Spin.Structured.Placement
