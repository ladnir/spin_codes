import SpinCodes.Structured.ConcretePlacementSimplexSum

/-! Sorted occupied blocks determine the exact empty gaps and a normalized simplex point. -/
noncomputable section
namespace Spin.Structured.Placement
open Finset Routing
attribute [local instance] Classical.propDecidable

def gapRight {R a : Nat} (B : BlockSubset R a) : Fin (a + 1) → Nat :=
  Fin.snoc (fun i => (orderedBlocks B i).val) R

def gapLeft {R a : Nat} (B : BlockSubset R a) : Fin (a + 1) → Nat :=
  Fin.cons 0 (fun i => (orderedBlocks B i).val + 1)

def emptyGaps {R a : Nat} (B : BlockSubset R a) (i : Fin (a + 1)) : Nat :=
  gapRight B i - gapLeft B i

theorem gapLeft_le_right {R a : Nat} (B : BlockSubset R a) (i : Fin (a + 1)) :
    gapLeft B i ≤ gapRight B i := by
  refine Fin.cases ?_ (fun j => ?_) i
  · simp [gapLeft]
  · simp only [gapLeft, Fin.cons_succ]
    by_cases hj : j.val + 1 = a
    · have he : j.succ = Fin.last a := Fin.ext (by simpa using hj)
      rw [he, gapRight, Fin.snoc_last]
      exact (orderedBlocks B j).isLt
    · have hjlt : j.val + 1 < a := by have := j.isLt; omega
      let k : Fin a := ⟨j.val + 1, hjlt⟩
      have he : j.succ = k.castSucc := Fin.ext rfl
      rw [he, gapRight, Fin.snoc_castSucc]
      have hm := orderedBlocks_strictMono B (show j < k from by change j.val < j.val + 1; omega)
      exact hm

theorem emptyGaps_sum {R a : Nat} (B : BlockSubset R a) :
    (∑ i, emptyGaps B i) + a = R := by
  have he : (∑ i, emptyGaps B i) + (∑ i, gapLeft B i) = ∑ i, gapRight B i := by
    rw [← sum_add_distrib]
    apply sum_congr rfl
    intro i _
    exact Nat.sub_add_cancel (gapLeft_le_right B i)
  have hl : (∑ i, gapLeft B i) = (∑ i, (orderedBlocks B i).val) + a := by
    simp [gapLeft, Fin.sum_cons, sum_add_distrib]
  have hr : (∑ i, gapRight B i) = (∑ i, (orderedBlocks B i).val) + R := by
    exact Fin.sum_snoc _ _
  rw [hl, hr] at he
  omega

theorem emptyGaps_le {R a : Nat} (B : BlockSubset R a) (i : Fin (a + 1)) :
    emptyGaps B i ≤ R := by
  have h := Finset.single_le_sum (fun j _ => Nat.zero_le (emptyGaps B j)) (Finset.mem_univ i)
  have hs := emptyGaps_sum B
  omega

def simplexGaps {R a : Nat} (B : BlockSubset R a) (i : Fin (a + 1)) : ℝ :=
  ((emptyGaps B i : ℝ) + if i = Fin.last a then (a : ℝ) else 0) / R

theorem simplexGaps_nonneg {R a : Nat} (B : BlockSubset R a) (i : Fin (a + 1)) :
    0 ≤ simplexGaps B i := by unfold simplexGaps; positivity

theorem simplexGaps_sum {R a : Nat} (hR : 0 < R) (B : BlockSubset R a) :
    ∑ i, simplexGaps B i = 1 := by
  unfold simplexGaps
  rw [← sum_div, sum_add_distrib]
  simp only [Finset.sum_ite_eq', mem_univ, ite_true]
  have hs : (∑ i, (emptyGaps B i : ℝ)) + a = (R : ℝ) := by
    exact_mod_cast emptyGaps_sum B
  rw [hs]
  exact div_self (by exact_mod_cast Nat.ne_of_gt hR)

theorem simplexGaps_distance {R a : Nat} (B : BlockSubset R a) :
    ∑ i, |simplexGaps B i - (emptyGaps B i : ℝ) / R| = (a : ℝ) / R := by
  have he (i : Fin (a + 1)) :
      |simplexGaps B i - (emptyGaps B i : ℝ) / R| =
        if i = Fin.last a then (a : ℝ) / R else 0 := by
    unfold simplexGaps
    split_ifs with hi
    · rw [add_div, add_sub_cancel_left, abs_of_nonneg (by positivity)]
    · simp
  simp_rw [he]
  simp

end Spin.Structured.Placement
