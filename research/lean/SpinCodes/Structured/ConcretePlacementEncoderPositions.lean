import SpinCodes.Structured.ConcretePlacementLimitNormalize
import Mathlib.Data.List.GetD

/-! Exact positions of singleton impulses in the gap-expanded input stream. -/
noncomputable section
namespace Spin.Structured.Placement
open Finset ConcreteEncoder

def impulsePosition {a : Nat} (gaps : Fin a → Nat) (i : Fin a) : Nat :=
  Fin.partialSum (fun j => gaps j + 1) i.castSucc + gaps i

@[simp] theorem impulsePosition_zero {a : Nat} (gaps : Fin (a + 1) → Nat) :
    impulsePosition gaps 0 = gaps 0 := by simp [impulsePosition]

theorem impulsePosition_succ {a : Nat} (gaps : Fin (a + 1) → Nat) (i : Fin a) :
    impulsePosition gaps i.succ = gaps 0 + 1 + impulsePosition (Fin.tail gaps) i := by
  unfold impulsePosition
  rw [Fin.castSucc_succ, Fin.partialSum_succ']
  simp only [Fin.tail_def, Nat.add_assoc]

theorem gaps_partialSum {R a : Nat} (B : BlockSubset R a) (i : Fin (a + 1)) :
    Fin.partialSum (fun j : Fin a => emptyGaps B j.castSucc + 1) i = gapLeft B i := by
  refine Fin.inductionOn i ?_ (fun j ih => ?_)
  · simp [gapLeft]
  · rw [Fin.partialSum_succ, ih]
    have hh := Nat.sub_add_cancel (gapLeft_le_right B j.castSucc)
    unfold emptyGaps at ⊢
    rw [show gapLeft B j.castSucc + (gapRight B j.castSucc - gapLeft B j.castSucc + 1) =
      gapRight B j.castSucc + 1 by omega]
    simp [gapRight, gapLeft]

theorem impulsePosition_emptyGaps {R a : Nat} (B : BlockSubset R a) (i : Fin a) :
    impulsePosition (fun j => emptyGaps B j.castSucc) i = (orderedBlocks B i).val := by
  rw [impulsePosition, gaps_partialSum]
  have h := Nat.sub_add_cancel (gapLeft_le_right B i.castSucc)
  have he : gapLeft B i.castSucc + emptyGaps B i.castSucc = gapRight B i.castSucc := by
    unfold emptyGaps
    omega
  rw [he]
  exact Fin.snoc_castSucc _ _ _

theorem impulseInputs_length {a : Nat} (gaps : Fin a → Nat) (coords : Fin a → Fin 128) (last : Nat) :
    (impulseInputs gaps coords last).length = (∑ i, gaps i) + a + last := by
  induction a with
  | zero => simp [impulseInputs]
  | succ a ih =>
    simp only [impulseInputs, List.length_append, List.length_replicate, List.length_cons, ih,
      Fin.sum_univ_succ]
    change gaps 0 + ((∑ i : Fin a, gaps i.succ) + a + last + 1) = _
    omega

theorem placement_impulseInputs_length {R a : Nat} (B : BlockSubset R a) (coords : Fin a → Fin 128) :
    (impulseInputs (fun i => emptyGaps B i.castSucc) coords (emptyGaps B (Fin.last a))).length = R := by
  rw [impulseInputs_length]
  have h := emptyGaps_sum B
  rw [Fin.sum_univ_castSucc] at h
  omega

end Spin.Structured.Placement

