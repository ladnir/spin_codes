import SpinCodes.UniformMixer.QuadraticRankSoundness

/-!
# Packed rows are the quadratic polar matrix

The fast evaluator stores each of the six polar equations as a six-bit
row mask. This module identifies that representation with the formal
polar linear map for every packed quadratic coefficient vector.
-/

namespace Spin.UniformMixer

theorem rowEval_onehot (j : Fin 6) (x : State6) :
    rowEval (2 ^ j.val) x = x j := by
  unfold rowEval
  simp only [Nat.testBit_two_pow, bitScalar_decide_eq]
  simp_rw [Fin.val_inj]
  simp

theorem pairLeft_ne_pairRight (p : Fin 15) :
    pairLeft p ≠ pairRight p := by
  fin_cases p <;> decide

theorem rowEval_foldl_xor (rows : List Nat) (initial : Nat)
    (x : State6) :
    rowEval (rows.foldl Nat.xor initial) x =
      rowEval initial x + (rows.map (fun row => rowEval row x)).sum := by
  induction rows generalizing initial with
  | nil => simp
  | cons row rest ih =>
      rw [List.foldl_cons, ih, rowEval_xor]
      simp [add_assoc]

theorem polarRowContribution_eval (q : Nat) (i : Fin 6)
    (p : Fin 15) (x : State6) :
    rowEval (polarRowContribution q i p) x =
      bitScalar (q.testBit p.val) *
        (x (pairLeft p) * basis6 i (pairRight p) +
          basis6 i (pairLeft p) * x (pairRight p)) := by
  by_cases hq : q.testBit p.val
  · by_cases hl : pairLeft p = i <;>
      by_cases hr : pairRight p = i
    · exact False.elim ((pairLeft_ne_pairRight p) (hl.trans hr.symm))
    all_goals simp [polarRowContribution, hq, hl, hr,
      basis6, rowEval_onehot, rowEval_zero, bitScalar, mul_comm]
  · simp [polarRowContribution, hq, rowEval_zero, bitScalar]

theorem polarRowMask_eval (q : Nat) (i : Fin 6) (x : State6) :
    rowEval (polarRowMask q i) x =
      quadraticPolarLinear
        (fun p : Fin 15 => bitScalar (q.testBit p.val))
        x (basis6 i) := by
  rw [polarRowMask, rowEval_foldl_xor]
  simp only [rowEval_zero, zero_add, List.map_ofFn, List.sum_ofFn,
    quadraticPolarLinear, LinearMap.coe_mk, AddHom.coe_mk]
  apply Finset.sum_congr rfl
  intro p _
  change rowEval (polarRowContribution q i p) x = _
  rw [polarRowContribution_eval]
  ring

theorem polarRows_vanish_iff (q : Fin (2 ^ 15)) (x : State6) :
    rowsVanish (polarRows q.val) x ↔
      x ∈ quadraticRadical (packedQuadratic q) := by
  rw [quadraticRadical_basis_mem]
  constructor
  · intro h i
    have hi : rowEval (polarRowMask q.val i) x = 0 :=
      h _ (List.mem_ofFn.mpr ⟨i, rfl⟩)
    rwa [polarRowMask_eval] at hi
  · intro h row hrow
    obtain ⟨i, rfl⟩ := List.mem_ofFn.mp hrow
    rw [polarRowMask_eval]
    exact h i

end Spin.UniformMixer
