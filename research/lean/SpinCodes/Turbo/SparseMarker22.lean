import SpinCodes.Turbo.MatrixMarker22
import SpinCodes.Turbo.SparseShellCoefficient22

/-! Audit: `MatrixMarker22.orderedMatrixCoefficient22_marker_le` discards
an exact degree constraint only after weighting each ordered choice.
For a common marker on the two selected groups, their separate shell
constraints imply the exact total degree u+v. Discarding the stronger
constraint is nonnegative, and reuses the existing marker theorem. -/
noncomputable section
namespace Spin.Turbo
open Finset
variable {I : Type*} [Fintype I] [DecidableEq I]

def sparseTotalRegionDegree22 (bits : Bool × Bool) : ℕ :=
  sparseBoolDegree22 bits.1 + sparseBoolDegree22 bits.2

private theorem sparse_list_nonneg22 (ms : List (Matrix I I ℝ))
    (hm : ∀ M ∈ ms, ∀ i j, 0 ≤ M i j) : ∀ i j, 0 ≤ ms.prod i j := by
  induction ms with
  | nil =>
    intro i j
    simp only [List.prod_nil, Matrix.one_apply]
    split_ifs <;> norm_num
  | cons M ms ih =>
    intro i j
    rw [List.prod_cons, Matrix.mul_apply]
    exact Finset.sum_nonneg fun k _ => mul_nonneg (hm M List.mem_cons_self i k)
      (ih (fun N hN => hm N (List.mem_cons_of_mem M hN)) k j)

theorem sparseTwoShellCoefficient22_le_total
    (W : Bool × Bool → Matrix I I ℝ) (hW : ∀ a i j, 0 ≤ W a i j)
    (u v : ℕ) (i j : I) :
    sparseTwoShellCoefficient22 W u v i j ≤
      orderedMatrixCoefficient sparseTotalRegionDegree22 W 64 (u + v) i j := by
  unfold sparseTwoShellCoefficient22 orderedMatrixCoefficient
  rw [Matrix.sum_apply, Matrix.sum_apply]
  apply Finset.sum_le_sum
  intro bits _
  by_cases hs : (∑ r, sparseBoolDegree22 (bits r).1) = u ∧
      (∑ r, sparseBoolDegree22 (bits r).2) = v
  · have ht : (∑ r, sparseTotalRegionDegree22 (bits r)) = u + v := by
      simp only [sparseTotalRegionDegree22, Finset.sum_add_distrib, hs.1, hs.2]
    simp only [hs, ht, and_self, ite_true, le_refl]
  · simp only [hs, ite_false, Matrix.zero_apply]
    split_ifs
    · apply sparse_list_nonneg22
      intro M hM
      obtain ⟨r, rfl⟩ := List.mem_ofFn.mp hM
      exact hW _
    · exact le_rfl

theorem sparseTwoShellCoefficient22_common_marker_le
    (W : Bool × Bool → Matrix I I ℝ) (hW : ∀ a i j, 0 ≤ W a i j)
    (t : ℝ) (ht : 0 ≤ t) (u v : ℕ) (i j : I) :
    t ^ (u + v) * sparseTwoShellCoefficient22 W u v i j ≤
      (orderedMatrixMarker22 sparseTotalRegionDegree22 W t ^ 64) i j := by
  apply (mul_le_mul_of_nonneg_left (sparseTwoShellCoefficient22_le_total W hW u v i j)
    (pow_nonneg ht (u + v))).trans
  exact orderedMatrixCoefficient22_marker_le sparseTotalRegionDegree22 W hW t ht 64 (u + v) i j

def sparseJointRegionChoice22 (R0 R1 R2 : Matrix I I ℝ) (bits : Bool × Bool) : Matrix I I ℝ :=
  if bits.1 then (if bits.2 then R2 else R1) else (if bits.2 then R1 else R0)

theorem sparseJointRegionChoice22_marker (R0 R1 R2 : Matrix I I ℝ) (t : ℝ) :
    orderedMatrixMarker22 sparseTotalRegionDegree22 (sparseJointRegionChoice22 R0 R1 R2) t =
      R0 + (2 * t) • R1 + (t ^ 2) • R2 := by
  simp only [orderedMatrixMarker22, Fintype.sum_prod_type, Fintype.sum_bool,
    sparseTotalRegionDegree22, sparseBoolDegree22, sparseJointRegionChoice22,
    Bool.false_eq_true, ite_false, ite_true, zero_add, add_zero, pow_zero, pow_one,
    one_smul]
  ext i j
  simp only [Matrix.add_apply, Matrix.smul_apply, smul_eq_mul]
  ring

#print axioms sparseTwoShellCoefficient22_common_marker_le
#print axioms sparseJointRegionChoice22_marker
end Spin.Turbo

