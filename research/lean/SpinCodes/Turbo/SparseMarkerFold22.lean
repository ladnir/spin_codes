import SpinCodes.Turbo.SparseMarker22

/-! A common positive marker bounds the complete two-group support fold.
The ordered(u,v) sum factors into one exact current-Cv weighted sum squared;
this covers every shell pair without a large numerical coefficient table. -/
noncomputable section
namespace Spin.Turbo
open Finset
variable {I : Type*} [Fintype I] [DecidableEq I]

theorem sparseTwoShell_common_marker_terminal_le
    (W : Bool × Bool → Matrix I I ℝ) (hW : ∀ a i j, 0 ≤ W a i j)
    (t : ℝ) (ht : 0 < t) (u v : ℕ) (i : I) :
    (∑ j, sparseTwoShellCoefficient22 W u v i j) ≤
      (∑ j, (orderedMatrixMarker22 sparseTotalRegionDegree22 W t ^ 64) i j) / t ^ (u + v) := by
  apply (le_div_iff₀ (pow_pos ht (u + v))).mpr
  rw [Finset.sum_mul]
  apply Finset.sum_le_sum
  intro j _
  simpa only [mul_comm] using
    sparseTwoShellCoefficient22_common_marker_le W hW t ht.le u v i j

theorem sparseTwoWeighted_common_marker_fold_le
    (W : Bool × Bool → Matrix I I ℝ) (hW : ∀ a i j, 0 ≤ W a i j)
    (t : ℝ) (ht : 0 < t) (V : Finset ℕ) (w : ℕ → ℝ)
    (hw : ∀ v ∈ V, 0 ≤ w v) (i : I) :
    (∑ u ∈ V, ∑ v ∈ V, (w u * w v) * ∑ j, sparseTwoShellCoefficient22 W u v i j) ≤
      (∑ j, (orderedMatrixMarker22 sparseTotalRegionDegree22 W t ^ 64) i j) *
        (∑ v ∈ V, w v / t ^ v) ^ 2 := by
  let T := ∑ j, (orderedMatrixMarker22 sparseTotalRegionDegree22 W t ^ 64) i j
  calc
    _ ≤ ∑ u ∈ V, ∑ v ∈ V, (w u * w v) * (T / t ^ (u + v)) := by
      apply Finset.sum_le_sum
      intro u hu
      apply Finset.sum_le_sum
      intro v hv
      exact mul_le_mul_of_nonneg_left
        (sparseTwoShell_common_marker_terminal_le W hW t ht u v i)
        (mul_nonneg (hw u hu) (hw v hv))
    _ = ∑ u ∈ V, ∑ v ∈ V, T * ((w u / t ^ u) * (w v / t ^ v)) := by
      apply Finset.sum_congr rfl
      intro u _
      apply Finset.sum_congr rfl
      intro v _
      simp only [div_eq_mul_inv, pow_add, mul_inv_rev]
      ring
    _ = T * (∑ u ∈ V, ∑ v ∈ V, (w u / t ^ u) * (w v / t ^ v)) := by
      rw [Finset.mul_sum]
      apply Finset.sum_congr rfl
      intro u _
      rw [Finset.mul_sum]
    _ = T * ((∑ v ∈ V, w v / t ^ v) * (∑ v ∈ V, w v / t ^ v)) := by
      rw [← Finset.sum_mul_sum]
    _ = _ := by rw [pow_two]

#print axioms sparseTwoWeighted_common_marker_fold_le
end Spin.Turbo
