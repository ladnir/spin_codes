import SpinCodes.Structured.SparseProgramDefs
import SpinCodes.Structured.SparseResidualsSound
noncomputable section
namespace Spin.Structured.SparsePolynomial
set_option maxRecDepth 100000
set_option maxHeartbeats 0
open Residuals
def programResidual (i : ℕ) (x : ℝ) : ℝ :=
  ((List.range 129).map (fun j => (contribution j (weightN j) i).eval x)).sum -
    (1 - 96 * (x / 10000)) * (witness i).eval x
theorem sum_chunks (chunks : List (List ℕ)) (f : ℕ → ℝ) :
    (chunks.map (fun ns => (ns.map f).sum)).sum = (chunks.flatten.map f).sum := by
  induction chunks with
  | nil => simp
  | cons ns chunks ih => simp [ih, List.sum_append]
theorem programResidual_eq_0 (x : ℝ) : programResidual 0 x = sourceResidual0 x := by
  have h := sum_chunks sourceChunks (fun j => (contribution j (weightN j) 0).eval x)
  rw [sourceChunks_cover] at h
  rw [programResidual, ← h]
  rfl
theorem programResidual_eq_1 (x : ℝ) : programResidual 1 x = sourceResidual1 x := by
  have h := sum_chunks sourceChunks (fun j => (contribution j (weightN j) 1).eval x)
  rw [sourceChunks_cover] at h
  rw [programResidual, ← h]
  rfl
theorem programResidual_eq_2 (x : ℝ) : programResidual 2 x = sourceResidual2 x := by
  have h := sum_chunks sourceChunks (fun j => (contribution j (weightN j) 2).eval x)
  rw [sourceChunks_cover] at h
  rw [programResidual, ← h]
  rfl
theorem programResidual_eq_3 (x : ℝ) : programResidual 3 x = sourceResidual3 x := by
  have h := sum_chunks sourceChunks (fun j => (contribution j (weightN j) 3).eval x)
  rw [sourceChunks_cover] at h
  rw [programResidual, ← h]
  rfl
theorem programResidual_eq_4 (x : ℝ) : programResidual 4 x = sourceResidual4 x := by
  have h := sum_chunks sourceChunks (fun j => (contribution j (weightN j) 4).eval x)
  rw [sourceChunks_cover] at h
  rw [programResidual, ← h]
  rfl
theorem programResidual_eq_5 (x : ℝ) : programResidual 5 x = sourceResidual5 x := by
  have h := sum_chunks sourceChunks (fun j => (contribution j (weightN j) 5).eval x)
  rw [sourceChunks_cover] at h
  rw [programResidual, ← h]
  rfl
theorem programResidual_eq_6 (x : ℝ) : programResidual 6 x = sourceResidual6 x := by
  have h := sum_chunks sourceChunks (fun j => (contribution j (weightN j) 6).eval x)
  rw [sourceChunks_cover] at h
  rw [programResidual, ← h]
  rfl
theorem programResidual_neg (i : Fin 7) {x : ℝ} (h0 : 0 < x) (h1 : x ≤ 1) :
    programResidual i x < 0 := by
  fin_cases i
  · rw [programResidual_eq_0]
    exact sourceResidual0_neg h0 h1
  · rw [programResidual_eq_1]
    exact sourceResidual1_neg h0 h1
  · rw [programResidual_eq_2]
    exact sourceResidual2_neg h0 h1
  · rw [programResidual_eq_3]
    exact sourceResidual3_neg h0 h1
  · rw [programResidual_eq_4]
    exact sourceResidual4_neg h0 h1
  · rw [programResidual_eq_5]
    exact sourceResidual5_neg h0 h1
  · rw [programResidual_eq_6]
    exact sourceResidual6_neg h0 h1
end Spin.Structured.SparsePolynomial
