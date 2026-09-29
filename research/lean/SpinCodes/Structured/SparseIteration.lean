import SpinCodes.Structured.SparseContraction
import SpinCodes.Structured.SparseNonneg

/-! The sparse certificate bounds every finite iterate of the numerical
occupation transfer, with the paper's explicit prefactor 2048. -/

noncomputable section
namespace Spin.Imt.Occupation.Sparse

theorem numericalMatrix_nonneg {β z : ℝ} (hβ0 : 0 ≤ β) (hβ1 : β ≤ 1) (hz : 0 ≤ z) :
    (numericalMatrix β z).Nonneg := by
  apply matrix_nonneg 128 count 524287 β _ hβ0 hβ1
  intro j
  exact fixed_row_nonneg (by omega) _
    (Spin.Structured.SparsePolynomial.checkData_ranges
      (Spin.Structured.SparsePolynomial.data_all_valid j)).1 hz

/-- The all-ones moment bound used after the sparse Collatz certificate. -/
theorem sparse_iterate_bound {α : ℝ} (h0 : 0 < α) (h1 : α ≤ 1 / 10000) (R : ℕ) :
    (((numericalMatrix ((4 / 5) * α) (1 - (8 / 5) * α)).apply)^[R]
      (Coords.eZ 5)).total ≤ 2048 * (1 - 96 * α) ^ R := by
  have hT := numericalMatrix_nonneg (β := (4 / 5) * α) (z := 1 - (8 / 5) * α)
    (by positivity) (by linarith) (by linarith)
  obtain ⟨hz, hd, hs⟩ := witness_lower h0.le h1
  have h := collatz_total hT (sparse_collatz h0 h1) (by norm_num : (0 : ℝ) < 1 / 2048)
    hz hd hs R
  simp only [witness, mul_one] at h
  linarith

end Spin.Imt.Occupation.Sparse
