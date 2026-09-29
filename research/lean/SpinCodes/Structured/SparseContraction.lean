import SpinCodes.Structured.SparseColumnBounds
import SpinCodes.Structured.SparseMaximumAll
import SpinCodes.Structured.SparseProgram

/-! The sparse Collatz inequality for the occupation matrix defined by the
frozen spectra and fiber data. Identifying these data with the concrete binary
maps and proving their one-step domination are separate obligations. -/

noncomputable section
namespace Spin.Structured.SparsePolynomial

open Spin.Imt Spin.Imt.Occupation

theorem sum_fin_range (f : ℕ → ℝ) (n : ℕ) :
    (∑ j : Fin n, f j) = ((List.range n).map f).sum := by
  rw [Fin.sum_univ_eq_sum_range]
  induction n with
  | zero => simp
  | succ n ih => simp [Finset.sum_range_succ, List.range_succ, ih]

theorem contribution_sum_lt (i : Fin 7) {x : ℝ} (h0 : 0 < x) (h1 : x ≤ 1) :
    (∑ j : Fin 129, (contribution j (Data.weight j) i).eval x) <
      (1 - 96 * (x / 10000)) * (witness i).eval x := by
  have h := programResidual_neg i h0 h1
  change _ - _ < 0 at h
  have he : (∑ j : Fin 129, (contribution j (Data.weight j) i).eval x) =
      ((List.range 129).map fun j => (contribution j (weightN j) i).eval x).sum :=
    sum_fin_range (fun j => (contribution j (weightN j) i).eval x) 129
  rw [he]
  linarith

def actionColumn (j : ℕ) (d : WeightData) (x : ℝ) : Coords 5 where
  Z := (action j d 0).eval x
  D := (action j d 1).eval x
  S i := (action j d (i + 2)).eval x

theorem fixed_col_le_actionColumn (j : Fin 129) {x : ℝ} (h0 : 0 ≤ x) (h1 : x ≤ 1) :
    ((fixed Sparse.count 524287 (Sparse.row j (Data.weight j) (1 - x / 6250))).applyCol
      (Sparse.witness (x / 10000))).le (actionColumn j (Data.weight j) x) := by
  have hj : (j : ℕ) ≤ 128 := by omega
  refine ⟨(eval_action_zero hj (Data.weight j) x).ge, ?_, fun i => ?_⟩
  · exact fixed_col_D_le_action hj _ x h0 (data_all_valid j)
      (maximum_all x h0 h1 j) (lowMaximum_all x h0 h1 j)
  · exact fixed_col_S_le_action hj _ i x (data_all_valid j)

theorem numericalMatrix_col_le {x : ℝ} (h0 : 0 < x) (h1 : x ≤ 1) :
    ((Sparse.numericalMatrix (x / 12500) (1 - x / 6250)).applyCol
      (Sparse.witness (x / 10000))).le
      (Coords.smul (1 - 96 * (x / 10000)) (Sparse.witness (x / 10000))) := by
  have hb0 : 0 ≤ x / 12500 := by positivity
  have hb1 : x / 12500 ≤ 1 := by linarith
  have hh := matrix_col_le 128 Sparse.count 524287 (x / 12500)
    (fun j => Sparse.row j (Data.weight j) (1 - x / 6250))
    (Sparse.witness (x / 10000)) (fun j => actionColumn j (Data.weight j) x)
    hb0 hb1 (fun j => fixed_col_le_actionColumn j h0.le h1)
  have he (j : Fin 129) (i : ℕ) :
      Occupation.probability 128 j (x / 12500) * (action j (Data.weight j) i).eval x =
        (contribution j (Data.weight j) i).eval x := by
    rw [contribution, RatPoly.eval_mul, eval_probability]
  simp only [Coords.le, Coords.sum, Coords.smul, actionColumn, he] at hh
  refine ⟨hh.1.trans ?_, hh.2.1.trans ?_, fun i => (hh.2.2 i).trans ?_⟩
  · simpa only [eval_witness_Z, Coords.smul] using (contribution_sum_lt ⟨0, by decide⟩ h0 h1).le
  · simpa only [eval_witness_D, Coords.smul] using (contribution_sum_lt ⟨1, by decide⟩ h0 h1).le
  · simpa only [eval_witness_S, Coords.smul] using
      (contribution_sum_lt ⟨i + 2, by omega⟩ h0 h1).le

end Spin.Structured.SparsePolynomial

namespace Spin.Imt.Occupation.Sparse

/-- Equation `imt-sparse-collatz`, for the numerical occupation matrix. -/
theorem sparse_collatz {α : ℝ} (h0 : 0 < α) (h1 : α ≤ 1 / 10000) :
    ((numericalMatrix ((4 / 5) * α) (1 - (8 / 5) * α)).applyCol (witness α)).le
      (Coords.smul (1 - 96 * α) (witness α)) := by
  have h := Spin.Structured.SparsePolynomial.numericalMatrix_col_le
    (x := 10000 * α) (by positivity) (by linarith)
  convert h using 1 <;> congr 2 <;> ring

end Spin.Imt.Occupation.Sparse
