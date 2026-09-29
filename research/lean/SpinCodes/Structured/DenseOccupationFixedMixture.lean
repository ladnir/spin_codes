import SpinCodes.Structured.DenseOccupationFixedColumn

noncomputable section
namespace Spin.Structured.DenseOccupationFixed
open Spin.Numeric Spin.Imt SparsePolynomial

lemma probability_mem {j : Nat} (hj : j≤128) (qn qd : Int) (hd : 0<qd) :
    Fix.Mem (probability j qn qd) (Occupation.probability 128 j ((qn:ℝ)/qd)) := by
  have h := Fix.ofFrac_mem (p := (polyChoose 128 j:Int)*qn^j*(qd-qn)^(128-j))
    (show 0<qd^128 by positivity)
  unfold probability
  convert h using 1
  unfold Occupation.probability
  push_cast [polyChoose_eq]
  have hdR : (qd:ℝ) ≠ 0 := by exact_mod_cast hd.ne'
  have he : (qd:ℝ)^128 = (qd:ℝ)^j*(qd:ℝ)^(128-j) := by
    rw [← pow_add, Nat.add_sub_of_le hj]
  rw [he]
  rw [show 1-(qn:ℝ)/qd = ((qd:ℝ)-qn)/qd by field_simp]
  simp only [div_pow]
  ring

lemma column_mem (qn qd : Int) (hd : 0<qd) (zp : Nat → Fix) (z : ℝ)
    (hz : ∀ n, Fix.Mem (zp n) (z^n)) {v : FCoords} {w : Coords 5} (hv : v.Mem w) :
    (column qn qd zp v).Mem
      ((Occupation.Sparse.numericalMatrix ((qn:ℝ)/qd) z).applyCol w) := by
  rw [Occupation.Sparse.numericalMatrix, Occupation.matrix_col]
  have hrow (j : Nat) (hj : j<129) :=
    fixedColumn_mem (by omega : j≤128) (Data.weights.getD j Data.weight0) zp z hz hv
  refine ⟨?_,?_,fun i => ?_⟩
  · change Fix.Mem _ (∑ j : Fin 129, Occupation.probability 128 j ((qn:ℝ)/qd) *
      ((Occupation.fixed Occupation.Sparse.count 524287
        (Occupation.Sparse.row j (Data.weights.getD j Data.weight0) z)).applyCol w).Z)
    rw [sum_fin_range (fun j => Occupation.probability 128 j ((qn:ℝ)/qd) *
      ((Occupation.fixed Occupation.Sparse.count 524287
        (Occupation.Sparse.row j (Data.weights.getD j Data.weight0) z)).applyCol w).Z) 129]
    apply sum_mem
    intro j hj
    exact Fix.mul_mem (probability_mem (by simpa only [List.mem_range] using Nat.le_of_lt_succ (List.mem_range.mp hj)) qn qd hd)
      (hrow j (List.mem_range.mp hj)).1
  · change Fix.Mem _ (∑ j : Fin 129, Occupation.probability 128 j ((qn:ℝ)/qd) *
      ((Occupation.fixed Occupation.Sparse.count 524287
        (Occupation.Sparse.row j (Data.weights.getD j Data.weight0) z)).applyCol w).D)
    rw [sum_fin_range (fun j => Occupation.probability 128 j ((qn:ℝ)/qd) *
      ((Occupation.fixed Occupation.Sparse.count 524287
        (Occupation.Sparse.row j (Data.weights.getD j Data.weight0) z)).applyCol w).D) 129]
    apply sum_mem
    intro j hj
    exact Fix.mul_mem (probability_mem (by have := List.mem_range.mp hj; omega) qn qd hd)
      (hrow j (List.mem_range.mp hj)).2.1
  · change Fix.Mem _ (∑ j : Fin 129, Occupation.probability 128 j ((qn:ℝ)/qd) *
      ((Occupation.fixed Occupation.Sparse.count 524287
        (Occupation.Sparse.row j (Data.weights.getD j Data.weight0) z)).applyCol w).S i)
    rw [sum_fin_range (fun j => Occupation.probability 128 j ((qn:ℝ)/qd) *
      ((Occupation.fixed Occupation.Sparse.count 524287
        (Occupation.Sparse.row j (Data.weights.getD j Data.weight0) z)).applyCol w).S i) 129]
    apply sum_mem
    intro j hj
    exact Fix.mul_mem (probability_mem (by have := List.mem_range.mp hj; omega) qn qd hd)
      ((hrow j (List.mem_range.mp hj)).2.2 i)

end Spin.Structured.DenseOccupationFixed
