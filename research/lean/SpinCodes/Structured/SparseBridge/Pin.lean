import SpinCodes.Structured.SparseMaximumAll
import SpinCodes.Structured.SparseProgram
import SpinCodes.Structured.SparseIteration

/-! Pins for the complete sparse polynomial program, its transfer-matrix
Collatz inequality, and the repeated-step bound. -/

noncomputable section

namespace Spin.Structured.SparseBridgePin

open SparsePolynomial

example : sourceChunks.flatten = List.range 129 := sourceChunks_cover

example (x : ℝ) (h0 : 0 ≤ x) (h1 : x ≤ 1) (j : Fin 129) :
    Spin.Imt.Occupation.Sparse.maximumMoment j (1 - x / 6250) ≤
      (arbitrary j (Data.weight j)).eval x := maximum_all x h0 h1 j

example (x : ℝ) (h0 : 0 ≤ x) (h1 : x ≤ 1) (j : Fin 129) :
    lowMaximum j (Data.weight j) (1 - x / 6250) ≤
      (selectedLowPolynomial j (Data.weight j)).eval x := lowMaximum_all x h0 h1 j

example (i : Fin 7) {x : ℝ} (h0 : 0 < x) (h1 : x ≤ 1) :
    ((List.range 129).map (fun j => (contribution j (weightN j) i).eval x)).sum <
      (1 - 96 * (x / 10000)) * (witness i).eval x := by
  have h := programResidual_neg i h0 h1
  change _ - _ < 0 at h
  linarith

example (β z : ℝ) : Spin.Imt.Occupation.Sparse.numericalMatrix β z =
    Spin.Imt.Occupation.matrix 128 Spin.Imt.Occupation.Sparse.count 524287 β
      (fun j => Spin.Imt.Occupation.Sparse.row j (Data.weight j) z) := rfl

end Spin.Structured.SparseBridgePin

namespace Spin.Imt.Occupation.Sparse

example {α : ℝ} (h0 : 0 < α) (h1 : α ≤ 1 / 10000) :
    ((numericalMatrix ((4 / 5) * α) (1 - (8 / 5) * α)).applyCol (witness α)).le
      (Coords.smul (1 - 96 * α) (witness α)) := sparse_collatz h0 h1

example {α : ℝ} (h0 : 0 < α) (h1 : α ≤ 1 / 10000) (R : ℕ) :
    (((numericalMatrix ((4 / 5) * α) (1 - (8 / 5) * α)).apply)^[R]
      (Coords.eZ 5)).total ≤ 2048 * (1 - 96 * α) ^ R := sparse_iterate_bound h0 h1 R

end Spin.Imt.Occupation.Sparse
