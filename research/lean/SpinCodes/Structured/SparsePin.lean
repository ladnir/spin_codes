import SpinCodes.Structured.SparseModel
import SpinCodes.Structured.SparseDataChecks
import SpinCodes.Structured.PolySignIdentity
import SpinCodes.Structured.SparseData.Contribution0

/-! Statement pins for the sparse bridge. These supplement the original pins.
They do not assert the still-open seven-coordinate Collatz inequality. -/

noncomputable section

namespace Spin.Structured.SparsePin

open Spin.Imt Spin.Imt.Occupation

example : Sparse.count = ![(5166 : ℝ), 110288, 293455, 110128, 5250] := rfl
example : Sparse.level = ![48, 56, 64, 72, 80] := rfl

example (α : ℝ) : Sparse.witness α =
    ⟨1, (1 + 900 * α) / 1024, fun i => (1 + Sparse.correction i * α) / 1024⟩ := rfl

example {k : ℕ} (count : Fin k → ℝ) (M : ℝ) (r : RowData k) (i : Fin k) :
    (fixed count M r).rowS i =
      ⟨r.cancelS i / 2 + min (r.momentS i) r.live / (2 * M),
       if r.live = 0 then 0 else r.momentS i / 2,
       fun h => count h * r.momentS i / (2 * M) +
         if r.live = 0 ∧ h = i then r.momentS i / 2 else 0⟩ := rfl

example (w j : ℕ) (z : ℝ) : hyperMoment w j z =
    ((List.range 129).map fun h => if h ≤ j then
      (w.choose h : ℝ) * ((128 - w).choose (j - h) : ℝ) * z ^ (w + j - 2 * h)
      else 0).sum / (128 : ℕ).choose j := rfl

example {j : ℕ} (hj : j ≤ 128) (d : SparsePolynomial.WeightData) (x : ℝ) :
    (SparsePolynomial.contribution j d 0).eval x =
      probability 128 j (x / 12500) *
      ((fixed Sparse.count 524287 (Sparse.row j d (1 - x / 6250))).applyCol
        (Sparse.witness (x / 10000))).Z := SparsePolynomial.eval_zero_contribution hj d x

example {α : ℝ} (h0 : 0 ≤ α) (h1 : α ≤ 1 / 10000) :
    1 / 2048 ≤ (Sparse.witness α).Z ∧ 1 / 2048 ≤ (Sparse.witness α).D ∧
      ∀ i, 1 / 2048 ≤ (Sparse.witness α).S i := Sparse.witness_lower h0 h1

end Spin.Structured.SparsePin
