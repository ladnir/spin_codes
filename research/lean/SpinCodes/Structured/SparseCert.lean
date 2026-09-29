/-
The sparse Collatz certificate (`eq:imt-sparse-collatz`), assembled.

The paper states

    T_occ(β,z) v(α) ≤ (1 - 96α) v(α)   for 0 < α ≤ 10⁻⁴,

with `β = (4/5)α`, `z = 1 - (8/5)α`, `v_Z = 1`, `v_i = 2⁻¹⁰(1 + h_i α)`, and
proves it by exhibiting degree-257 rational polynomial upper bounds for the
seven residuals in the variable `x = 10000α ∈ [0,1]`.

The seven polynomial sign conditions are now checked by the Lean kernel.
Each coefficient list was regenerated from `certify_imt_sparse.py` and matched
against the `coefficients_sha256` recorded in `SPARSE_EXACT.json`.

What this does *not* yet establish: that these polynomials really are upper
bounds for the residuals of `T_occ`.  That needs `T_occ` itself — the finite
IMT transfer construction — which is T9.  The numerical core is done; the
bridge from `T_occ` to these polynomials is still owed.
-/
import SpinCodes.Structured.PolyCert
import SpinCodes.Structured.SparseData.Row0
import SpinCodes.Structured.SparseData.Row1
import SpinCodes.Structured.SparseData.Row2
import SpinCodes.Structured.SparseData.Row3
import SpinCodes.Structured.SparseData.Row4
import SpinCodes.Structured.SparseData.Row5
import SpinCodes.Structured.SparseData.Row6

namespace Spin.Structured

open SparseData

/-- The seven certified residual rows, indexed as `(Z, D, S₁, …, S₅)`. -/
def sparseRows : List (List ℤ) := [row0, row1, row2, row3, row4, row5, row6]

/-- **Every residual row satisfies the positive-power-tail certificate.** -/
theorem sparseRows_cert : ∀ r ∈ sparseRows, tailBound r < 0 := by
  intro r hr
  simp only [sparseRows, List.mem_cons, List.not_mem_nil, or_false] at hr
  rcases hr with rfl | rfl | rfl | rfl | rfl | rfl | rfl
  · exact row0_cert
  · exact row1_cert
  · exact row2_cert
  · exact row3_cert
  · exact row4_cert
  · exact row5_cert
  · exact row6_cert

/-- **Every residual row is strictly negative on `[0,1]`**, hence — after
restoring the order-one zero at the origin — the residuals of the sparse
Collatz inequality are strictly negative for `0 < α ≤ 10⁻⁴`. -/
theorem sparseRows_neg {x : ℝ} (hx0 : 0 ≤ x) (hx1 : x ≤ 1) :
    ∀ r ∈ sparseRows, listEval r x < 0 :=
  fun r hr => listEval_neg_of_tailBound r (sparseRows_cert r hr) hx0 hx1

/-- With the order-one zero restored: strictly negative on `(0,1]`. -/
theorem sparseRows_residual_neg {x : ℝ} (hx0 : 0 < x) (hx1 : x ≤ 1) :
    ∀ r ∈ sparseRows, x * listEval r x < 0 :=
  fun r hr => residual_neg_of_tailBound r (sparseRows_cert r hr) hx0 hx1

end Spin.Structured
