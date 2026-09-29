import SpinCodes.Structured.ConcreteFixedInsertionSites

/-! Exact state-path expansion of a matrix product with diagonal gap weights. -/
noncomputable section
namespace Spin.Structured.Placement
open Finset ConcreteEncoder Matrix
attribute [local instance] Classical.propDecidable

def diagonalPathProduct {κ : Type*} [Fintype κ] [DecidableEq κ] : {Q : Nat} →
    (Fin (Q+1) → κ → ℝ) → (Fin Q → Matrix κ κ ℝ) → Matrix κ κ ℝ
  | 0, D, _ => Matrix.diagonal (D 0)
  | _+1, D, M => Matrix.diagonal (D 0) * M 0 * diagonalPathProduct (Fin.tail D) (Fin.tail M)

theorem diagonalPathProduct_expansion {κ : Type*} [Fintype κ] [DecidableEq κ] {Q : Nat}
    (D : Fin (Q+1) → κ → ℝ) (M : Fin Q → Matrix κ κ ℝ) (w : κ → ℝ) (i : κ) :
    (diagonalPathProduct D M *ᵥ w) i = ∑ ys : Fin Q → κ,
      (∏ k : Fin (Q+1), D k ((Fin.cons i ys : Fin (Q+1) → κ) k)) *
      (∏ k : Fin Q, M k ((Fin.cons i ys : Fin (Q+1) → κ) k.castSucc) (ys k)) * w ((Fin.cons i ys : Fin (Q+1) → κ) (Fin.last Q)) := by
  induction Q generalizing i with
  | zero => simp [diagonalPathProduct, Matrix.mulVec_diagonal]
  | succ Q ih =>
    simp only [diagonalPathProduct, ← Matrix.mulVec_mulVec, Matrix.mulVec_diagonal]
    change D 0 i * (∑ s, M 0 i s * (diagonalPathProduct (Fin.tail D) (Fin.tail M) *ᵥ w) s) = _
    simp_rw [ih]
    rw [sum_tuple_succ]
    simp only [Fin.prod_univ_succ, Fin.cons_zero, Fin.cons_succ, Fin.castSucc_zero,
      Fin.castSucc_succ, ← Fin.succ_last, Fin.tail, Finset.mul_sum]
    apply sum_congr rfl
    intro s hs
    apply sum_congr rfl
    intro ys hys
    ring

theorem diagonalPathProduct_list {κ : Type*} [Fintype κ] [DecidableEq κ] {Q : Nat}
    (D : Fin (Q+1) → κ → ℝ) (M : Fin Q → Matrix κ κ ℝ) :
    diagonalPathProduct D M =
      (List.ofFn (fun k => Matrix.diagonal (D k.castSucc) * M k)).prod * Matrix.diagonal (D (Fin.last Q)) := by
  induction Q with
  | zero => simp [diagonalPathProduct]
  | succ Q ih =>
    simp only [diagonalPathProduct, List.ofFn_succ, List.prod_cons, Fin.castSucc_zero, ih,
      Fin.castSucc_succ, ← Fin.succ_last, Fin.tail, Matrix.mul_assoc]

end Spin.Structured.Placement

