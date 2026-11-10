import SpinCodes.FiniteProductLaw
import Mathlib.Algebra.BigOperators.Fin

/-! The audited `ScalableRouteAverage.matrixExpect22_pi_product` proof with
the finite matrix index abstracted, so the identical chronological
averaging argument applies to both full states and eleven coordinates. -/
noncomputable section
namespace Spin.Turbo.SparseMatrixAverage22
open Finset
variable {I : Type*} [Fintype I] [DecidableEq I]

def average {A : Type*} [Fintype A] (P : Spin.FinPMF A) (F : A → Matrix I I ℝ) :
    Matrix I I ℝ := fun i j => P.expect (fun a => F a i j)

private theorem expect_const_mul {A : Type*} [Fintype A] (P : Spin.FinPMF A)
    (c : ℝ) (f : A → ℝ) : P.expect (fun a => c * f a) = c * P.expect f := by
  unfold Spin.FinPMF.expect
  rw [Finset.mul_sum]
  apply Finset.sum_congr rfl
  intro a _
  ring

private theorem expect_mul_const {A : Type*} [Fintype A] (P : Spin.FinPMF A)
    (c : ℝ) (f : A → ℝ) : P.expect (fun a => f a * c) = P.expect f * c := by
  simpa only [mul_comm] using expect_const_mul P c f

theorem average_const {A : Type*} [Fintype A] (P : Spin.FinPMF A) (K : Matrix I I ℝ) :
    average P (fun _ => K) = K := by
  ext i j
  exact Spin.FinPMF.expect_const _ _

theorem average_mul_left {A : Type*} [Fintype A]
    (P : Spin.FinPMF A) (K : Matrix I I ℝ) (F : A → Matrix I I ℝ) :
    average P (fun a => K * F a) = K * average P F := by
  ext i j
  simp only [average, Matrix.mul_apply, Spin.FinPMF.expect_sum, expect_const_mul]

theorem average_mul_right {A : Type*} [Fintype A]
    (P : Spin.FinPMF A) (F : A → Matrix I I ℝ) (K : Matrix I I ℝ) :
    average P (fun a => F a * K) = average P F * K := by
  ext i j
  simp only [average, Matrix.mul_apply, Spin.FinPMF.expect_sum, expect_mul_const]

theorem average_pi_succ {A : Type*} [Fintype A] {n : ℕ}
    (P : Fin (n + 1) → Spin.FinPMF A) (K : (Fin (n + 1) → A) → Matrix I I ℝ) :
    average (Spin.piPMF P) K = average (P 0) (fun a =>
      average (Spin.piPMF (Fin.tail P)) (fun xs => K (Fin.cons a xs))) := by
  ext i j
  unfold average Spin.FinPMF.expect
  rw [← Equiv.sum_comp (Fin.consEquiv (fun _ : Fin (n + 1) => A))]
  simp only [Fintype.sum_prod_type, Fin.consEquiv, Equiv.coe_fn_mk,
    Spin.piPMF_apply, Fin.prod_univ_succ,
    Fin.cons_zero, Fin.cons_succ, Fin.tail, Finset.mul_sum, mul_assoc]

theorem average_pi_product {A : Type*} [Fintype A] {n : ℕ}
    (P : Fin n → Spin.FinPMF A) (K : Fin n → A → Matrix I I ℝ) :
    average (Spin.piPMF P) (fun xs => (List.ofFn (fun i => K i (xs i))).prod) =
      (List.ofFn (fun i => average (P i) (K i))).prod := by
  induction n with
  | zero => simp only [List.ofFn_zero, List.prod_nil, average_const]
  | succ n ih =>
    rw [average_pi_succ]
    simp only [List.ofFn_succ, List.prod_cons, Fin.cons_zero, Fin.cons_succ]
    have he (a : A) : average (Spin.piPMF (Fin.tail P))
        (fun xs : Fin n → A => K 0 a * (List.ofFn (fun i : Fin n => K i.succ (xs i))).prod) =
        K 0 a * (List.ofFn (fun i : Fin n => average (P i.succ) (K i.succ))).prod := by
      rw [average_mul_left, ih]
      rfl
    simp only [he, average_mul_right]

#print axioms average_pi_product
end Spin.Turbo.SparseMatrixAverage22

