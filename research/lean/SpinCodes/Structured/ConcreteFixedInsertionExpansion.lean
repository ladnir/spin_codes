import SpinCodes.Structured.ConcreteFixedTwoMoment

/-! Exact finite subset expansion for impulses inserted at all potential sites. -/
noncomputable section
namespace Spin.Structured.Placement
open Finset ConcreteEncoder ConcreteMarked
attribute [local instance] Classical.propDecidable

theorem finite_matrix_sum_product {κ α : Type*} [Fintype κ] [DecidableEq κ] [Fintype α] {b : Nat}
    (K : Fin b → α → Matrix κ κ ℝ) :
    (∑ xs : Fin b → α, (List.ofFn (fun i => K i (xs i))).prod) =
      (List.ofFn (fun i => ∑ x, K i x)).prod := by
  induction b with
  | zero => simp
  | succ b ih =>
    have hs (f : (Fin (b+1) → α) → Matrix κ κ ℝ) :
        (∑ xs, f xs) = ∑ x, ∑ xs : Fin b → α, f (Fin.cons x xs) := by
      rw [← Equiv.sum_comp (Fin.consEquiv (fun _ : Fin (b+1) => α)), Fintype.sum_prod_type]
      rfl
    rw [hs]
    simp only [List.ofFn_succ, List.prod_cons, Fin.cons_zero, Fin.cons_succ]
    simp_rw [← Finset.mul_sum, ih]
    rw [← Finset.sum_mul]

theorem finite_matrix_smul_product {κ : Type*} [Fintype κ] [DecidableEq κ] {b : Nat}
    (c : Fin b → ℝ) (K : Fin b → Matrix κ κ ℝ) :
    (List.ofFn (fun i => c i • K i)).prod = (∏ i, c i) • (List.ofFn K).prod := by
  induction b with
  | zero => simp
  | succ b ih =>
    simp only [List.ofFn_succ, List.prod_cons, Fin.prod_univ_succ]
    rw [ih, smul_mul_smul_comm]

def fullGapProduct {Q : Nat} (γ : ℝ) (P : Matrix (Fin 2) (Fin 2) ℝ)
    (gaps u : Fin Q → ℝ) (last : ℝ) : Matrix (Fin 2) (Fin 2) ℝ :=
  (List.ofFn (fun i => timeEmpty γ (gaps i) * (1 + u i • P))).prod * timeEmpty γ last

def selectedGapProduct {Q : Nat} (γ : ℝ) (P : Matrix (Fin 2) (Fin 2) ℝ)
    (gaps : Fin Q → ℝ) (last : ℝ) (A : Finset (Fin Q)) : Matrix (Fin 2) (Fin 2) ℝ :=
  (List.ofFn (fun i => timeEmpty γ (gaps i) * (if i ∈ A then P else 1))).prod * timeEmpty γ last

theorem fullGapProduct_subset_expansion {Q : Nat} (γ : ℝ) (P : Matrix (Fin 2) (Fin 2) ℝ)
    (gaps u : Fin Q → ℝ) (last : ℝ) :
    fullGapProduct γ P gaps u last = ∑ A : Finset (Fin Q), markCoefficient u A • selectedGapProduct γ P gaps last A := by
  have he (i : Fin Q) : timeEmpty γ (gaps i) * (1 + u i • P) =
      ∑ bit : Bool, timeEmpty γ (gaps i) * (if bit then u i • P else 1) := by
    simp [Matrix.mul_add, add_comm]
  unfold fullGapProduct
  simp only [he]
  rw [← finite_matrix_sum_product, Finset.sum_mul]
  rw [← Equiv.sum_comp (supportEquiv Q)]
  apply sum_congr rfl
  intro bits _
  have hf (i : Fin Q) : timeEmpty γ (gaps i) * (if bits i then u i • P else 1) =
      (if bits i then u i else 1) • (timeEmpty γ (gaps i) * (if bits i then P else 1)) := by
    cases h : bits i <;> simp [h, Matrix.mul_smul]
  simp only [hf, finite_matrix_smul_product, Matrix.smul_mul]
  have hc : (∏ i, if bits i then u i else 1) = markCoefficient u (support bits) := by
    rw [markCoefficient_eq_prod_ite]
    apply prod_congr rfl
    intro i _
    simp [support]
  rw [hc]
  simp only [supportEquiv, Equiv.coe_fn_mk, selectedGapProduct, mem_support]

theorem timeEmpty_add (γ s t : ℝ) : timeEmpty γ (s+t) = timeEmpty γ s * timeEmpty γ t := by
  ext i j
  fin_cases i <;> fin_cases j <;> simp [timeEmpty, Matrix.mul_apply, Fin.sum_univ_two]
  rw [← Real.exp_add]
  congr 1
  ring

theorem timeEmpty_zero (γ : ℝ) : timeEmpty γ 0 = 1 := by
  ext i j
  fin_cases i <;> fin_cases j <;> simp [timeEmpty, Matrix.one_apply]

theorem timeEmpty_list_prod (γ : ℝ) (gaps : List ℝ) :
    (gaps.map (timeEmpty γ)).prod = timeEmpty γ gaps.sum := by
  induction gaps with
  | nil => simp [timeEmpty_zero]
  | cons g gaps ih => simp only [List.map_cons, List.prod_cons, List.sum_cons, ih, timeEmpty_add]

end Spin.Structured.Placement

