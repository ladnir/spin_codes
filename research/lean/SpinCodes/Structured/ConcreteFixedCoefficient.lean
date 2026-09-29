import SpinCodes.Structured.ConcreteFixedFairProduct

/-! Exact fugacity expansion for deterministic marked-region input patterns. -/
noncomputable section
namespace Spin.Structured.Placement
open Finset Routing ConcreteEncoder FiniteKernel Filter
attribute [local instance] Classical.propDecidable

theorem matrix_sum_tuple_succ {α : Type*} [Fintype α] {b : Nat}
    (f : (Fin (b + 1) → α) → Matrix State State ℝ) :
    (∑ xs, f xs) = ∑ x, ∑ xs : Fin b → α, f (Fin.cons x xs) := by
  rw [← Equiv.sum_comp (Fin.consEquiv (fun _ : Fin (b + 1) => α))]
  rw [Fintype.sum_prod_type]
  rfl

theorem matrix_sum_ofFn_prod {α : Type*} [Fintype α] {b : Nat}
    (K : Fin b → α → Matrix State State ℝ) :
    (∑ xs : Fin b → α, (List.ofFn (fun i => K i (xs i))).prod) =
      (List.ofFn (fun i => ∑ x, K i x)).prod := by
  induction b with
  | zero => simp
  | succ b ih =>
    rw [matrix_sum_tuple_succ]
    simp only [List.ofFn_succ, List.prod_cons, Fin.cons_zero, Fin.cons_succ]
    simp_rw [← Finset.mul_sum, ih]
    rw [← Finset.sum_mul]

theorem matrix_ofFn_smul_prod {b : Nat} (c : Fin b → ℝ) (K : Fin b → Matrix State State ℝ) :
    (List.ofFn (fun i => c i • K i)).prod = (∏ i, c i) • (List.ofFn K).prod := by
  induction b with
  | zero => simp
  | succ b ih =>
    simp only [List.ofFn_succ, List.prod_cons, Fin.prod_univ_succ]
    rw [ih, smul_mul_smul_comm]

theorem fugacity_product_expansion {R Q b : Nat} (θ : ℝ)
    (marks : Fin b → (Fin Q ↪ Fin (128 * R))) (u : Fin Q → ℝ) :
    (List.ofFn (fun j => fugacityKernel θ (marks j) u)).prod =
      ∑ A : Fin b → Finset (Fin Q), (∏ j, markCoefficient u (A j)) •
        (List.ofFn (fun j => Matrix.of (shuffledRegionKernel θ ((A j).map (marks j))))).prod := by
  have he (j : Fin b) : fugacityKernel θ (marks j) u = ∑ A : Finset (Fin Q), markCoefficient u A • Matrix.of (shuffledRegionKernel θ (A.map (marks j))) := by
    ext q r; simp only [fugacityKernel, Matrix.sum_apply, Matrix.smul_apply, smul_eq_mul]; rfl
  simp only [he]
  rw [← matrix_sum_ofFn_prod (fun j A => markCoefficient u A • Matrix.of (shuffledRegionKernel θ (A.map (marks j))))]
  · apply sum_congr rfl
    intro A _
    exact matrix_ofFn_smul_prod _ _

theorem fugacity_regionStream_moment {R Q b : Nat} (θ : ℝ)
    (marks : Fin b → (Fin Q ↪ Fin (128 * R))) (u : Fin Q → ℝ) (q : State) :
    (∑ A : Fin b → Finset (Fin Q), (∏ j, markCoefficient u (A j)) *
      (Spin.piPMF (fun j => shuffleLaw ((A j).map (marks j)))).expect
        (fun regions => inputMoment (Real.exp (-(θ / (128 * R)))) (regionStream regions).get q)) =
      ∑ r, (List.ofFn (fun j => fugacityKernel θ (marks j) u)).prod q r := by
  rw [fugacity_product_expansion]
  simp only [shuffled_regionStream_moment, Matrix.sum_apply, Matrix.smul_apply, smul_eq_mul, Finset.mul_sum]
  exact sum_comm

end Spin.Structured.Placement

