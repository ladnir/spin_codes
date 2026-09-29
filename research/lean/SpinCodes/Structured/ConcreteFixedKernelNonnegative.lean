import SpinCodes.Structured.ConcreteFixedOneMoment

/-! All continuum kernels inherited from the actual encoder are nonnegative. -/
noncomputable section
namespace Spin.Structured.Placement
open Finset ConcreteEncoder MeasureTheory

theorem matrix_nonneg_mul {α : Type*} [Fintype α] (K L : Matrix α α ℝ)
    (hK : ∀ i j, 0 ≤ K i j) (hL : ∀ i j, 0 ≤ L i j) (i j : α) : 0 ≤ (K*L) i j :=
  sum_nonneg (fun s _ => mul_nonneg (hK i s) (hL s j))

theorem timeEmpty_nonneg (γ t : ℝ) (i j : Fin 2) : 0 ≤ timeEmpty γ t i j := by
  fin_cases i <;> fin_cases j <;> simp [timeEmpty, (Real.exp_pos _).le]

theorem timeProduct_nonneg {a : Nat} (γ : ℝ) (gaps : Fin a → ℝ) (last : ℝ) (i j : Fin 2) :
    0 ≤ timeProduct γ gaps last i j := by
  induction a generalizing i j with
  | zero => exact timeEmpty_nonneg γ last i j
  | succ a ih =>
    change 0 ≤ ((timeEmpty γ (gaps 0) * impulseMatrix) * timeProduct γ (Fin.tail gaps) last) i j
    exact matrix_nonneg_mul _ _
      (matrix_nonneg_mul _ _ (timeEmpty_nonneg _ _) impulseMatrix_substochastic.nonneg)
      (fun k l => ih (Fin.tail gaps) k l) i j

theorem continuumRegionKernel_nonneg (θ : ℝ) (a : Nat) (i j : Fin 2) :
    0 ≤ continuumRegionKernel θ a i j := by
  unfold continuumRegionKernel
  apply mul_nonneg (Nat.cast_nonneg _)
  apply setIntegral_nonneg (orderedSiteDomain_measurable _)
  intro x hx
  exact timeProduct_nonneg _ _ _ i j

theorem fugacityContinuum_nonneg {Q : Nat} (θ : ℝ) (u : Fin Q → ℝ) (hu : ∀ i, 0 ≤ u i)
    (i j : Fin 2) : 0 ≤ fugacityContinuum θ u i j := by
  simp only [fugacityContinuum, Matrix.sum_apply, Matrix.smul_apply, smul_eq_mul]
  exact sum_nonneg (fun A _ => mul_nonneg (markCoefficient_nonneg hu A) (continuumRegionKernel_nonneg θ A.card i j))

end Spin.Structured.Placement

