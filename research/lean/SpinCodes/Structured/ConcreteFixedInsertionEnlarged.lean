import SpinCodes.Structured.ConcreteFixedInsertionLimit

/-! Enlarging P0 to Pplus gives a genuine upper bound for the actual continuum kernel. -/
noncomputable section
namespace Spin.Structured.Placement
open Finset ConcreteEncoder MeasureTheory
attribute [local instance] Classical.propDecidable

theorem matrix_mul_mono {κ : Type*} [Fintype κ] {K K' L L' : Matrix κ κ ℝ}
    (hK : ∀ i j, K i j ≤ K' i j) (hL : ∀ i j, L i j ≤ L' i j)
    (hL0 : ∀ i j, 0 ≤ L i j) (hK0 : ∀ i j, 0 ≤ K' i j) (i j : κ) :
    (K*L) i j ≤ (K'*L') i j := by
  apply sum_le_sum
  intro s hs
  exact mul_le_mul (hK i s) (hL s j) (hL0 s j) (hK0 i s)

theorem matrix_list_prod_nonneg {κ : Type*} [Fintype κ] [DecidableEq κ]
    (Ks : List (Matrix κ κ ℝ)) (hK : ∀ K ∈ Ks, ∀ i j, 0 ≤ K i j) (i j : κ) : 0 ≤ Ks.prod i j := by
  induction Ks generalizing i j with
  | nil => exact FiniteKernel.Substochastic.one.nonneg i j
  | cons K Ks ih =>
    exact matrix_nonneg_mul _ _ (hK K (by simp)) (ih (fun L hL => hK L (by simp [hL]))) i j

theorem matrix_ofFn_prod_mono {κ : Type*} [Fintype κ] [DecidableEq κ] {a : Nat}
    (K L : Fin a → Matrix κ κ ℝ) (hK0 : ∀ k i j, 0 ≤ K k i j) (hL0 : ∀ k i j, 0 ≤ L k i j)
    (hKL : ∀ k i j, K k i j ≤ L k i j) (i j : κ) :
    (List.ofFn K).prod i j ≤ (List.ofFn L).prod i j := by
  induction a generalizing i j with
  | zero => simp
  | succ a ih =>
    simp only [List.ofFn_succ, List.prod_cons]
    apply matrix_mul_mono (hKL 0) (ih (fun k => K k.succ) (fun k => L k.succ)
      (fun k => hK0 k.succ) (fun k => hL0 k.succ) (fun k => hKL k.succ)) _ (hL0 0)
    intro i j
    apply matrix_list_prod_nonneg
    intro M hM
    obtain ⟨k,rfl⟩ := List.mem_ofFn.mp hM
    exact hK0 k.succ

theorem inserted_factor_nonneg (γ t u : ℝ) (P : Matrix (Fin 2) (Fin 2) ℝ)
    (hu : 0 ≤ u) (hP : ∀ i j, 0 ≤ P i j) (i j : Fin 2) :
    0 ≤ (timeEmpty γ t * (1 + u • P)) i j := by
  apply matrix_nonneg_mul _ _ (timeEmpty_nonneg _ _)
  intro k l
  exact add_nonneg (FiniteKernel.Substochastic.one.nonneg k l) (mul_nonneg hu (hP k l))

theorem fullGapProduct_nonneg {Q : Nat} (γ : ℝ) (P : Matrix (Fin 2) (Fin 2) ℝ)
    (gaps u : Fin Q → ℝ) (last : ℝ) (hu : ∀ k, 0 ≤ u k) (hP : ∀ i j, 0 ≤ P i j) (i j : Fin 2) :
    0 ≤ fullGapProduct γ P gaps u last i j := by
  apply matrix_nonneg_mul _ _ _ (timeEmpty_nonneg _ _) i j
  apply matrix_list_prod_nonneg
  intro M hM
  obtain ⟨k,rfl⟩ := List.mem_ofFn.mp hM
  exact inserted_factor_nonneg _ _ _ _ (hu k) hP

theorem fullGapProduct_mono {Q : Nat} (γ : ℝ) (P P' : Matrix (Fin 2) (Fin 2) ℝ)
    (gaps u : Fin Q → ℝ) (last : ℝ) (hu : ∀ k, 0 ≤ u k)
    (hP : ∀ i j, 0 ≤ P i j) (hP' : ∀ i j, 0 ≤ P' i j) (hle : ∀ i j, P i j ≤ P' i j)
    (i j : Fin 2) : fullGapProduct γ P gaps u last i j ≤ fullGapProduct γ P' gaps u last i j := by
  apply matrix_mul_mono _ (fun _ _ => le_rfl) (timeEmpty_nonneg _ _) _ i j
  · apply matrix_ofFn_prod_mono
      _ _ (fun k => inserted_factor_nonneg _ _ _ _ (hu k) hP)
      (fun k => inserted_factor_nonneg _ _ _ _ (hu k) hP')
    intro k i j
    apply matrix_mul_mono (fun _ _ => le_rfl) _ _ (timeEmpty_nonneg _ _) i j
    · intro s t
      exact add_le_add le_rfl (mul_le_mul_of_nonneg_left (hle s t) (hu k))
    · intro s t
      exact add_nonneg (FiniteKernel.Substochastic.one.nonneg s t) (mul_nonneg (hu k) (hP s t))
  · apply matrix_list_prod_nonneg
    intro M hM
    obtain ⟨k,rfl⟩ := List.mem_ofFn.mp hM
    exact inserted_factor_nonneg _ _ _ _ (hu k) hP'

def enlargedInsertionKernel {Q : Nat} (θ : ℝ) (u : Fin Q → ℝ) : Matrix (Fin 2) (Fin 2) ℝ :=
  fun i j => ∑ π : Equiv.Perm (Fin Q), ∫ x in orderedSiteDomain Q,
    fullSiteProduct (θ * epochMean / 128) (Spin.Pplus (1/524287)) x (fun k => u (π.symm k)) i j

theorem fugacityContinuum_le_enlarged {Q : Nat} (θ : ℝ) (u : Fin Q → ℝ) (hu : ∀ k, 0 ≤ u k)
    (i j : Fin 2) : fugacityContinuum θ u i j ≤ enlargedInsertionKernel θ u i j := by
  rw [fugacityContinuum_eq_insertion]
  apply sum_le_sum
  intro π hπ
  have hi (P : Matrix (Fin 2) (Fin 2) ℝ) : IntegrableOn
      (fun x => fullSiteProduct (θ * epochMean / 128) P x (fun k => u (π.symm k)) i j) (orderedSiteDomain Q) := by
    have hc := (continuous_apply j).comp ((continuous_apply i).comp (continuous_fullSiteProduct (θ*epochMean/128) P (fun k => u (π.symm k))))
    apply hc.integrableOn_Icc.mono_set
    intro x hx
    exact ⟨fun i => (hx.1 i).1,fun i => (hx.1 i).2.le⟩
  apply setIntegral_mono_on (hi _) (hi _) (orderedSiteDomain_measurable Q)
  intro x hx
  apply fullGapProduct_mono _ _ _ _ _ _ (fun k => hu (π.symm k)) impulseMatrix_substochastic.nonneg
  · intro i j; fin_cases i <;> fin_cases j <;> norm_num [Spin.Pplus]
  · intro i j; fin_cases i <;> fin_cases j <;> norm_num [Spin.Pplus,impulseMatrix]

theorem fugacity_norm_of_enlarged {Q : Nat} {v c : ℝ} (hv : 0 ≤ v) (θ : ℝ)
    (u : Fin Q → ℝ) (hu : ∀ k, 0 ≤ u k) (h : Spin.RowNormLe v c (enlargedInsertionKernel θ u)) :
    Spin.RowNormLe v c (fugacityContinuum θ u) := by
  unfold Spin.RowNormLe
  simp only [abs_of_nonneg (fugacityContinuum_nonneg θ u hu _ _)]
  have hh (i j : Fin 2) := (fugacityContinuum_le_enlarged θ u hu i j).trans (le_abs_self _)
  constructor
  · exact (add_le_add (hh 0 0) (mul_le_mul_of_nonneg_right (hh 0 1) hv)).trans h.1
  · exact (add_le_add (hh 1 0) (mul_le_mul_of_nonneg_right (hh 1 1) hv)).trans h.2

end Spin.Structured.Placement

