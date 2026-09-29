import SpinCodes.Structured.ConcreteFixedInsertionPathBound

/-! The actual continuum fugacity norm follows from the scalar spacing estimate. -/
noncomputable section
namespace Spin.Structured.Placement
open Finset ConcreteEncoder Matrix MeasureTheory
attribute [local instance] Classical.propDecidable

theorem integrated_row_action {Q : Nat} (γ : ℝ) (P : Matrix (Fin 2) (Fin 2) ℝ)
    (u : Fin Q → ℝ) (w : Fin 2 → ℝ) (i : Fin 2) :
    (∑ j, (∫ x in orderedSiteDomain Q, fullSiteProduct γ P x u i j) * w j) =
      ∫ x in orderedSiteDomain Q, (fullSiteProduct γ P x u *ᵥ w) i := by
  symm
  change (∫ x in orderedSiteDomain Q, ∑ j, fullSiteProduct γ P x u i j * w j) = _
  have hint (j : Fin 2) : IntegrableOn
      (fun x : Fin Q → ℝ => fullSiteProduct γ P x u i j * w j) (orderedSiteDomain Q) := by
    apply continuous_integrable_ordered
    exact ((continuous_apply j).comp ((continuous_apply i).comp (continuous_fullSiteProduct γ P u))).mul continuous_const
  rw [integral_finsetSum univ (fun j hj => hint j)]
  simp only [integral_mul_const]

theorem enlargedInsertionKernel_action {Q : Nat} (θ : ℝ) (u : Fin Q → ℝ) (w : Fin 2 → ℝ) (i : Fin 2) :
    (enlargedInsertionKernel θ u *ᵥ w) i =
      ∑ π : Equiv.Perm (Fin Q), ∫ x in orderedSiteDomain Q,
        (fullSiteProduct (θ * epochMean / 128) (Spin.Pplus (1/524287)) x (fun k => u (π.symm k)) *ᵥ w) i := by
  simp only [Matrix.mulVec, dotProduct, enlargedInsertionKernel, Finset.sum_mul]
  rw [sum_comm]
  apply sum_congr rfl
  intro π hπ
  exact integrated_row_action _ _ _ _ _

theorem enlargedInsertionKernel_nonneg {Q : Nat} (θ : ℝ) (u : Fin Q → ℝ) (hu : ∀ k, 0 ≤ u k)
    (i j : Fin 2) : 0 ≤ enlargedInsertionKernel θ u i j := by
  apply sum_nonneg
  intro π hπ
  apply setIntegral_nonneg (orderedSiteDomain_measurable Q)
  intro x hx
  apply fullGapProduct_nonneg _ _ _ _ _ (fun k => hu (π.symm k))
  intro k l
  fin_cases k <;> fin_cases l <;> norm_num [Spin.Pplus]

theorem enlargedInsertionKernel_action_bound {Q : Nat} {θ σ v C : ℝ}
    (hσ0 : 0 ≤ σ) (hσ1 : σ ≤ 1) (hv : 0 < v) (hC : 0 ≤ C)
    (hscalar : ∀ S : Finset (Fin (Q+1)), (Q.factorial : ℝ) *
      (∫ x in orderedSiteDomain Q, Real.exp (-(θ * epochMean / 128) * liveDuration S x)) ≤ C * σ^S.card)
    (u : Fin Q → ℝ) (hu : ∀ k, 0 ≤ u k) (i : Fin 2) :
    (enlargedInsertionKernel θ u *ᵥ twoWeight v) i ≤
      C * (∏ k, Spin.mv σ v (1/524287) (u k)) * twoWeight v i := by
  rw [enlargedInsertionKernel_action]
  have hf : (0:ℝ) < Q.factorial := by positivity
  have hp (π : Equiv.Perm (Fin Q)) :
      (∏ k, Spin.mv σ v (1/524287) (u (π.symm k))) = ∏ k, Spin.mv σ v (1/524287) (u k) :=
    Equiv.prod_comp π.symm (fun k => Spin.mv σ v (1/524287) (u k))
  have hπ (π : Equiv.Perm (Fin Q)) :
      (∫ x in orderedSiteDomain Q,
        (fullSiteProduct (θ * epochMean / 128) (Spin.Pplus (1/524287)) x (fun k => u (π.symm k)) *ᵥ twoWeight v) i) ≤
        (C * (∏ k, Spin.mv σ v (1/524287) (u k)) * twoWeight v i) / (Q.factorial : ℝ) := by
    apply (le_div_iff₀ hf).mpr
    rw [mul_comm]
    have h := fullSiteProduct_integrated_path_bound hσ0 hσ1 hv (by norm_num : (0:ℝ) ≤ 1/524287)
      hC hscalar (fun k => u (π.symm k)) (fun k => hu (π.symm k)) i
    simpa only [hp] using h
  calc
    _ ≤ ∑ _π : Equiv.Perm (Fin Q),
        (C * (∏ k, Spin.mv σ v (1/524287) (u k)) * twoWeight v i) / (Q.factorial : ℝ) := sum_le_sum (fun π h => hπ π)
    _ = _ := by
      simp only [sum_const, card_univ, Fintype.card_perm, Fintype.card_fin, nsmul_eq_mul]
      field_simp

theorem fugacityContinuum_rowNorm_of_spacing {Q : Nat} {θ σ v C : ℝ}
    (hσ0 : 0 ≤ σ) (hσ1 : σ ≤ 1) (hv : 0 < v) (hC : 0 ≤ C)
    (hscalar : ∀ S : Finset (Fin (Q+1)), (Q.factorial : ℝ) *
      (∫ x in orderedSiteDomain Q, Real.exp (-(θ * epochMean / 128) * liveDuration S x)) ≤ C * σ^S.card)
    (u : Fin Q → ℝ) (hu : ∀ k, 0 ≤ u k) :
    Spin.RowNormLe v (C * ∏ k, Spin.mv σ v (1/524287) (u k)) (fugacityContinuum θ u) := by
  apply fugacity_norm_of_enlarged hv.le θ u hu
  have h0 := enlargedInsertionKernel_action_bound hσ0 hσ1 hv hC hscalar u hu 0
  have h1 := enlargedInsertionKernel_action_bound hσ0 hσ1 hv hC hscalar u hu 1
  unfold Spin.RowNormLe
  simp only [abs_of_nonneg (enlargedInsertionKernel_nonneg θ u hu _ _)]
  constructor
  · simpa [Matrix.mulVec, dotProduct, Fin.sum_univ_two, twoWeight] using h0
  · simpa [Matrix.mulVec, dotProduct, Fin.sum_univ_two, twoWeight] using h1

end Spin.Structured.Placement

