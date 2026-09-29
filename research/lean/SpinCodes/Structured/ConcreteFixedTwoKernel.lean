import SpinCodes.Structured.ConcreteFixedTriangle

/-! Two-site continuum kernel comparison with the explicit enlarged kernel. -/
noncomputable section
namespace Spin.Structured.Placement
open Finset ConcreteEncoder MeasureTheory
attribute [local instance] Classical.propDecidable

theorem siteProduct_two (γ : ℝ) (x : Fin 2 → ℝ) :
    siteProduct γ x = timeEmpty γ (x 0) * impulseMatrix * timeEmpty γ (x 1-x 0) *
      impulseMatrix * timeEmpty γ (1-x 1) := by
  have h0 : positionGaps x (0 : Fin 3) = x 0 := by norm_num [positionGaps, Fin.snoc, Fin.cons]
  have h1 : positionGaps x (1 : Fin 3) = x 1-x 0 := by
    norm_num [positionGaps, Fin.snoc, Fin.cons]
    rfl
  have h2 : positionGaps x (Fin.last 2) = 1-x 1 := by
    norm_num [positionGaps, Fin.snoc, Fin.cons]
    rfl
  change timeEmpty γ (positionGaps x (0:Fin 3)) * impulseMatrix *
    (timeEmpty γ (positionGaps x (1:Fin 3)) * impulseMatrix * timeEmpty γ (positionGaps x (Fin.last 2))) = _
  rw [h0,h1,h2]
  simp only [Matrix.mul_assoc]

def plusSiteTwo (γ : ℝ) (x : Fin 2 → ℝ) : Matrix (Fin 2) (Fin 2) ℝ :=
  timeEmpty γ (x 0) * Spin.Pplus (1/524287) * timeEmpty γ (x 1-x 0) *
    Spin.Pplus (1/524287) * timeEmpty γ (1-x 1)

theorem continuous_plusSiteTwo (γ : ℝ) : Continuous (plusSiteTwo γ) := by
  exact (((((continuous_timeEmpty γ).comp (continuous_apply 0)).mul continuous_const).mul
    ((continuous_timeEmpty γ).comp ((continuous_apply 1).sub (continuous_apply 0)))).mul continuous_const).mul
    ((continuous_timeEmpty γ).comp (continuous_const.sub (continuous_apply 1)))

theorem plusSiteTwo_gaps (γ u v : ℝ) :
    plusSiteTwo γ ![u,u+v] = Spin.Dg γ u * Spin.Pplus (1/524287) * Spin.Dg γ v *
      Spin.Pplus (1/524287) * Spin.Dg γ (1-u-v) := by
  simp only [plusSiteTwo, Matrix.cons_val_zero, Matrix.cons_val_one, add_sub_cancel_left]
  rw [show (1:ℝ)-(u+v)=1-u-v by ring]
  rfl

theorem siteProduct_two_le (γ : ℝ) (x : Fin 2 → ℝ) (i j : Fin 2) :
    siteProduct γ x i j ≤ plusSiteTwo γ x i j := by
  rw [siteProduct_two]
  have h0 := (Real.exp_pos (-γ * x 0)).le
  have h1 := (Real.exp_pos (-γ * (x 1-x 0))).le
  have h2 := (Real.exp_pos (-γ * (1-x 1))).le
  have h01 := mul_nonneg h0 h1
  have h12 := mul_nonneg h1 h2
  have h02 := mul_nonneg h0 h2
  have h012 := mul_nonneg h01 h2
  simp only [neg_mul] at h0 h1 h2 h01 h12 h02 h012
  fin_cases i <;> fin_cases j <;>
    simp [plusSiteTwo, timeEmpty, Spin.Pplus, impulseMatrix, Matrix.mul_apply, Fin.sum_univ_two] <;>
    nlinarith

def enlargedKernelTwo (γ : ℝ) : Matrix (Fin 2) (Fin 2) ℝ :=
  !![(1/524287) * (2 * (γ-1+Real.exp (-γ)) / γ^2), 2*(1-(1+γ)*Real.exp (-γ))/γ^2;
    (1/524287) * (2*(1-(1+γ)*Real.exp (-γ))/γ^2),
    (1/524287) * (2*(1-(1+γ)*Real.exp (-γ))/γ^2) + Real.exp (-γ)]

theorem plusSiteTwo_integral {γ : ℝ} (hγ : γ ≠ 0) (i j : Fin 2) :
    2 * (∫ x in orderedSiteDomain 2, plusSiteTwo γ x i j) = enlargedKernelTwo γ i j := by
  rw [orderedSite_integral_two_gaps (fun x => plusSiteTwo γ x i j)
    ((continuous_apply j).comp ((continuous_apply i).comp (continuous_plusSiteTwo γ)))]
  simp only [plusSiteTwo_gaps]
  fin_cases i <;> fin_cases j
  · exact Spin.K2_entry_00 _ hγ
  · exact Spin.K2_entry_01 _ hγ
  · exact Spin.K2_entry_10 _ hγ
  · exact Spin.K2_entry_11 _ hγ

theorem continuumRegionKernel_two_le {θ : ℝ} (hγ : θ * epochMean / 128 ≠ 0) (i j : Fin 2) :
    continuumRegionKernel θ 2 i j ≤ enlargedKernelTwo (θ * epochMean / 128) i j := by
  rw [← plusSiteTwo_integral hγ i j]
  unfold continuumRegionKernel
  norm_num only [Nat.factorial_succ, Nat.factorial_zero, Nat.reduceAdd, Nat.reduceMul, Nat.cast_ofNat]
  apply mul_le_mul_of_nonneg_left _ (by norm_num)
  have hi (F : (Fin 2 → ℝ) → ℝ) (hF : Continuous F) : IntegrableOn F (orderedSiteDomain 2) := by
    apply hF.integrableOn_Icc.mono_set
    intro x hx
    exact ⟨fun i => (hx.1 i).1, fun i => (hx.1 i).2.le⟩
  apply setIntegral_mono_on
    (hi _ ((continuous_apply j).comp ((continuous_apply i).comp (continuous_siteProduct _ 2))))
    (hi _ ((continuous_apply j).comp ((continuous_apply i).comp (continuous_plusSiteTwo _))))
    (orderedSiteDomain_measurable 2)
  intro x hx
  exact siteProduct_two_le _ x i j

end Spin.Structured.Placement

