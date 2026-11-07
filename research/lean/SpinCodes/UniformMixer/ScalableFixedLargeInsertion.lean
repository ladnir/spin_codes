import SpinCodes.UniformMixer.ScalableFixedLargeSite

/-! The scalable fugacity continuum equals a full-site insertion average. -/

noncomputable section

namespace Spin.UniformMixer

open Spin.Structured.Placement Finset Filter MeasureTheory
open scoped Topology

attribute [local instance] Classical.propDecidable

theorem continuous_scalableSiteProduct (γ : ℝ) (a : ℕ) :
    Continuous (fun x : Fin a → ℝ => scalableSiteProduct γ x) := by
  apply continuous_pi
  intro i
  apply continuous_pi
  intro j
  have he : (fun x : Fin a → ℝ => scalableSiteProduct γ x i j) =
      (fun x => continuumSiteProduct γ x (finTwoEquiv i) (finTwoEquiv j)) := by
    funext x
    rw [← scalableSiteProduct_eq, toFinTwo_apply]
  rw [he]
  exact (continuous_apply (finTwoEquiv j)).comp
    ((continuous_apply (finTwoEquiv i)).comp
      (continuous_continuumSiteProduct γ a))

theorem uniform_scalableSiteProduct_sum_tendsto
    (γ : ℝ) (a : ℕ) (i j : Fin 2) :
    Tendsto (fun R : ℕ =>
      (∑ B : BlockSubset R a,
        scalableSiteProduct γ (normalizedSites B) i j) /
          (R.choose a : ℝ)) atTop
      (𝓝 ((a.factorial : ℝ) *
        ∫ x in orderedSiteDomain a,
          scalableSiteProduct γ x i j)) :=
  uniform_orderedSite_sum_tendsto _
    ((continuous_apply j).comp
      ((continuous_apply i).comp
        (continuous_scalableSiteProduct γ a)))

theorem scalableFiniteInsertionAverage_tendsto_fugacity
    {Q : ℕ} (θ : ℝ) (u : Fin Q → ℝ)
    (i j : Fin 2) :
    Tendsto (fun R : ℕ =>
      scalableFiniteInsertionAverage R (regionGamma θ) u i j)
        atTop (𝓝 (toFinTwo (scalableFugacityTwo θ u) i j)) := by
  have h := tendsto_finsetSum univ (fun A (_ : A ∈ univ) =>
    (uniform_scalableSiteProduct_sum_tendsto
      (regionGamma θ) A.card i j).const_mul (markCoefficient u A))
  have hh : Tendsto (fun R : ℕ =>
      ∑ A : Finset (Fin Q), markCoefficient u A *
        ((∑ C : BlockSubset R A.card,
          scalableSiteProduct (regionGamma θ)
            (normalizedSites C) i j) /
            (R.choose A.card : ℝ))) atTop
      (𝓝 (toFinTwo (scalableFugacityTwo θ u) i j)) := by
    rw [toFinTwo_apply]
    simpa only [scalableFugacityTwo, Matrix.sum_apply,
      Matrix.smul_apply, smul_eq_mul,
      ← scalableContinuumTwo_site, toFinTwo_apply] using h
  apply hh.congr'
  filter_upwards [eventually_ge_atTop Q] with R hR
  exact (scalableFiniteInsertionAverage_eq hR
    (regionGamma θ) u i j).symm

def scalableAveragedInsertionKernel {Q : ℕ}
    (γ : ℝ) (u : Fin Q → ℝ) : Matrix (Fin 2) (Fin 2) ℝ :=
  fun i j => ∑ π : Equiv.Perm (Fin Q),
    ∫ x in orderedSiteDomain Q,
      fullSiteProduct γ scalableImpulseFin x
        (fun k => u (π.symm k)) i j

theorem scalableFiniteInsertionAverage_tendsto_insertion
    {Q : ℕ} (γ : ℝ) (u : Fin Q → ℝ)
    (i j : Fin 2) :
    Tendsto (fun R : ℕ =>
      scalableFiniteInsertionAverage R γ u i j)
        atTop (𝓝 (scalableAveragedInsertionKernel γ u i j)) := by
  have hπ (π : Equiv.Perm (Fin Q)) :
      Tendsto (fun R : ℕ =>
        ((∑ B : BlockSubset R Q,
          fullSiteProduct γ scalableImpulseFin
            (normalizedSites B)
            (fun k => u (π.symm k)) i j) /
              (R.choose Q : ℝ)) /
                (Q.factorial : ℝ)) atTop
        (𝓝 (∫ x in orderedSiteDomain Q,
          fullSiteProduct γ scalableImpulseFin x
            (fun k => u (π.symm k)) i j)) := by
    have h := (uniform_orderedSite_sum_tendsto
      (fun x => fullSiteProduct γ scalableImpulseFin x
        (fun k => u (π.symm k)) i j)
      ((continuous_apply j).comp
        ((continuous_apply i).comp
          (continuous_fullSiteProduct _ _ _)))).div_const
            (Q.factorial : ℝ)
    simpa only [mul_div_cancel_left₀ _
      (by positivity : (Q.factorial : ℝ) ≠ 0)] using h
  have hh := tendsto_finsetSum univ
    (fun π (_ : π ∈ univ) => hπ π)
  convert hh using 1
  · funext R
    unfold scalableFiniteInsertionAverage
    rw [sum_comm, sum_div]
    apply sum_congr rfl
    intro π hπ
    exact (div_div _ _ _).symm
  · rfl

theorem scalableFugacityTwo_eq_insertion {Q : ℕ}
    (θ : ℝ) (u : Fin Q → ℝ) :
    toFinTwo (scalableFugacityTwo θ u) =
      scalableAveragedInsertionKernel (regionGamma θ) u := by
  ext i j
  exact tendsto_nhds_unique
    (scalableFiniteInsertionAverage_tendsto_fugacity θ u i j)
    (scalableFiniteInsertionAverage_tendsto_insertion
      (regionGamma θ) u i j)

end Spin.UniformMixer
