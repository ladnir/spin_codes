import SpinCodes.Structured.ConcreteFixedInsertionFinite

/-! Full-site insertion and the actual continuum fugacity kernel agree exactly. -/
noncomputable section
namespace Spin.Structured.Placement
open Finset ConcreteEncoder MeasureTheory Filter
open scoped Topology
attribute [local instance] Classical.propDecidable

theorem uniform_orderedSite_sum_tendsto {a : Nat} (F : (Fin a → ℝ) → ℝ) (hF : Continuous F) :
    Tendsto (fun R : Nat => (∑ B : BlockSubset R a, F (normalizedSites B)) / (R.choose a : ℝ))
      atTop (𝓝 ((a.factorial : ℝ) * ∫ x in orderedSiteDomain a, F x)) := by
  have h := (pow_div_choose_tendsto a).mul (orderedSite_sum_tendsto F hF)
  apply h.congr'
  filter_upwards [eventually_gt_atTop 0] with R hR
  have hRp : (R:ℝ) ≠ 0 := by exact_mod_cast Nat.ne_of_gt hR
  field_simp

theorem uniform_siteProduct_sum_tendsto (θ : ℝ) (a : Nat) (i j : Fin 2) :
    Tendsto (fun R : Nat => (∑ B : BlockSubset R a,
      siteProduct (θ * epochMean / 128) (normalizedSites B) i j) / (R.choose a : ℝ))
      atTop (𝓝 (continuumRegionKernel θ a i j)) := by
  exact uniform_orderedSite_sum_tendsto (fun x => siteProduct (θ * epochMean / 128) x i j)
    ((continuous_apply j).comp ((continuous_apply i).comp (continuous_siteProduct _ a)))

theorem finiteInsertionAverage_tendsto_fugacity {Q : Nat} (θ : ℝ) (u : Fin Q → ℝ) (i j : Fin 2) :
    Tendsto (fun R : Nat => finiteInsertionAverage R θ u i j) atTop (𝓝 (fugacityContinuum θ u i j)) := by
  have h := tendsto_finset_sum univ (fun A (_ : A ∈ univ) =>
    (uniform_siteProduct_sum_tendsto θ A.card i j).const_mul (markCoefficient u A))
  have hh : Tendsto (fun R : Nat => ∑ A : Finset (Fin Q), markCoefficient u A *
      ((∑ C : BlockSubset R A.card, siteProduct (θ * epochMean / 128) (normalizedSites C) i j) /
        (R.choose A.card : ℝ))) atTop (𝓝 (fugacityContinuum θ u i j)) := by
    simpa only [fugacityContinuum, Matrix.sum_apply, Matrix.smul_apply, smul_eq_mul] using h
  apply hh.congr'
  filter_upwards [eventually_ge_atTop Q] with R hR
  exact (finiteInsertionAverage_eq hR θ u i j).symm

def averagedInsertionKernel {Q : Nat} (θ : ℝ) (u : Fin Q → ℝ) : Matrix (Fin 2) (Fin 2) ℝ :=
  fun i j => ∑ π : Equiv.Perm (Fin Q), ∫ x in orderedSiteDomain Q,
    fullSiteProduct (θ * epochMean / 128) impulseMatrix x (fun k => u (π.symm k)) i j

theorem finiteInsertionAverage_tendsto_insertion {Q : Nat} (θ : ℝ) (u : Fin Q → ℝ) (i j : Fin 2) :
    Tendsto (fun R : Nat => finiteInsertionAverage R θ u i j) atTop (𝓝 (averagedInsertionKernel θ u i j)) := by
  have hπ (π : Equiv.Perm (Fin Q)) :
      Tendsto (fun R : Nat => ((∑ B : BlockSubset R Q,
        fullSiteProduct (θ * epochMean / 128) impulseMatrix (normalizedSites B) (fun k => u (π.symm k)) i j) /
          (R.choose Q : ℝ)) / (Q.factorial : ℝ)) atTop
        (𝓝 (∫ x in orderedSiteDomain Q,
          fullSiteProduct (θ * epochMean / 128) impulseMatrix x (fun k => u (π.symm k)) i j)) := by
    have h := (uniform_orderedSite_sum_tendsto
      (fun x => fullSiteProduct (θ * epochMean / 128) impulseMatrix x (fun k => u (π.symm k)) i j)
      ((continuous_apply j).comp ((continuous_apply i).comp (continuous_fullSiteProduct _ _ _)))).div_const
        (Q.factorial : ℝ)
    simpa only [mul_div_cancel_left₀ _ (by positivity : (Q.factorial : ℝ) ≠ 0)] using h
  have hh := tendsto_finset_sum univ (fun π (_ : π ∈ univ) => hπ π)
  convert hh using 1
  · funext R
    unfold finiteInsertionAverage
    rw [sum_comm, sum_div]
    apply sum_congr rfl
    intro π hπ
    exact (div_div _ _ _).symm
  · rfl

theorem fugacityContinuum_eq_insertion {Q : Nat} (θ : ℝ) (u : Fin Q → ℝ) :
    fugacityContinuum θ u = averagedInsertionKernel θ u := by
  ext i j
  exact tendsto_nhds_unique (finiteInsertionAverage_tendsto_fugacity θ u i j)
    (finiteInsertionAverage_tendsto_insertion θ u i j)

end Spin.Structured.Placement
