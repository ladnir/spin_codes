import SpinCodes.UniformMixer.ScalableFixedFugacity
import SpinCodes.Structured.ConcreteFixedInsertionExpansion
import SpinCodes.Structured.ConcreteFixedProfile

/-! Coefficient extraction for fixed active rows under the actual Scalable
regional kernel. The numerical continuum margin remains a premise. -/

noncomputable section

namespace Spin.UniformMixer

open Spin.Structured Spin.Structured.ConcreteRoute
  Spin.Structured.Placement Filter

attribute [local instance] Classical.propDecidable

theorem scalable_fugacity_product_expansion {m Q b : ℕ}
    (θ : ℝ) (marks : Fin b → (Fin Q ↪ Fin (Lsched m)))
    (u : Fin Q → ℝ) :
    (List.ofFn (fun j => scalableFugacityKernel θ (marks j) u)).prod =
      ∑ A : Fin b → Finset (Fin Q),
        (∏ j, markCoefficient u (A j)) •
          (List.ofFn (fun j => scalableAveragedRegion m
            (Real.exp (-(θ / (Lsched m : ℝ))))
            ((A j).map (marks j)))).prod := by
  have he (j : Fin b) :
      scalableFugacityKernel θ (marks j) u =
        ∑ A : Finset (Fin Q), markCoefficient u A •
          scalableAveragedRegion m
            (Real.exp (-(θ / (Lsched m : ℝ))))
            (A.map (marks j)) := by
    ext q r
    simp only [scalableFugacityKernel, Matrix.sum_apply,
      Matrix.smul_apply, smul_eq_mul]
  simp only [he]
  rw [← finite_matrix_sum_product (fun j A =>
    markCoefficient u A • scalableAveragedRegion m
      (Real.exp (-(θ / (Lsched m : ℝ)))) (A.map (marks j)))]
  apply Finset.sum_congr rfl
  intro A _
  exact finite_matrix_smul_product _ _

theorem scalable_fugacity_region_moment {m Q : ℕ} (θ : ℝ)
    (marks : Fin (bsched m) → (Fin Q ↪ Fin (Lsched m)))
    (u : Fin Q → ℝ) (q : State22) :
    (∑ A : Fin (bsched m) → Finset (Fin Q),
      (∏ j, markCoefficient u (A j)) *
        (regionKernel (fun j => (A j).map (marks j))).expect
          (fun regions =>
            ∑ r : State22,
              inputListKernel (Real.exp (-(θ / (Lsched m : ℝ))))
                (List.ofFn (scalableFlatInput m regions)) q r)) =
      ∑ r : State22,
        (List.ofFn (fun j => scalableFugacityKernel θ (marks j) u)).prod q r := by
  have hpoint (A : Fin (bsched m) → Finset (Fin Q)) :
      (regionKernel (fun j => (A j).map (marks j))).expect
          (fun regions =>
            ∑ r : State22,
              inputListKernel (Real.exp (-(θ / (Lsched m : ℝ))))
                (List.ofFn (scalableFlatInput m regions)) q r) =
      ∑ r : State22,
        (List.ofFn (fun j => scalableAveragedRegion m
          (Real.exp (-(θ / (Lsched m : ℝ))))
          ((A j).map (marks j)))).prod q r := by
    rw [Spin.FinPMF.expect_sum]
    apply Finset.sum_congr rfl
    intro r _
    exact congrFun (congrFun
      (scalable_region_shuffle_product m
        (Real.exp (-(θ / (Lsched m : ℝ))))
        (fun j => (A j).map (marks j))) q) r
  rw [scalable_fugacity_product_expansion]
  simp only [hpoint, Matrix.sum_apply, Matrix.smul_apply,
    smul_eq_mul, Finset.mul_sum]
  exact Finset.sum_comm

theorem scalable_fixed_profile_margin {Q : ℕ}
    {θ c δ m₀ : ℝ} (hθ : 0 ≤ θ) (hc : 0 ≤ c) (hδ : 0 < δ)
    (u : Fin Q → ℝ) (hu : ∀ i, 0 < u i)
    (weight : State22 → ℝ) (hw : ∀ r, 0 < weight r)
    (hm₀ : 0 < m₀) (hmw : ∀ r, m₀ ≤ weight r)
    (hK : ∀ q : State22,
      (∑ r, scalableFugacityContinuum θ u q r * weight r) ≤
        c * weight q) :
    ∀ᶠ m : ℕ in atTop,
      ∀ marks : Fin (bsched m) → (Fin Q ↪ Fin (Lsched m)),
        ∀ rowWeights : Fin Q → ℕ,
          (∑ A : Fin (bsched m) → Finset (Fin Q),
            if rowCount A = rowWeights then
              (regionKernel (fun j => (A j).map (marks j))).expect
                (fun regions =>
                  ∑ r : State22,
                    inputListKernel (Real.exp (-(θ / (Lsched m : ℝ))))
                      (List.ofFn (scalableFlatInput m regions)) 0 r)
            else 0) ≤
          ((c + δ) ^ bsched m * weight 0 / m₀) /
            profileCoefficient u rowWeights := by
  filter_upwards [scalable_fugacity_product_total hθ hc hδ
    u (fun i => (hu i).le) weight hw hm₀ hmw hK] with m hm
  intro marks rowWeights
  apply (le_div_iff₀ (profileCoefficient_pos hu rowWeights)).mpr
  rw [mul_comm]
  calc
    _ ≤ ∑ A : Fin (bsched m) → Finset (Fin Q),
        (∏ j, markCoefficient u (A j)) *
          (regionKernel (fun j => (A j).map (marks j))).expect
            (fun regions =>
              ∑ r : State22,
                inputListKernel (Real.exp (-(θ / (Lsched m : ℝ))))
                  (List.ofFn (scalableFlatInput m regions)) 0 r) := by
      apply profile_coefficient_bound u (fun i => (hu i).le)
      intro A
      apply Spin.FinPMF.expect_nonneg
      intro regions
      apply Finset.sum_nonneg
      intro r _
      exact (inputListKernel_substochastic _ (Real.exp_pos _).le
        (by
          have hL : (0 : ℝ) < (Lsched m : ℝ) := by
            unfold Lsched
            positivity
          exact Real.exp_le_one_iff.mpr
            (neg_nonpos.mpr (div_nonneg hθ hL.le)))
        (List.ofFn (scalableFlatInput m regions))).1 0 r
    _ = ∑ r : State22,
        (List.ofFn (fun j => scalableFugacityKernel θ (marks j) u)).prod 0 r :=
      scalable_fugacity_region_moment θ marks u 0
    _ ≤ _ := hm (bsched m) marks

end Spin.UniformMixer
