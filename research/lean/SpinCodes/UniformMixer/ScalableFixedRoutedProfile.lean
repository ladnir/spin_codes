import SpinCodes.UniformMixer.ScalableFixedTwoNumeric
import SpinCodes.UniformMixer.ScalableRouteDomination
import SpinCodes.UniformMixer.ScalableFullMoment
import SpinCodes.Structured.ConcreteFixedProfileRouted

/-! Fixed-row profile bounds after the actual two-stage Scalable route. -/

noncomputable section

namespace Spin.UniformMixer

open Spin.Structured Spin.Structured.ConcreteRoute
  Spin.Structured.Placement Filter

theorem scalable_fixed_routed_profile_margin {Q : ℕ}
    {θ v c δ : ℝ} (hθ : 0 ≤ θ) (hv : 0 < v)
    (hc : 0 ≤ c) (hδ : 0 < δ)
    (u : Fin Q → ℝ) (hu : ∀ i, 0 < u i)
    (hK : scalableTwoRowNormLe v c (scalableFugacityTwo θ u)) :
    ∀ᶠ m : ℕ in atTop,
      ∀ e : Fin Q ↪ Fin (Lsched m),
        ∀ rows : Fin Q → Finset (Fin (bsched m)),
          (law (embedRows e rows)).expect
            (scalableRoutedMoment m
              (Real.exp (-(θ / (Lsched m : ℝ))))) ≤
          (((c + δ) ^ bsched m / min 1 v) /
            (profileChoices rows *
              profileCoefficient u (fun i => (rows i).card))) := by
  filter_upwards [scalable_fixed_profile_two_state_margin
    hθ hv hc hδ u hu hK] with m hm
  intro e rows
  rw [route_embed_expect]
  change (rowLaw rows).expect (fun shuffled =>
    (regionKernel (fun j => ((transpose shuffled) j).map e)).expect
      (fun regions => scalableRoutedMoment m
        (Real.exp (-(θ / (Lsched m : ℝ)))) regions)) ≤ _
  rw [rowLaw_pattern_expect rows (fun A =>
    (regionKernel (fun j => (A j).map e)).expect
      (fun regions => scalableRoutedMoment m
        (Real.exp (-(θ / (Lsched m : ℝ)))) regions))]
  have hh := div_le_div_of_nonneg_right
    (hm (fun _ => e) (fun i => (rows i).card))
    (profileChoices_pos rows).le
  simp only [scalableRoutedMoment] at hh ⊢
  convert hh using 1
  field_simp [ne_of_gt (profileChoices_pos rows),
    ne_of_gt (profileCoefficient_pos hu (fun i => (rows i).card))]

theorem scalable_fixed_routed_probability {Q : ℕ}
    {θ v c δ : ℝ} (hθ : 0 ≤ θ) (hv : 0 < v)
    (hc : 0 ≤ c) (hδ : 0 < δ)
    (u : Fin Q → ℝ) (hu : ∀ i, 0 < u i)
    (hK : scalableTwoRowNormLe v c (scalableFugacityTwo θ u)) :
    ∀ᶠ m : ℕ in atTop,
      ∀ e : Fin Q ↪ Fin (Lsched m),
        ∀ rows : Fin Q → Finset (Fin (bsched m)),
          ∀ d : ℕ,
            ((law (embedRows e rows)).prod
              (scalableStateMapLaw m)).prob
              (fun ω =>
                sampledOutputWeight (scalableFlatInput m ω.1)
                  (scalableFlatMaps m ω.2) 0 ≤ d) ≤
            ((((c + δ) ^ bsched m / min 1 v) /
                (profileChoices rows *
                  profileCoefficient u (fun i => (rows i).card))) /
              Real.exp (-(θ / (Lsched m : ℝ))) ^ d) := by
  filter_upwards [scalable_fixed_routed_profile_margin
    hθ hv hc hδ u hu hK] with m hm
  intro e rows d
  apply FinPMF.prob_weight_le_of_moment _ _ d
    (Real.exp_pos _)
    (Real.exp_le_one_iff.mpr (by
      apply neg_nonpos.mpr
      have hL : (0 : ℝ) < Lsched m := by
        exact_mod_cast Lsched_pos m
      exact div_nonneg hθ hL.le))
  rw [FinPMF.expect_prod]
  simp_rw [scalable_fixed_route_moment]
  exact hm e rows

end Spin.UniformMixer
