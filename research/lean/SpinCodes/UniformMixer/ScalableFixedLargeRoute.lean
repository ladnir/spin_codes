import SpinCodes.UniformMixer.ScalableFixedLargeGate
import SpinCodes.UniformMixer.ScalableFixedRoutedProfile

/-! Uniform fixed-Q profile bound under the actual Scalable route and maps. -/

noncomputable section

namespace Spin.UniformMixer

open Spin.Structured Spin.Structured.ConcreteRoute
  Spin.Structured.ConcreteNativeFamily Spin.Structured.Placement
  Filter

theorem scalable_fixed_large_routed {Q : ℕ} (hQ : 3 ≤ Q) :
    ∀ᶠ m : ℕ in atTop,
      ∀ u : FixedFugacityChoice Q,
        ∀ e : Fin Q ↪ Fin (Lsched m),
          ∀ rows : Fin Q → Finset (Fin (bsched m)),
            ∀ d : ℕ,
              ((law (embedRows e rows)).prod
                (scalableStateMapLaw m)).prob
                (fun ω => sampledOutputWeight
                  (scalableFlatInput m ω.1)
                  (scalableFlatMaps m ω.2) 0 ≤ d) ≤
              ((((fixedNormCost u + fixedNormSlack u) ^
                    bsched m / min 1 (3 / 1600 : ℝ)) /
                  (profileChoices rows *
                    profileCoefficient (fixedFugacity u)
                      (fun i => (rows i).card))) /
                Real.exp (-(scalableLargeTheta Q /
                  (Lsched m : ℝ))) ^ d) := by
  classical
  apply Filter.eventually_all.mpr
  intro u
  have h := scalable_fixed_routed_probability
    (Q := Q) (θ := scalableLargeTheta Q)
    (v := 3 / 1600) (c := fixedNormCost u)
    (δ := fixedNormSlack u)
    (by unfold scalableLargeTheta scalableThetaOne; positivity)
    (by norm_num)
    (fixedNormCost_pos u).le
    (fixedNormSlack_pos (by omega) u)
    (fixedFugacity u) (fixedFugacity_pos u)
    (scalableLargeNorm_proved hQ u)
  filter_upwards [h] with m hm
  intro e rows d
  exact hm e rows d

end Spin.UniformMixer
