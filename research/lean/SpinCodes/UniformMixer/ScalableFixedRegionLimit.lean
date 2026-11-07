import SpinCodes.UniformMixer.ScalableCastShuffle

/-! Transfer the checked fixed-occupation regional limit to the actual
Scalable SPIN region length and shuffle law. -/

noncomputable section

namespace Spin.UniformMixer

open Spin.Structured Filter

theorem scalable_fixed_region_continuum_upper (a : ℕ)
    {θ : ℝ} (hθ : 0 ≤ θ) {ε : ℝ} (hε : 0 < ε) :
    ∀ᶠ m : ℕ in atTop,
      ∀ S : Finset (Fin (Lsched m)), S.card = a →
        ∀ q r : State22,
          scalableAveragedRegion m
            (Real.exp (-(θ / (Lsched m : ℝ)))) S q r <
          continuumEndpoint θ a q r + ε := by
  obtain ⟨R₀, hR₀⟩ := Filter.eventually_atTop.1
    (averagedRegionKernel_uniform_continuum_upper a hθ hε)
  apply Filter.eventually_atTop.2
  refine ⟨R₀, ?_⟩
  intro m hm S hS q r
  have hR : R₀ ≤ scalableStepsPerRegion m := by
    unfold scalableStepsPerRegion
    omega
  have hcard : (scalableRegionSupport m S).card = a := by
    simpa [scalableRegionSupport_eq_congr] using hS
  have h := hR₀ (scalableStepsPerRegion m) hR
    (scalableRegionSupport m S) hcard q r
  have hlen : (Lsched m : ℝ) =
      64 * (scalableStepsPerRegion m : ℝ) := by
    exact_mod_cast scalable_region_length m
  rw [scalableAveragedRegion_eq]
  change averagedRegionKernel (scalableRegionSupport m S)
    (Real.exp (-(θ / (Lsched m : ℝ)))) q r <
      continuumEndpoint θ a q r + ε
  simpa only [hlen] using h

end Spin.UniformMixer
