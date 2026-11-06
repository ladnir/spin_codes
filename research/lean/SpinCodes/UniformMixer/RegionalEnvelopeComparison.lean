import SpinCodes.UniformMixer.RegionalBoundaryError

/-!
# Pointwise comparison of the production regional envelope

For a good routed support, the actual coarse envelope differs from the
continuum integrand by a uniform `O(a/R)` additive term, multiplied by
the already proved finite gap-comparison factor.
-/

noncomputable section

namespace Spin.UniformMixer

open Spin.Structured.Placement

theorem regionGamma_nonneg {θ : ℝ} (hθ : 0 ≤ θ) :
    0 ≤ regionGamma θ := by
  have hμ : 0 ≤ liveWeightAverage coordMap := by
    rw [production_liveWeightAverage]
    norm_num
  unfold regionGamma
  positivity

theorem good_coarseRegionalEnvelope_continuum {R a : ℕ} (hR : 0 < R)
    (B : BlockSubset (64 * R) a) (hgood : ¬ badRegionalSites B.val)
    (θ : ℝ) (hθ : 0 ≤ θ) (q r : State22) :
    let L : ℝ := 64 * R
    let C := ((64 * θ + 64 * θ ^ 2) / L) / Real.exp (-θ)
    coarseRegionalEnvelope θ B.val q r ≤
      (1 + C) ^ (a + 1) *
        ((liftM * continuumSiteProduct (regionGamma θ)
          (normalizedSites B) * projectM) q r +
          regionGamma θ * (3 * (a + 1 : ℕ) / (R : ℝ))) := by
  have hfun :
      (fun g : ℕ => Real.exp
        (-((θ / (64 * (R : ℝ))) * liveWeightAverage coordMap)) ^ g) =
      (fun g : ℕ => Real.exp
        (-(regionGamma θ) * ((g : ℝ) / R))) := by
    funext g
    exact production_coarse_gap_factor θ R g
  have hpoint := good_support_boundary_error hR B hgood
    (regionGamma_nonneg hθ) q r
  rw [← hfun] at hpoint
  have hC : 0 ≤ ((64 * θ + 64 * θ ^ 2) /
      (64 * (R : ℝ))) / Real.exp (-θ) := by
    have hRp : (0 : ℝ) < R := by exact_mod_cast hR
    positivity
  have hpow : 0 ≤ (1 + ((64 * θ + 64 * θ ^ 2) /
      (64 * (R : ℝ))) / Real.exp (-θ)) ^ (a + 1) := by
    positivity
  have hmul := mul_le_mul_of_nonneg_left hpoint hpow
  simpa only [coarseRegionalEnvelope, B.property,
    Nat.cast_mul, Nat.cast_ofNat] using hmul

end Spin.UniformMixer
