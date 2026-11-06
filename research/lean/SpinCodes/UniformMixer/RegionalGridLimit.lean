import SpinCodes.UniformMixer.RegionalSiteSum
import SpinCodes.Structured.ConcretePlacementLimitRiemann

/-!
# Uniform support grid limit for 64-bit-step regions

The exact support sum from `RegionalSiteSum` is a Riemann sum over sorted
bit positions. This module proves its limit for any continuous integrand.
The separate task is to compare the discrete two-state product to such an
integrand uniformly over supports.
-/

noncomputable section

namespace Spin.UniformMixer

open Spin.Structured.Placement Filter MeasureTheory
open scoped Topology

theorem regional_grid_limit {a : ℕ}
    (F : (Fin a → ℝ) → ℝ) (hF : Continuous F) :
    Tendsto (fun R : ℕ =>
      (∑ B : BlockSubset (64 * R) a, F (normalizedSites B)) /
        (((64 * R).choose a : ℕ) : ℝ)) atTop
      (𝓝 ((a.factorial : ℝ) * ∫ x in orderedSiteDomain a, F x)) := by
  have hn : Tendsto (fun R : ℕ => 64 * R) atTop atTop :=
    tendsto_atTop_mono (fun R => by change R ≤ 64 * R; omega) tendsto_id
  have hchoose := (pow_div_choose_tendsto a).comp hn
  have hsite := (orderedSite_sum_tendsto F hF).comp hn
  have h := hchoose.mul hsite
  apply h.congr'
  filter_upwards [eventually_gt_atTop 0] with R hR
  have hN : (((64 * R : ℕ) : ℝ) ^ a) ≠ 0 := by
    positivity
  change (((64 * R : ℕ) : ℝ) ^ a /
      (((64 * R).choose a : ℕ) : ℝ)) *
      ((∑ B : BlockSubset (64 * R) a, F (normalizedSites B)) /
        (((64 * R : ℕ) : ℝ) ^ a)) =
    (∑ B : BlockSubset (64 * R) a, F (normalizedSites B)) /
      (((64 * R).choose a : ℕ) : ℝ)
  by_cases hc : ((((64 * R).choose a : ℕ) : ℝ)) = 0
  · simp [hc]
  · field_simp [hN, hc]

end Spin.UniformMixer
