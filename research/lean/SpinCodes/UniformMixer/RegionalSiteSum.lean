import SpinCodes.UniformMixer.RegionalAverage
import SpinCodes.Structured.ConcreteFixedSubsetMarginal

/-!
# Exact finite support sum for a uniform regional shuffle

The previous module bounds the average under the actual shuffle law.
This module rewrites that expectation as a uniform sum over fixed-size
subsets of the 64-bit-step region. The resulting sum is the finite grid
that must be compared with the paper's ordered-simplex integral.
-/

noncomputable section

namespace Spin.UniformMixer

open Spin.Structured.Routing Spin.Structured.Placement
attribute [local instance] Classical.propDecidable

theorem averagedRegionKernel_layer {R a : ℕ}
    (S : Finset (Fin (64 * R))) (hS : S.card = a)
    (z : ℝ) (q r : State22) :
    averagedRegionKernel S z q r =
      (∑ T : BlockSubset (64 * R) a,
        chronologicalRegionKernel z T.val q r) /
        (((64 * R).choose a : ℕ) : ℝ) := by
  exact shuffleLaw_expect_layer S hS
    (fun T => chronologicalRegionKernel z T q r)

theorem averagedRegionKernel_coarse_layer {R a : ℕ} (hR : 0 < R)
    (S : Finset (Fin (64 * R))) (hS : S.card = a)
    (θ : ℝ) (hθ : 0 ≤ θ) (q r : State22) :
    let L : ℝ := 64 * R
    averagedRegionKernel S (Real.exp (-(θ / L))) q r ≤
      (∑ T : BlockSubset (64 * R) a,
        if ¬ badRegionalSites T.val then
          coarseRegionalEnvelope θ T.val q r else 0) /
        (((64 * R).choose a : ℕ) : ℝ) +
      128 * (a : ℝ) / L +
      257 * (a : ℝ) ^ 2 / (((64 * R : ℕ) : ℝ) - 1) := by
  have h := averagedRegionKernel_coarse_upper hR S θ hθ q r
  rw [shuffleLaw_expect_layer S hS] at h
  simpa only [hS] using h

end Spin.UniformMixer
