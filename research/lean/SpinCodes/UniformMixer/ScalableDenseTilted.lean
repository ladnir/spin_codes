import SpinCodes.UniformMixer.ScalableIidWiring
import SpinCodes.UniformMixer.ScalableTransferPaperOrientation
import SpinCodes.Structured.PositiveFiniteTilted
import SpinCodes.UniformMixer.ScalableFixedRoutedProfile

/-!
# Dense occupation: actual route at a tilted reference density

For a fixed outer word, the route preserves total weight. This lets the
existing exact iid likelihood ratio change the Bernoulli density before
applying the actual 64-bit transfer bound. The numerical dense-box cover
is a separate obligation.
-/

noncomputable section

namespace Spin.UniformMixer

open Spin.Structured Spin.Structured.ConcreteRoute
  Spin.Structured.PositiveFinite

theorem scalable_dense_tilted_moment (m : ℕ)
    (rows : Fin (Lsched m) → Finset (Fin (bsched m)))
    (hw0 : 0 < totalWeight rows)
    (hw1 : totalWeight rows < Lsched m * bsched m)
    {β z c ρ : ℝ} (hβ0 : 0 < β) (hβ1 : β < 1)
    (hz0 : 0 ≤ z) (hz1 : z ≤ 1)
    (v : Fin 3 → ℝ) (hc : 0 ≤ c) (hρ : 0 ≤ ρ)
    (hones : ∀ i, (1 : ℝ) ≤ c * v i)
    (hv : ∀ i,
      (paperSpectrumTransfer β z).mulVec v i ≤ ρ * v i) :
    (law rows).expect (scalableRoutedMoment m z) ≤
      ((((bsched m : ℝ) + 1) ^ activeRows rows *
        ((Lsched m : ℝ) + 1) ^ bsched m) *
        Real.exp (((Lsched m : ℝ) * bsched m) *
          binKL (density rows) β)) *
        (c * ρ ^ (bsched m * scalableStepsPerRegion m) * v 0) := by
  have hdom := routed_dominates_tilted
    (Lsched_pos m) (bsched_pos m) rows hw0 hw1 hβ0 hβ1
  have hiid := scalable_iid_moment_power m β z hβ0.le hβ1.le
  calc
    _ ≤ ((((bsched m : ℝ) + 1) ^ activeRows rows *
          ((Lsched m : ℝ) + 1) ^ bsched m) *
          Real.exp (((Lsched m : ℝ) * bsched m) *
            binKL (density rows) β)) *
          (iidLaw (Lsched m) (bsched m) β hβ0.le hβ1.le).expect
            (scalableRoutedMoment m z) := by
        rw [← permutation_law_eq rows]
        exact hdom.expect_le (scalableRoutedMoment_nonneg m hz0)
    _ = _ := by rw [hiid]
    _ ≤ _ := by
      apply mul_le_mul_of_nonneg_left
      · exact scalable_bernoulli_paper_collatz_bound β z
          hβ0.le hβ1.le hz0 hz1 v c ρ hc hρ hones hv _
      · positivity

theorem scalable_dense_tilted_probability (m d : ℕ)
    (rows : Fin (Lsched m) → Finset (Fin (bsched m)))
    (hw0 : 0 < totalWeight rows)
    (hw1 : totalWeight rows < Lsched m * bsched m)
    {β z c ρ : ℝ} (hβ0 : 0 < β) (hβ1 : β < 1)
    (hz0 : 0 < z) (hz1 : z ≤ 1)
    (v : Fin 3 → ℝ) (hc : 0 ≤ c) (hρ : 0 ≤ ρ)
    (hones : ∀ i, (1 : ℝ) ≤ c * v i)
    (hv : ∀ i,
      (paperSpectrumTransfer β z).mulVec v i ≤ ρ * v i) :
    ((law rows).prod (scalableStateMapLaw m)).prob
      (fun ω => sampledOutputWeight
        (scalableFlatInput m ω.1) (scalableFlatMaps m ω.2) 0 ≤ d) ≤
      (((((bsched m : ℝ) + 1) ^ activeRows rows *
          ((Lsched m : ℝ) + 1) ^ bsched m) *
          Real.exp (((Lsched m : ℝ) * bsched m) *
            binKL (density rows) β)) *
          (c * ρ ^ (bsched m * scalableStepsPerRegion m) * v 0)) /
        z ^ d := by
  apply FinPMF.prob_weight_le_of_moment _ _ d hz0 hz1
  rw [FinPMF.expect_prod]
  simp_rw [scalable_fixed_route_moment]
  exact scalable_dense_tilted_moment m rows hw0 hw1 hβ0 hβ1
    hz0.le hz1 v hc hρ hones hv

#print axioms scalable_dense_tilted_moment
#print axioms scalable_dense_tilted_probability

end Spin.UniformMixer
