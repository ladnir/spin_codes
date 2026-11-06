import SpinCodes.UniformMixer.RegionalAverageContinuum

/-!
# Uniform asymptotic regional comparison

At fixed support size, the finite rounded-site error, the multiplicative
empty-gap comparison, and the discarded-placement probability vanish as
the region grows. The ideal support grid converges to the paper's
ordered-simplex endpoint.
-/

noncomputable section

namespace Spin.UniformMixer

open Spin.Structured.Placement Filter
open scoped Topology

def regionalUpper (θ : ℝ) (a R : ℕ) (q r : State22) : ℝ :=
  (1 + ((64 * θ + 64 * θ ^ 2) / (64 * (R : ℝ))) /
      Real.exp (-θ)) ^ (a + 1) *
    ((∑ B : BlockSubset (64 * R) a,
      (liftM * continuumSiteProduct (regionGamma θ)
        (normalizedSites B) * projectM) q r) /
        (((64 * R).choose a : ℕ) : ℝ) +
      regionGamma θ * (3 * (a + 1 : ℕ) / (R : ℝ))) +
    128 * (a : ℝ) / (64 * (R : ℝ)) +
    257 * (a : ℝ) ^ 2 / (((64 * R : ℕ) : ℝ) - 1)

theorem regionalUpper_tendsto (θ : ℝ) (a : ℕ) (q r : State22) :
    Tendsto (fun R : ℕ => regionalUpper θ a R q r) atTop
      (𝓝 (continuumEndpoint θ a q r)) := by
  have hn : Tendsto (fun R : ℕ => (R : ℝ)) atTop atTop :=
    tendsto_natCast_atTop_atTop
  have hd : Tendsto (fun R : ℕ => 64 * (R : ℝ)) atTop atTop :=
    hn.const_mul_atTop (by norm_num)
  have hd' : Tendsto (fun R : ℕ => 64 * (R : ℝ) - 1) atTop atTop := by
    simpa only [sub_eq_add_neg] using
      tendsto_atTop_add_const_right atTop (-1 : ℝ) hd
  have hC : Tendsto (fun R : ℕ =>
      ((64 * θ + 64 * θ ^ 2) / (64 * (R : ℝ))) /
        Real.exp (-θ)) atTop (𝓝 0) := by
    simpa only [zero_div] using
      (hd.const_div_atTop (64 * θ + 64 * θ ^ 2)).div_const
        (Real.exp (-θ))
  have hpow : Tendsto (fun R : ℕ =>
      (1 + ((64 * θ + 64 * θ ^ 2) / (64 * (R : ℝ))) /
        Real.exp (-θ)) ^ (a + 1)) atTop (𝓝 1) := by
    simpa only [add_zero, one_pow] using (hC.const_add 1).pow (a + 1)
  have hδ : Tendsto (fun R : ℕ =>
      regionGamma θ * (3 * (a + 1 : ℕ) / (R : ℝ))) atTop (𝓝 0) := by
    have hh := hn.const_div_atTop
      (regionGamma θ * (3 * (a + 1 : ℕ)))
    apply hh.congr'
    filter_upwards [eventually_gt_atTop 0] with R hR
    ring
  have hbad1 : Tendsto (fun R : ℕ =>
      128 * (a : ℝ) / (64 * (R : ℝ))) atTop (𝓝 0) :=
    hd.const_div_atTop (128 * (a : ℝ))
  have hbad2 : Tendsto (fun R : ℕ =>
      257 * (a : ℝ) ^ 2 / (((64 * R : ℕ) : ℝ) - 1))
      atTop (𝓝 0) := by
    simpa only [Nat.cast_mul, Nat.cast_ofNat] using
      hd'.const_div_atTop (257 * (a : ℝ) ^ 2)
  have hgrid := idealRegionalGrid_tendsto θ a q r
  have h := ((hpow.mul (hgrid.add hδ)).add hbad1).add hbad2
  simpa only [regionalUpper, one_mul, add_zero] using h

theorem averagedRegionKernel_uniform_continuum_upper (a : ℕ)
    {θ : ℝ} (hθ : 0 ≤ θ) {ε : ℝ} (hε : 0 < ε) :
    ∀ᶠ R : ℕ in atTop,
      ∀ S : Finset (Fin (64 * R)), S.card = a →
        ∀ q r : State22,
          averagedRegionKernel S
            (Real.exp (-(θ / (64 * (R : ℝ))))) q r <
          continuumEndpoint θ a q r + ε := by
  have hlim : ∀ᶠ R : ℕ in atTop, ∀ q r : State22,
      regionalUpper θ a R q r < continuumEndpoint θ a q r + ε := by
    apply Filter.eventually_all.mpr
    intro q
    apply Filter.eventually_all.mpr
    intro r
    exact (regionalUpper_tendsto θ a q r).eventually
      (gt_mem_nhds (by linarith :
        continuumEndpoint θ a q r < continuumEndpoint θ a q r + ε))
  filter_upwards [hlim, eventually_gt_atTop 0,
    eventually_ge_atTop a] with R hlim hR ha
  intro S hS q r
  have ha' : a ≤ 64 * R := by omega
  have hfinite := averagedRegionKernel_continuum_finite
    hR ha' S hS θ hθ q r
  exact hfinite.trans_lt (hlim q r)

end Spin.UniformMixer
