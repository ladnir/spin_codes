import SpinCodes.Structured.ConcretePlacementLimitKernel

/-! Uniform componentwise convergence of the actual finite shuffled IMT region experiment. -/
noncomputable section
namespace Spin.Structured.Placement
open Finset Routing ConcreteEncoder Filter MeasureTheory
open scoped Topology

def kernelComparisonError (a R H : Nat) (θ : ℝ) : ℝ :=
  (placementOmissionError a R H + regionKernelError a R H θ) +
    ((θ * epochMean / 128) * ((a : ℝ) / R) + placementOmissionError a R H)

theorem kernelComparisonError_tendsto (a H : Nat) (θ : ℝ) :
    Tendsto (fun R : Nat => kernelComparisonError a R H θ) atTop
      (𝓝 (((a : ℝ) + 1) * 524288 * (1 / 2 : ℝ)^H)) := by
  have hn : Tendsto (fun R : Nat => (R : ℝ)) atTop atTop := tendsto_natCast_atTop_atTop
  have hd := (hn.const_div_atTop (a : ℝ)).const_mul (θ * epochMean / 128)
  simpa only [kernelComparisonError, mul_zero, add_zero, zero_add] using
    (totalPlacementError_tendsto a H θ).add (hd.add (placementOmissionError_tendsto a H))

theorem shuffledRegionKernel_uniform_limit (a : Nat) {θ : ℝ} (hθ : 0 ≤ θ)
    {ε : ℝ} (hε : 0 < ε) :
    ∀ᶠ R : Nat in atTop, ∀ S : Finset (Fin (128 * R)), S.card = a → ∀ q r : State,
      |shuffledRegionKernel θ S q r - (liveLift * continuumRegionKernel θ a * liveProjection) q r| < ε := by
  have hg : Tendsto (fun H : Nat => ((a : ℝ) + 1) * 524288 * (1 / 2 : ℝ)^H) atTop (𝓝 0) := by
    simpa only [mul_zero] using
      (tendsto_pow_atTop_nhds_zero_of_lt_one (by norm_num : (0 : ℝ) ≤ 1 / 2)
        (by norm_num : (1 / 2 : ℝ) < 1)).const_mul (((a : ℝ) + 1) * 524288)
  obtain ⟨H, hH⟩ := (hg.eventually (gt_mem_nhds (half_pos hε))).exists
  have he := (kernelComparisonError_tendsto a H θ).eventually (gt_mem_nhds hH)
  have hf : ∀ᶠ R : Nat in atTop, ∀ q r : State,
      |finiteSiteKernel θ R a q r - (liveLift * continuumRegionKernel θ a * liveProjection) q r| < ε / 2 := by
    apply Filter.eventually_all.mpr
    intro q
    apply Filter.eventually_all.mpr
    intro r
    have hh : Tendsto (fun R : Nat =>
        |finiteSiteKernel θ R a q r - (liveLift * continuumRegionKernel θ a * liveProjection) q r|)
        atTop (𝓝 0) := by
      simpa only [sub_self, abs_zero] using ((finiteSiteKernel_tendsto θ a q r).sub_const
        ((liveLift * continuumRegionKernel θ a * liveProjection) q r)).abs
    exact hh.eventually (gt_mem_nhds (half_pos hε))
  filter_upwards [he, hf, eventually_gt_atTop 0, eventually_ge_atTop H] with R he hf hR hHR
  intro S hS q r
  have h1 := (shuffledRegionKernel_site_error hR hθ S hS H hHR q r).trans_lt he
  have h2 := hf q r
  exact (abs_sub_le _ (finiteSiteKernel θ R a q r) _).trans_lt (by linarith)

theorem shuffledRegionKernel_tendsto {a : Nat} {θ : ℝ} (hθ : 0 ≤ θ)
    (S : (R : Nat) → Finset (Fin (128 * R))) (hS : ∀ᶠ R : Nat in atTop, (S R).card = a)
    (q r : State) :
    Tendsto (fun R => shuffledRegionKernel θ (S R) q r) atTop
      (𝓝 ((liveLift * continuumRegionKernel θ a * liveProjection) q r)) := by
  apply Metric.tendsto_nhds.mpr
  intro ε hε
  filter_upwards [shuffledRegionKernel_uniform_limit a hθ hε, hS] with R hR hSR
  simpa only [Real.dist_eq] using hR (S R) hSR q r

end Spin.Structured.Placement
