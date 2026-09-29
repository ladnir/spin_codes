import SpinCodes.Structured.ConcretePlacementLimitError

/-! The actual finite encoder has arbitrarily small error against its coarse region sum. -/
noncomputable section
namespace Spin.Structured.Placement
open Routing ConcreteEncoder Filter
open scoped Topology

theorem actual_region_coarse_eventually (a : Nat) {θ : ℝ} (hθ : 0 ≤ θ)
    {ε : ℝ} (hε : 0 < ε) :
    ∃ H : Nat, ∀ᶠ R : Nat in atTop, ∀ S : Finset (Fin (128 * R)), S.card = a → ∀ q : State,
      |(shuffleLaw S).expect (fun T => regionMoment (Real.exp (-(θ / (128 * R)))) T q) -
        finiteGoodMoment θ R a H q| < ε := by
  have hg : Tendsto (fun H : Nat => ((a : ℝ) + 1) * 524288 * (1 / 2 : ℝ)^H) atTop (𝓝 0) := by
    simpa only [mul_zero] using
      (tendsto_pow_atTop_nhds_zero_of_lt_one (by norm_num : (0 : ℝ) ≤ 1 / 2)
        (by norm_num : (1 / 2 : ℝ) < 1)).const_mul (((a : ℝ) + 1) * 524288)
  obtain ⟨H, hH⟩ := (hg.eventually (gt_mem_nhds hε)).exists
  refine ⟨H, ?_⟩
  have he := (totalPlacementError_tendsto a H θ).eventually (gt_mem_nhds hH)
  filter_upwards [he, eventually_gt_atTop 0, eventually_ge_atTop H] with R he hR hHR
  intro S hS q
  exact (shuffled_regionMoment_coarse_error hR hθ S hS H hHR q).trans_lt he

end Spin.Structured.Placement
