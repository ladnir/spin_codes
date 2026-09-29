import SpinCodes.Structured.ConcretePlacementEncoderApprox

/-! Vanishing actual-to-finite-product error at each fixed occupation. -/
noncomputable section
namespace Spin.Structured.Placement
open Finset Routing ConcreteEncoder Filter
open scoped Topology

theorem regionKernelError_sqrt {a R H : Nat} (hR : 0 < R) (θ : ℝ) :
    regionKernelError a R H θ =
      ((a : ℝ) + 1) * 524288 * (θ * Real.sqrt 3 / Real.sqrt R + (1 / 2 : ℝ)^H) +
        (a : ℝ) * θ / R := by
  have hR0 : (R : ℝ) ≠ 0 := by exact_mod_cast Nat.ne_of_gt hR
  have hs0 : Real.sqrt (R : ℝ) ≠ 0 := (Real.sqrt_pos.mpr (by exact_mod_cast hR)).ne'
  have hs := Real.sq_sqrt (Nat.cast_nonneg R : (0 : ℝ) ≤ R)
  unfold regionKernelError
  rw [Real.sqrt_mul (by norm_num : (0 : ℝ) ≤ 3)]
  field_simp
  linear_combination (((a : ℝ) + 1) * 524288 * θ * Real.sqrt 3) * hs

theorem placementOmissionError_tendsto (a H : Nat) :
    Tendsto (fun R : Nat => placementOmissionError a R H) atTop (𝓝 0) := by
  have hR : Tendsto (fun R : Nat => (R : ℝ)) atTop atTop := tendsto_natCast_atTop_atTop
  have hd : Tendsto (fun R : Nat => 128 * (R : ℝ)) atTop atTop := hR.const_mul_atTop (by norm_num)
  have hd' : Tendsto (fun R : Nat => 128 * (R : ℝ) - 1) atTop atTop := by
    simpa only [sub_eq_add_neg] using tendsto_atTop_add_const_right atTop (-1 : ℝ) hd
  have h1 := hd.const_div_atTop (256 * (H : ℝ) * a)
  have h2 := hd'.const_div_atTop ((256 * ((H : ℝ) + 1) + 1) * (a : ℝ)^2)
  simpa only [placementOmissionError, add_zero] using h1.add h2

theorem regionKernelError_tendsto (a H : Nat) (θ : ℝ) :
    Tendsto (fun R : Nat => regionKernelError a R H θ) atTop
      (𝓝 (((a : ℝ) + 1) * 524288 * (1 / 2 : ℝ)^H)) := by
  have hR : Tendsto (fun R : Nat => (R : ℝ)) atTop atTop := tendsto_natCast_atTop_atTop
  have hs := (Real.tendsto_sqrt_atTop.comp hR).const_div_atTop (θ * Real.sqrt 3)
  have hl := hR.const_div_atTop ((a : ℝ) * θ)
  have h := ((hs.add_const ((1 / 2 : ℝ)^H)).const_mul (((a : ℝ) + 1) * 524288)).add hl
  simp only [zero_add, add_zero] at h
  apply h.congr'
  filter_upwards [eventually_gt_atTop 0] with R hRpos
  exact (regionKernelError_sqrt hRpos θ).symm

theorem totalPlacementError_tendsto (a H : Nat) (θ : ℝ) :
    Tendsto (fun R : Nat => placementOmissionError a R H + regionKernelError a R H θ) atTop
      (𝓝 (((a : ℝ) + 1) * 524288 * (1 / 2 : ℝ)^H)) := by
  simpa only [zero_add] using (placementOmissionError_tendsto a H).add (regionKernelError_tendsto a H θ)

end Spin.Structured.Placement

