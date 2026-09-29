import SpinCodes.Structured.ConcreteOuterTailDensePath

noncomputable section
namespace Spin.Structured.ConcreteOuter
open Filter

/-- Every fixed polynomial factor is dominated by the certified dense-tail gap. -/
theorem dense_polynomial_exp_tendsto (d : ℕ) :
    Tendsto (fun b : ℕ => ((b : ℝ)+1)^d *
      Real.exp (-(768/10^10 : ℝ)*(b : ℝ)+denseLogError b)) atTop (nhds 0) := by
  let c : ℝ := 768/10^10
  have hc : 0 < c := by norm_num [c]
  have hb : Tendsto (fun b : ℕ => (b : ℝ)+1) atTop atTop :=
    tendsto_atTop_add_const_right atTop 1 tendsto_natCast_atTop_atTop
  have hh := (tendsto_rpow_mul_exp_neg_mul_atTop_nhds_zero (d+48 : ℕ) c hc).comp hb
  simp only [Real.rpow_natCast] at hh
  have h := hh.mul_const (Real.exp (c+298))
  rw [zero_mul] at h
  apply h.congr
  intro b
  change (((b : ℝ)+1)^(d+48)*Real.exp (-c*((b : ℝ)+1)))*Real.exp (c+298) = _
  have hl : Real.exp (48*Real.log ((b : ℝ)+1)) = ((b : ℝ)+1)^48 := by
    rw [show (48 : ℝ) = (48 : ℕ) by norm_num, Real.exp_nat_mul,
      Real.exp_log (by positivity : (0 : ℝ) < (b : ℝ)+1)]
  have he : -(768/10^10 : ℝ)*(b : ℝ)+denseLogError b =
      48*Real.log ((b : ℝ)+1)+(-c*((b : ℝ)+1)+(c+298)) := by
    unfold denseLogError c
    ring
  rw [he]
  simp only [Real.exp_add, hl, pow_add]
  ring

theorem dense_nativeBlocks_exp_tendsto (d : ℕ) :
    Tendsto (fun k : ℕ => (((k*24 : ℕ) : ℝ)+1)^d *
      Real.exp (-(768/10^10 : ℝ)*((k*24 : ℕ) : ℝ)+denseLogError (k*24))) atTop (nhds 0) := by
  apply (dense_polynomial_exp_tendsto d).comp
  exact tendsto_atTop_mono (fun k => by omega : ∀ k : ℕ, k ≤ k*24) tendsto_id

#print axioms dense_polynomial_exp_tendsto
#print axioms dense_nativeBlocks_exp_tendsto

end Spin.Structured.ConcreteOuter
