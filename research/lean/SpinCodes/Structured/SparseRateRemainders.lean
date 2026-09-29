import SpinCodes.Structured.SparseFiniteExponentScalar

noncomputable section
namespace Spin.Structured.SparseRate

theorem log_conditioning_le {Q : ℝ} (hQ : 4096 ≤ Q) :
    Real.log (8*Real.sqrt Q) ≤ (1/500)*Q := by
  have hQp : 0 < Q := by linarith
  have hs : 0 < Real.sqrt Q := Real.sqrt_pos.2 hQp
  have hr := Real.log_le_sub_one_of_pos (by positivity : 0 < Q/4096)
  rw [Real.log_div hQp.ne' (by norm_num)] at hr
  have h4096 : Real.log (4096:ℝ) = 12*Real.log 2 := by
    have hh := Real.log_pow (2:ℝ) 12
    norm_num at hh
    exact hh
  have h8 : Real.log (8:ℝ) = 3*Real.log 2 := by
    have hh := Real.log_pow (2:ℝ) 3
    norm_num at hh
    exact hh
  rw [Real.log_mul (by norm_num) hs.ne', Real.log_sqrt hQp.le, h8]
  rw [h4096] at hr
  have h2 := log_two_lt_seven_tenths
  linarith

theorem log_2048_le {Q : ℝ} (hQ : 4096 ≤ Q) : Real.log 2048 ≤ Q := by
  have hh : Real.log (2048:ℝ) = 11*Real.log 2 := by
    have h := Real.log_pow (2:ℝ) 11
    norm_num at h
    exact h
  rw [hh]
  have h2 := log_two_lt_seven_tenths
  linarith

/-- Explicit uniform bound for the lower-order sparse terms. -/
theorem remainder_le {Q b ε : ℝ} (hQ : 4096 ≤ Q) (hb : 1000 ≤ b)
    (hε : ε ≤ 1/500) : ε + 1/b + Real.log 2048/(Q*b) ≤ 1/250 := by
  have hQp : 0 < Q := by linarith
  have hbp : 0 < b := by linarith
  have hi : 1/b ≤ (1/1000:ℝ) := by
    apply (div_le_iff₀ hbp).mpr
    linarith
  have hl : Real.log 2048/(Q*b) ≤ 1/b := by
    apply (div_le_div_iff₀ (mul_pos hQp hbp) hbp).mpr
    have hh := mul_le_mul_of_nonneg_right (log_2048_le hQ) hbp.le
    nlinarith
  linarith

end Spin.Structured.SparseRate
