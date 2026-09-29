import SpinCodes.Structured.SparseRateRemainders

noncomputable section
namespace Spin.Structured.SparseRate

/-- Logarithm of the finite sparse estimate after the elementary bound on choose(L,Q). -/
def exponent (L Q b R d α ε : ℝ) : ℝ :=
  Real.log 2048 + Q*Real.log L + Q + Q*b*(Real.log 2/2+1281/100000+ε) +
    b*Real.log (8*Real.sqrt Q) + L*b*Spin.binKL α ((8/5)*α) +
    R*Real.log (1-96*α) - d*Real.log (1-(8/5)*α)

/-- Explicit finite sparse exponent. The only spectral input is ε≤.002. -/
theorem exponent_le {L Q b R d α ε : ℝ}
    (hL : 0 < L) (hQ : 4096 ≤ Q) (hb : 1000 ≤ b)
    (hR : 0 ≤ R) (hd : 0 ≤ d) (hα0 : 0 < α) (hα1 : α ≤ 1/10000)
    (hα : L*α = Q) (hrounds : 128*R = L*b) (hdmax : d ≤ (11/100)*L*b)
    (hschedule : Real.log L ≤ (4*Real.log 2/39)*b) (hε : ε ≤ 1/500) :
    exponent L Q b R d α ε ≤ -(3/500)*Q*b := by
  have hQp : 0 < Q := by linarith
  have hbp : 0 < b := by linarith
  have hQb : 0 < Q*b := mul_pos hQp hbp
  have hNα : L*b*α = Q*b := by
    calc L*b*α = (L*α)*b := by ring
         _ = Q*b := by rw [hα]
  have hRα : R*α = Q*b/128 := by nlinarith [congrArg (fun x : ℝ => x*α) hrounds]
  have hchoose : Q*Real.log L ≤ Q*b*(4*Real.log 2/39) := by
    simpa only [mul_assoc, mul_comm, mul_left_comm] using
      mul_le_mul_of_nonneg_left hschedule hQp.le
  have hcond : b*Real.log (8*Real.sqrt Q) ≤ (1/500)*Q*b := by
    have hh := mul_le_mul_of_nonneg_left (log_conditioning_le hQ) hbp.le
    nlinarith
  have hkl : L*b*Spin.binKL α ((8/5)*α) ≤ Q*b*(Real.log (5/8:ℝ)+(3/5)/zmin) := by
    have hh := mul_le_mul_of_nonneg_left (sparse_kl_le hα0 hα1) (mul_pos hL hbp).le
    simpa only [← mul_assoc, hNα] using hh
  have hcontract : R*Real.log (1-96*α) ≤ -(96/128)*Q*b := by
    have hh := mul_le_mul_of_nonneg_left (sparse_log_contraction_le hα0 hα1) hR
    nlinarith
  have hden : 0 ≤ ((8/5:ℝ)*α)/zmin := by unfold zmin; positivity
  have hmarkov : -d*Real.log (1-(8/5)*α) ≤ Q*b*((11/100)*(8/5)/zmin) := by
    have hh := mul_le_mul_of_nonneg_left (sparse_neg_log_le hα0 hα1) hd
    have hh2 := mul_le_mul_of_nonneg_right hdmax hden
    have he : ((11/100)*L*b)*(((8/5)*α)/zmin) = Q*b*((11/100)*(8/5)/zmin) := by
      calc ((11/100)*L*b)*(((8/5)*α)/zmin) = (L*b*α)*((11/100)*(8/5)/zmin) := by ring
        _ = _ := by rw [hNα]
    rw [he] at hh2
    nlinarith
  have hlog : Real.log 2048 ≤ Q := log_2048_le hQ
  have hsmall : ε*(Q*b)+Q+Real.log 2048 ≤ (1/250)*(Q*b) := by
    have he := mul_le_mul_of_nonneg_right hε hQb.le
    have hbQ := mul_le_mul_of_nonneg_left hb hQp.le
    nlinarith
  have hlead : leadingConstant*(Q*b) < (-(1340384:ℝ)/100000000)*(Q*b) :=
    mul_lt_mul_of_pos_right leadingConstant_lt hQb
  unfold exponent
  unfold leadingConstant at hlead
  norm_num [zmin] at hkl hmarkov hlead
  nlinarith

end Spin.Structured.SparseRate


