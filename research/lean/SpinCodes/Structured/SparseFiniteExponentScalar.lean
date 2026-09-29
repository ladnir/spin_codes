import SpinCodes.Structured.KLChain
import SpinCodes.Numeric.RatLog

noncomputable section
namespace Spin.Structured.SparseRate

def zmin : ℝ := 6249/6250

theorem log_five_eighths_upper : Real.log (5/8:ℝ) ≤ -(470003629:ℝ)/1000000000 := by
  have h := Spin.Numeric.real_log_le_logUpper (q := (5/8:ℚ)) (by norm_num)
    (k := 0) (n := 26) (by norm_num)
  have hn : Spin.Numeric.logUpper (5/8) 0 26 ≤ -(470003629:ℚ)/1000000000 := by
    norm_num [Spin.Numeric.logUpper, Spin.Numeric.ratLogSeries, Spin.Numeric.ratLogRem,
      Finset.sum_range_succ]
  have hnR : ((Spin.Numeric.logUpper (5/8) 0 26 : ℚ):ℝ) ≤
      -(470003629:ℝ)/1000000000 := by exact_mod_cast hn
  have hR : Real.log (5/8:ℝ) ≤ ((Spin.Numeric.logUpper (5/8) 0 26:ℚ):ℝ) := by simpa using h
  exact hR.trans hnR

def leadingConstant : ℝ := Real.log 2/2 + Real.log (5/8:ℝ) +
  (3/5 + (11/100)*(8/5))/zmin - 96/128 + 1281/100000 + 4*Real.log 2/39

theorem leadingConstant_lt : leadingConstant < -(1340384:ℝ)/100000000 := by
  have h2 := Spin.Numeric.lt_log2Hi
  norm_num [Spin.Numeric.log2Hi] at h2
  have h5 := log_five_eighths_upper
  unfold leadingConstant zmin
  norm_num
  linarith

theorem log_two_lt_seven_tenths : Real.log 2 < (7/10:ℝ) := by
  have hh := Spin.Numeric.lt_log2Hi
  norm_num [Spin.Numeric.log2Hi] at hh
  linarith

theorem sparse_kl_le {α : ℝ} (hα0 : 0 < α) (hα1 : α ≤ 1/10000) :
    Spin.binKL α ((8/5)*α) ≤ α*(Real.log (5/8:ℝ)+(3/5)/zmin) := by
  have ha1 : 0 < 1-α := by linarith
  have hz : 0 < 1-(8/5:ℝ)*α := by linarith
  have hzmin : zmin ≤ 1-(8/5:ℝ)*α := by unfold zmin; linarith
  have hratio : 0 < (1-α)/(1-(8/5:ℝ)*α) := div_pos ha1 hz
  have hl := Real.log_le_sub_one_of_pos hratio
  have ha : α/((8/5:ℝ)*α) = 5/8 := by field_simp
  have he : (1-α)/(1-(8/5:ℝ)*α)-1 = ((3/5)*α)/(1-(8/5:ℝ)*α) := by
    field_simp [show (5:ℝ)-α*8 ≠ 0 by linarith]
    <;> ring
  rw [he] at hl
  have hd : ((3/5)*α)/(1-(8/5:ℝ)*α) ≤ ((3/5)*α)/zmin :=
    div_le_div_of_nonneg_left (by positivity) (by norm_num [zmin]) hzmin
  have hnn : 0 ≤ ((3/5:ℝ)*α)/zmin := by unfold zmin; positivity
  have hh := mul_le_mul_of_nonneg_left (hl.trans hd) ha1.le
  have hh2 : (1-α)*(((3/5:ℝ)*α)/zmin) ≤ ((3/5)*α)/zmin := by nlinarith
  unfold Spin.binKL
  rw [ha]
  calc α*Real.log (5/8:ℝ)+(1-α)*Real.log ((1-α)/(1-(8/5)*α))
      ≤ α*Real.log (5/8:ℝ)+((3/5)*α)/zmin := add_le_add le_rfl (hh.trans hh2)
    _ = _ := by ring

theorem sparse_neg_log_le {α : ℝ} (hα0 : 0 < α) (hα1 : α ≤ 1/10000) :
    -Real.log (1-(8/5:ℝ)*α) ≤ ((8/5)*α)/zmin := by
  have hz : 0 < 1-(8/5:ℝ)*α := by linarith
  have hzmin : zmin ≤ 1-(8/5:ℝ)*α := by unfold zmin; linarith
  have hh := Real.one_sub_inv_le_log_of_pos hz
  have he : (1-(8/5:ℝ)*α)⁻¹-1 = ((8/5)*α)/(1-(8/5:ℝ)*α) := by
    apply (eq_div_iff hz.ne').mpr
    rw [sub_mul, inv_mul_cancel₀ hz.ne']
    ring
  have hd : ((8/5)*α)/(1-(8/5:ℝ)*α) ≤ ((8/5)*α)/zmin :=
    div_le_div_of_nonneg_left (by positivity) (by norm_num [zmin]) hzmin
  linarith

theorem sparse_log_contraction_le {α : ℝ} (hα0 : 0 < α) (hα1 : α ≤ 1/10000) :
    Real.log (1-96*α) ≤ -96*α := by
  have hh := Real.log_le_sub_one_of_pos (by linarith : 0 < 1-96*α)
  linarith

end Spin.Structured.SparseRate


