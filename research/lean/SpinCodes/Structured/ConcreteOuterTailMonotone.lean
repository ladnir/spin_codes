import SpinCodes.Structured.ConcreteOuterTailMoment

noncomputable section
namespace Spin.Structured.ConcreteOuter
open Spin.Numeric

/-- The logarithmic sparse moment decreases throughout the paper's finite range. -/
theorem sparse_log_moment_antitone {x y : ℝ} (hy : 0 < y) (hyx : y ≤ x)
    (hx : x ≤ 3/250) :
    x*Real.log (kappa*x/(1-2*x)) ≤ y*Real.log (kappa*y/(1-2*y)) := by
  have hx0 : 0<x := hy.trans_le hyx
  have hdx : 0<1-2*x := by linarith
  have hdy : 0<1-2*y := by linarith
  have hk := kappa_pos
  have hAx : 0<kappa*x/(1-2*x) := by positivity
  have hAy : 0<kappa*y/(1-2*y) := by positivity
  have hAmax : kappa*x/(1-2*x) ≤ logArg := by
    unfold logArg
    apply (div_le_div_iff₀ hdx (by norm_num)).mpr
    nlinarith
  have hlogmax := Real.log_le_log hAx hAmax
  have hratio : 2*y/(1-2*x) ≤ 2*(3/250)/(1-2*(3/250:ℝ)) := by
    apply (div_le_div_iff₀ hdx (by norm_num)).mpr
    nlinarith
  have hcoeff : Real.log (kappa*x/(1-2*x))+1+2*y/(1-2*x) ≤ 0 := by
    have hh := logDeriv_lt
    linarith
  have hlog₁ := Real.log_le_sub_one_of_pos (show 0<x/y by positivity)
  have hlog₂ := Real.log_le_sub_one_of_pos (show 0<(1-2*y)/(1-2*x) by positivity)
  have hm₁ := mul_le_mul_of_nonneg_left hlog₁ hy.le
  have hm₂ := mul_le_mul_of_nonneg_left hlog₂ hy.le
  have he₁ : y*(x/y-1)=x-y := by field_simp
  have he₂ : y*((1-2*y)/(1-2*x)-1)=(x-y)*(2*y/(1-2*x)) := by field_simp; ring
  rw [he₁] at hm₁
  rw [he₂] at hm₂
  have hlogid : Real.log (kappa*x/(1-2*x)) - Real.log (kappa*y/(1-2*y)) =
      Real.log (x/y)+Real.log ((1-2*y)/(1-2*x)) := by
    rw [Real.log_div (by positivity) (ne_of_gt hdx), Real.log_div (by positivity) (ne_of_gt hdy),
      Real.log_mul (ne_of_gt hk) (ne_of_gt hx0), Real.log_mul (ne_of_gt hk) (ne_of_gt hy),
      Real.log_div (ne_of_gt hx0) (ne_of_gt hy), Real.log_div (ne_of_gt hdy) (ne_of_gt hdx)]
    ring
  have hm := mul_nonpos_of_nonneg_of_nonpos (sub_nonneg.mpr hyx) hcoeff
  nlinarith

end Spin.Structured.ConcreteOuter
