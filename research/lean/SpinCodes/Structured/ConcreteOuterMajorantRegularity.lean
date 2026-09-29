import SpinCodes.Majorant.Refined

noncomputable section
namespace Spin.Structured.ConcreteOuter
open Spin.Numeric

lemma refined_slopes_checked : ∀ p ∈ Spin.Majorant.supports, |p.1| ≤ (2:ℚ) := by
  decide +kernel

/-- The certified affine envelope changes by at most twice the input displacement. -/
theorem refined_step (x y : ℝ) :
    Spin.Majorant.refined.toFun x ≤ Spin.Majorant.refined.toFun y + 2*|x-y| := by
  obtain ⟨p,hp,he⟩ := Finset.exists_mem_eq_inf' Spin.Majorant.refined.nonempty
    (fun p : ℚ×ℚ => (p.1:ℝ)*y+(p.2:ℝ))
  have hs : |(p.1:ℝ)| ≤ 2 := by exact_mod_cast refined_slopes_checked p hp
  have ht := Spin.Majorant.refined.toFun_le_support hp x
  have hm : (p.1:ℝ)*(x-y) ≤ 2*|x-y| :=
    (le_abs_self _).trans ((abs_mul _ _).le.trans (mul_le_mul_of_nonneg_right hs (abs_nonneg _)))
  change Spin.Majorant.refined.toFun y = (p.1:ℝ)*y+(p.2:ℝ) at he
  rw [he]
  linarith

lemma refined_nonneg {w : ℝ} (hw : w ∈ Set.Icc ((13:ℝ)/125) (112/125)) :
    0 ≤ Spin.Majorant.refined.toFun w := by
  have hh := Spin.Majorant.majorant_pointwise (a:=0) (b:=0)
    (by norm_num) (by norm_num) hw (by constructor <;> norm_num; constructor <;> linarith [hw.1,hw.2])
  have hg := gBA_nonneg (a:=0) (by norm_num) (by norm_num)
  simp only [piBA, zero_div, hEnt_zero, mul_zero, sub_zero, add_zero] at hh
  exact hg.trans hh

end Spin.Structured.ConcreteOuter
