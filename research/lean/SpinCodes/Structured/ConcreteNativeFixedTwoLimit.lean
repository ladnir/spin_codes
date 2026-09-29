import SpinCodes.Structured.ConcreteNativeFixedTwoBound
import Mathlib.Analysis.Complex.ExponentialBounds

noncomputable section
namespace Spin.Structured.ConcreteNativeFamily
open Filter

theorem fixed_two_log_gap :
    8*Real.log 2/39+2*(Real.log 2/2+1281/100000)+Real.log (1011/4000)+11/25 ≤ -(7/100:ℝ) := by
  have hh := Real.log_le_sub_one_of_pos (by norm_num : (0:ℝ)<1011/1000)
  have he : Real.log (1011/4000:ℝ)=Real.log (1011/1000)-2*Real.log 2 := by
    rw [show (2:ℝ)*Real.log 2=Real.log 4 by rw [show (4:ℝ)=2^2 by norm_num,Real.log_pow]; norm_num]
    rw [←Real.log_div (by norm_num : (1011/1000:ℝ)≠0) (by norm_num : (4:ℝ)≠0)]
    congr 1
    norm_num
  rw [he]
  linarith [Real.log_two_gt_d9]

theorem fixedTwoBound_exp (m : ℕ) :
    fixedTwoBound m = 1000*Real.exp (2*Real.log (Lsched m:ℝ)+
      2*(bsched m:ℝ)*(Real.log 2/2+1281/100000+selectedRemainder m (nativeShellRemainder m))+
      (bsched m:ℝ)*Real.log (1011/4000)+(4/(Lsched m:ℝ))*(threshold m:ℝ)) := by
  have hL : (0:ℝ)<Lsched m := by exact_mod_cast Lsched_pos m
  have hp : Real.exp ((bsched m:ℝ)*Real.log (1011/4000))=(1011/4000:ℝ)^(bsched m) := by
    rw [Real.exp_nat_mul,Real.exp_log (by norm_num : (0:ℝ)<1011/4000)]
  have hz : Real.exp (-(4/(Lsched m:ℝ)))^(threshold m)=
      (Real.exp ((4/(Lsched m:ℝ))*(threshold m:ℝ)))⁻¹ := by
    rw [←Real.exp_nat_mul,←Real.exp_neg]
    congr 1
    ring
  unfold fixedTwoBound
  rw [hz,div_inv_eq_mul]
  have hLp : Real.exp (2*Real.log (Lsched m:ℝ))=(Lsched m:ℝ)^2 := by rw [show (2:ℝ)=((2:ℕ):ℝ) by norm_num,Real.exp_nat_mul,Real.exp_log hL]
  simp only [Real.exp_add,hLp,hp]
  ring

theorem fixedTwoBound_eventually : ∀ᶠ m in atTop,
    fixedTwoBound m≤1000*Real.exp (-(1/20:ℝ)*(bsched m:ℝ)) := by
  have hr := (selectedRemainder_tendsto nativeShellRemainder_tendsto).eventually
    (gt_mem_nhds (by norm_num : (0:ℝ)<1/100))
  filter_upwards [hr] with m hm
  rw [fixedTwoBound_exp]
  apply mul_le_mul_of_nonneg_left _ (by norm_num)
  apply Real.exp_le_exp.mpr
  have hL : (0:ℝ)<Lsched m := by exact_mod_cast Lsched_pos m
  have hb : (0:ℝ)≤bsched m := Nat.cast_nonneg _
  have hth := mul_le_mul_of_nonneg_left (threshold_le m) (show 0≤4/(Lsched m:ℝ) by positivity)
  have hc : (4/(Lsched m:ℝ))*((11/100)*(Lsched m:ℝ)*(bsched m:ℝ))=(11/25)*(bsched m:ℝ) := by
    field_simp
    ring
  rw [hc] at hth
  have hgap := mul_le_mul_of_nonneg_right fixed_two_log_gap hb
  have hs := SparseRate.native_log_schedule m
  have he := mul_le_mul_of_nonneg_left hm.le hb
  nlinarith

theorem native_two_EZ_tendsto : Tendsto (fun m => concreteFamily.EZ m 2) atTop (nhds 0) := by
  have hb : Tendsto (fun m => (bsched m:ℝ)) atTop atTop :=
    tendsto_natCast_atTop_atTop.comp bsched_tendsto
  have he := (Real.tendsto_exp_atBot.comp
    (hb.const_mul_atTop_of_neg (by norm_num : -(1/20:ℝ)<0))).const_mul 1000
  simp only [mul_zero] at he
  apply tendsto_of_tendsto_of_tendsto_of_le_of_le' tendsto_const_nhds he
  · exact Eventually.of_forall (fun m => concreteFamily.EZ_nonneg m 2)
  · filter_upwards [native_two_EZ_bound,fixedTwoBound_eventually] with m hm hbound
    exact hm.trans hbound

#print axioms fixed_two_log_gap
#print axioms fixedTwoBound_eventually
#print axioms native_two_EZ_tendsto
end Spin.Structured.ConcreteNativeFamily
