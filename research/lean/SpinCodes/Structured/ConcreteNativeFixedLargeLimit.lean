import SpinCodes.Structured.ConcreteNativeFixedLargeSpectrum
import SpinCodes.Structured.DenseNativeRemainder

noncomputable section
namespace Spin.Structured.ConcreteNativeFamily
open Filter ConcreteOuter

def fixedLargeRemainder (m : ℕ) : ℝ :=
  selectedRemainder m (nativeShellRemainder m)+Real.log ((bsched m:ℝ)+1)/(bsched m:ℝ)

def fixedLargeBound (Q m : ℕ) : ℝ :=
  (1600/3)*Real.exp ((Q:ℝ)*(Real.log ((bsched m:ℝ)+1)+2*Real.log (bsched m:ℝ)+shellLogError (bsched m))-
    (7679/10000000)*(Q:ℝ)*(bsched m:ℝ))

theorem fixedLargeRemainder_tendsto : Tendsto fixedLargeRemainder atTop (nhds 0) := by
  have h := (selectedRemainder_tendsto nativeShellRemainder_tendsto).add
    (DenseOccupationFixed.log_succ_div_tendsto.comp bsched_tendsto)
  change Tendsto (fun m => selectedRemainder m (nativeShellRemainder m)+Real.log ((bsched m:ℝ)+1)/(bsched m:ℝ)) atTop (nhds 0)
  simpa only [Function.comp_apply,zero_add] using h

theorem fixedLargeBound_exp (Q m : ℕ) : fixedLargeBound Q m=
    (1600/3)*Real.exp ((Q:ℝ)*(bsched m:ℝ)*(fixedLargeRemainder m-7679/10000000)) := by
  have hb : (bsched m:ℝ)≠0 := by exact_mod_cast (bsched_pos m).ne'
  unfold fixedLargeBound fixedLargeRemainder selectedRemainder nativeShellRemainder
  congr 2
  field_simp
  ring

theorem fixedLargeBound_eventually (Q : ℕ) : ∀ᶠ m in atTop,
    fixedLargeBound Q m≤(1600/3)*Real.exp (-(1/2000:ℝ)*(Q:ℝ)*(bsched m:ℝ)) := by
  filter_upwards [fixedLargeRemainder_tendsto.eventually (gt_mem_nhds (by norm_num : (0:ℝ)<1/5000))] with m hm
  rw [fixedLargeBound_exp]
  apply mul_le_mul_of_nonneg_left _ (by norm_num)
  apply Real.exp_le_exp.mpr
  have h := mul_le_mul_of_nonneg_left (show fixedLargeRemainder m-7679/10000000≤-(1/2000:ℝ) by linarith) (show 0≤(Q:ℝ)*(bsched m:ℝ) by positivity)
  nlinarith

theorem fixedLargeBound_tendsto {Q : ℕ} (hQ : 0<Q) : Tendsto (fixedLargeBound Q) atTop (nhds 0) := by
  have hb : Tendsto (fun m => (bsched m:ℝ)) atTop atTop := tendsto_natCast_atTop_atTop.comp bsched_tendsto
  have he := (Real.tendsto_exp_atBot.comp (hb.const_mul_atTop_of_neg
    (show -(1/2000:ℝ)*(Q:ℝ)<0 by
      have hQr : (0:ℝ)<Q := by exact_mod_cast hQ
      nlinarith))).const_mul (1600/3)
  simp only [mul_zero] at he
  apply tendsto_of_tendsto_of_tendsto_of_le_of_le' tendsto_const_nhds he
  · exact Eventually.of_forall (fun m => by unfold fixedLargeBound; positivity)
  · exact fixedLargeBound_eventually Q

#print axioms fixedLargeRemainder_tendsto
#print axioms fixedLargeBound_eventually
#print axioms fixedLargeBound_tendsto
end Spin.Structured.ConcreteNativeFamily
