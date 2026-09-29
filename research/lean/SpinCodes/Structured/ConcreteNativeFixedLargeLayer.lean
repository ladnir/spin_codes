import SpinCodes.Structured.ConcreteNativeFixedLargeLimit
import SpinCodes.Structured.ConcreteNativeFixedLargeCount

noncomputable section
namespace Spin.Structured.ConcreteNativeFamily
open Finset ConcreteOuter ConcreteFixedNumeric

theorem fixedLarge_counting_cost_le (Q m : ℕ) :
    ((Lsched m).choose Q:ℝ)*(((bsched m:ℝ)+1)^Q*
      fixedProfileCountBound Q (bsched m) (Lsched m) (threshold m))≤fixedLargeBound Q m := by
  have hL : (0:ℝ)<Lsched m := by exact_mod_cast Lsched_pos m
  have hb : (0:ℝ)<bsched m := by exact_mod_cast bsched_pos m
  have hc : ((Lsched m).choose Q:ℝ)≤(Lsched m:ℝ)^Q := by exact_mod_cast Nat.choose_le_pow (Lsched m) Q
  have hLp : (Lsched m:ℝ)^Q=Real.exp ((Q:ℝ)*Real.log (Lsched m:ℝ)) := by
    rw [Real.exp_nat_mul,Real.exp_log hL]
  have hbp : ((bsched m:ℝ)+1)^Q=Real.exp ((Q:ℝ)*Real.log ((bsched m:ℝ)+1)) := by
    rw [Real.exp_nat_mul,Real.exp_log (by positivity)]
  have hsc : fixedShellCost (bsched m)^Q=Real.exp ((Q:ℝ)*(2*Real.log (bsched m:ℝ)+shellLogError (bsched m))) := by
    rw [Real.exp_nat_mul,Real.exp_add]
    unfold fixedShellCost
    rw [show Real.exp (2*Real.log (bsched m:ℝ))=(bsched m:ℝ)^2 by
      rw [show (2:ℝ)=((2:ℕ):ℝ) by norm_num,Real.exp_nat_mul,Real.exp_log hb]]
  have he : (Lsched m:ℝ)^Q*(((bsched m:ℝ)+1)^Q*fixedProfileCountBound Q (bsched m) (Lsched m) (threshold m))=
      (1600/3)*Real.exp ((Q:ℝ)*Real.log (Lsched m:ℝ)+(Q:ℝ)*Real.log ((bsched m:ℝ)+1)+
        (Q:ℝ)*(2*Real.log (bsched m:ℝ)+shellLogError (bsched m))-
        (bsched m:ℝ)*(Q:ℝ)*fixedRowCharge+(bsched m:ℝ)*(Q:ℝ)*(gThree+1/10000)+
        (((Q:ℝ)*(133/125)/(262144/524287))/(Lsched m:ℝ))*(threshold m:ℝ)) := by
    rw [hLp,hbp]
    unfold fixedProfileCountBound
    rw [hsc]
    simp only [Real.exp_add,Real.exp_sub]
    ring_nf
    simp only [Real.exp_neg]
    ring
  calc
    _ ≤ (Lsched m:ℝ)^Q*(((bsched m:ℝ)+1)^Q*fixedProfileCountBound Q (bsched m) (Lsched m) (threshold m)) := by
      apply mul_le_mul_of_nonneg_right hc
      unfold fixedProfileCountBound fixedShellCost
      positivity
    _ ≤ _ := by
      rw [he]
      unfold fixedLargeBound
      apply mul_le_mul_of_nonneg_left _ (by norm_num)
      apply Real.exp_le_exp.mpr
      have hsch := mul_le_mul_of_nonneg_left (SparseRate.native_log_schedule m) (Nat.cast_nonneg Q : (0:ℝ)≤Q)
      have hth := mul_le_mul_of_nonneg_left (threshold_le m)
        (show 0≤((Q:ℝ)*(133/125)/(262144/524287))/(Lsched m:ℝ) by positivity)
      have hcancel : (((Q:ℝ)*(133/125)/(262144/524287))/(Lsched m:ℝ))*
          ((11/100)*(Lsched m:ℝ)*(bsched m:ℝ))=
          (Q:ℝ)*((11/100)*(133/125)/(262144/524287))*(bsched m:ℝ) := by
        field_simp
      rw [hcancel] at hth
      unfold fixedRowCharge
      nlinarith

#print axioms fixedLarge_counting_cost_le
end Spin.Structured.ConcreteNativeFamily

