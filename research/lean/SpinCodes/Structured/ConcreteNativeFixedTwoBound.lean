import SpinCodes.Structured.ConcreteFixedRoutedTwo
import SpinCodes.Structured.ConcreteNativeFixedOneBound
import Mathlib.Data.Nat.Choose.Bounds

noncomputable section
namespace Spin.Structured.ConcreteNativeFamily
open Finset Filter ConcreteOuter ConcreteRoute ConcreteMarked Placement

def fixedTwoBound (m : ℕ) : ℝ :=
  (Lsched m:ℝ)^2*Real.exp (2*(bsched m:ℝ)*(Real.log 2/2+1281/100000+
    selectedRemainder m (nativeShellRemainder m)))*
    (1000*(1011/4000:ℝ)^(bsched m))/(Real.exp (-(4/(Lsched m:ℝ))))^(threshold m)

/-- Actual selected native two-row message layer, with no kernel inequality left as a premise. -/
theorem native_two_tuple_bound :
    ∀ᶠ m in atTop, ∀ seed, nativeGood m seed →
      (∑ x ∈ univ.filter (fun x : NativeMessage m => occupation x=2),
        failureProbability (tupleWiring m) (rowSupports seed x) (threshold m)) ≤ fixedTwoBound m := by
  have hlim : Tendsto (fun m : ℕ => m+1) atTop atTop := tendsto_add_atTop_nat 1
  filter_upwards [hlim.eventually actual_two_fair_probability_regionMajor] with m hm
  intro seed hg
  have hB := native_good_spectrum_bound m seed hg
    (B:=spectrumRatioBound (bsched m) (nativeShellRemainder m)) (by unfold spectrumRatioBound; positivity)
    (native_expected_spectrum_uniform m)
  have hp : ∀ S : Finset (Fin (Lsched m)), S.card=2 →
      (fairExperimentLaw S (nativeBlocks m*24) (rounds m)).prob
        (fun ω => ConcreteRoutedEncoder.weight (tupleWiring m) ω.1 ω.2≤threshold m) ≤
          (1000*(1011/4000:ℝ)^(bsched m))/(Real.exp (-(4/(Lsched m:ℝ))))^(threshold m) := by
    intro S hS
    rw [tupleWiring_eq_regionMajor]
    have hh := hm (nativeBlocks m*24) (rounds m)
      (by change (nativeBlocks m*24)*Lsched m=rounds m*128
          rw [native_width]; simpa only [Nat.mul_comm] using (rounds_length m).symm) S hS (threshold m)
    simp only [native_width,Lsched,Nat.cast_mul,Nat.cast_ofNat,Nat.cast_add,Nat.cast_one] at hh ⊢
    convert hh using 1
    apply FinPMF.prob_congr
    intro ω
    rfl
  have hh := occupation_failure_le_fair seed 2 (tupleWiring m) (threshold m)
    (by simpa only [native_width] using hB) hp
  have hc : ((Lsched m).choose 2:ℝ) ≤ (Lsched m:ℝ)^2 := by exact_mod_cast Nat.choose_le_pow (Lsched m) 2
  have hh' := hh.trans (mul_le_mul_of_nonneg_right hc (by positivity))
  simpa only [native_width,selected_row_cost,fixedTwoBound,←Real.exp_nat_mul,
    Nat.cast_ofNat,mul_div_assoc,mul_assoc] using hh' 

theorem native_two_EZ_bound : ∀ᶠ m in atTop, concreteFamily.EZ m 2 ≤ fixedTwoBound m := by
  filter_upwards [native_good_positive_eventually,native_two_tuple_bound] with m hm hb
  exact EZ_le_tuple_layer_bound m 2 hm hb

#print axioms native_two_tuple_bound
#print axioms native_two_EZ_bound
end Spin.Structured.ConcreteNativeFamily
