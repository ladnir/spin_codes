import SpinCodes.Structured.ConcreteNativeFixedLargeDefs

noncomputable section
namespace Spin.Structured.ConcreteNativeFamily
open Finset ConcreteOuter ConcreteFixedNumeric Spin.Numeric

def fixedShellCost (b : ℕ) : ℝ := (b:ℝ)^2*Real.exp (shellLogError b)

def fixedRowCharge : ℝ := gThree+(11/100)*(133/125)/(262144/524287)+
  (4/39)*Real.log 2+8679/10000000

theorem selected_shell_ratio {k w : ℕ} (hk : 0 < k) (seed : Seed k) (W : Finset ℕ)
    (hg : Spin.Good (seedLaw k) spectrum (k*24) W seed) (hw : w∈W)
    (hr : (w:ℝ)/(k*24:ℕ)∈Set.Icc (13/125) (112/125)) :
    (spectrum seed w:ℝ)/((k*24).choose w:ℝ) ≤
      fixedShellCost (k*24)*Real.exp (((k*24:ℕ):ℝ)*
        (Spin.Majorant.refined.toFun ((w:ℝ)/(k*24:ℕ))-hEnt ((w:ℝ)/(k*24:ℕ)))) := by
  have hb : (0:ℝ)<(k*24:ℕ) := by positivity
  have hwb : w≤k*24 := by
    have hh : (w:ℝ)/(k*24:ℕ)≤1 := by linarith [hr.2]
    exact_mod_cast (div_le_one hb).mp hh
  have hc : (0:ℝ)<(k*24).choose w := by exact_mod_cast Nat.choose_pos hwb
  have hs := (hg.2 w hw).trans (mul_le_mul_of_nonneg_left (expected_spectrum_shell_refined hk hr) (sq_nonneg _))
  apply (div_le_iff₀ hc).mpr
  convert hs using 1
  unfold fixedShellCost
  rw [Real.exp_add]
  ring

theorem weighted_shell_bound {Q k w : ℕ} (hk : 0 < k) (seed : Seed k) (W : Finset ℕ)
    (hg : Spin.Good (seedLaw k) spectrum (k*24) W seed) (hw : w∈W)
    (hr : (w:ℝ)/(k*24:ℕ)∈Set.Icc (13/125) (112/125))
    (u : FixedFugacityChoice Q) (i : Fin Q)
    (hu : ConcreteFixedNumeric.rate ((w:ℝ)/(k*24:ℕ)) (fixedFugacity u i) ≤ -(8679/10000000)) :
    ((spectrum seed w:ℝ)/((k*24).choose w:ℝ))*
      ((Spin.mv (127/250) (3/1600) (1/524287) (fixedFugacity u i))^(k*24)/(fixedFugacity u i)^w) ≤
        fixedShellCost (k*24)*Real.exp (-((k*24:ℕ):ℝ)*fixedRowCharge) := by
  have hb : (0:ℝ)<(k*24:ℕ) := by positivity
  have hp : (Spin.mv (127/250) (3/1600) (1/524287) (fixedFugacity u i))^(k*24)/(fixedFugacity u i)^w =
      Real.exp (((k*24:ℕ):ℝ)*Real.log (Spin.mv (127/250) (3/1600) (1/524287) (fixedFugacity u i))-
        (w:ℝ)*Real.log (fixedFugacity u i)) := by
    rw [Real.exp_sub,Real.exp_nat_mul,Real.exp_nat_mul,Real.exp_log (fixedMv_pos u i),Real.exp_log (fixedFugacity_pos u i)]
  rw [hp]
  calc
    _ ≤ (fixedShellCost (k*24)*Real.exp (((k*24:ℕ):ℝ)*
        (Spin.Majorant.refined.toFun ((w:ℝ)/(k*24:ℕ))-hEnt ((w:ℝ)/(k*24:ℕ)))))*
        Real.exp (((k*24:ℕ):ℝ)*Real.log (Spin.mv (127/250) (3/1600) (1/524287) (fixedFugacity u i))-
          (w:ℝ)*Real.log (fixedFugacity u i)) :=
      mul_le_mul_of_nonneg_right (selected_shell_ratio hk seed W hg hw hr) (Real.exp_pos _).le
    _ ≤ _ := by
      rw [mul_assoc,←Real.exp_add]
      apply mul_le_mul_of_nonneg_left _ (by unfold fixedShellCost; positivity)
      apply Real.exp_le_exp.mpr
      have h := mul_le_mul_of_nonneg_left hu hb.le
      unfold ConcreteFixedNumeric.rate fixedRowCharge at *
      have hx : ((k*24:ℕ):ℝ)*((w:ℝ)/(k*24:ℕ))=(w:ℝ) := by field_simp
      have hxl := congrArg (fun t : ℝ => t*Real.log (fixedFugacity u i)) hx
      nlinarith

#print axioms selected_shell_ratio
#print axioms weighted_shell_bound
end Spin.Structured.ConcreteNativeFamily
