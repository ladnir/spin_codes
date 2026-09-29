import SpinCodes.Structured.ConcreteOuterMajorantBridge

noncomputable section
namespace Spin.Structured.ConcreteOuter
open Spin.Numeric

/-- Uniform explicit finite-size remainder for the expected spectrum. -/
def spectrumLogError (b : ℕ) : ℝ := 16*Real.log ((b:ℝ)+1)+13
/-- The additional log(b+1) is exactly the lower type-bound denominator. -/
def shellLogError (b : ℕ) : ℝ := 17*Real.log ((b:ℝ)+1)+13

theorem expected_spectrum_refined {k w : ℕ} (hk : 0 < k)
    (hw : (w:ℝ)/(k*24:ℕ) ∈ Set.Icc ((13:ℝ)/125) (112/125)) :
    (seedLaw k).expect (fun seed => (spectrum seed w:ℝ)) ≤
      Real.exp (((k*24:ℕ):ℝ)*Spin.Majorant.refined.toFun ((w:ℝ)/(k*24:ℕ)) +
        spectrumLogError (k*24)) := by
  refine (expected_spectrum_explicit_log_bound hk w).trans (Real.exp_le_exp.mpr ?_)
  have hb : (0:ℝ)<(k*24:ℕ) := by exact_mod_cast (by omega : 0<k*24)
  have hl := Real.log_le_log hb (show ((k*24:ℕ):ℝ) ≤ ((k*24:ℕ):ℝ)+1 by linarith)
  have hh := finiteEnvelope_le_refined hk hw
  unfold spectrumLogError
  linarith

/-- The actual expected shell ratio is controlled by the certified majorant
minus binary entropy, with no omitted parity or boundary remainder. -/
theorem expected_spectrum_shell_refined {k w : ℕ} (hk : 0 < k)
    (hw : (w:ℝ)/(k*24:ℕ) ∈ Set.Icc ((13:ℝ)/125) (112/125)) :
    (seedLaw k).expect (fun seed => (spectrum seed w:ℝ)) ≤
      Real.exp (((k*24:ℕ):ℝ)*(Spin.Majorant.refined.toFun ((w:ℝ)/(k*24:ℕ)) -
        hEnt ((w:ℝ)/(k*24:ℕ))) + shellLogError (k*24)) * ((k*24).choose w:ℝ) := by
  have hb : (0:ℝ)<(k*24:ℕ) := by exact_mod_cast (by omega : 0<k*24)
  have hwb : w ≤ k*24 := by
    have hh : (w:ℝ)/(k*24:ℕ) ≤ 1 := by linarith [hw.2]
    exact_mod_cast (div_le_one hb).mp hh
  have ht := exp_entropy_div_le_choose hwb
  let X := ((k*24:ℕ):ℝ)*Spin.Majorant.refined.toFun ((w:ℝ)/(k*24:ℕ)) + spectrumLogError (k*24)
  let H := layerEntropy (k*24) w
  let t : ℝ := ((k*24:ℕ):ℝ)+1
  have htpos : 0<t := by dsimp [t]; positivity
  have he : Real.exp (X-H+Real.log t) * (Real.exp H/t) = Real.exp X := by
    rw [Real.exp_add, Real.exp_sub, Real.exp_log htpos]
    field_simp
  have hbound : Real.exp X ≤ Real.exp (X-H+Real.log t)*((k*24).choose w:ℝ) := by
    rw [← he]
    exact mul_le_mul_of_nonneg_left ht (Real.exp_pos _).le
  have he' : X-H+Real.log t = ((k*24:ℕ):ℝ)*(Spin.Majorant.refined.toFun ((w:ℝ)/(k*24:ℕ)) -
      hEnt ((w:ℝ)/(k*24:ℕ)))+shellLogError (k*24) := by
    dsimp [X,H,t,layerEntropy,spectrumLogError,shellLogError]
    ring
  rw [he'] at hbound
  exact (expected_spectrum_refined hk hw).trans hbound

theorem expected_shell_ratio_refined {k w : ℕ} (hk : 0 < k)
    (hw : (w:ℝ)/(k*24:ℕ) ∈ Set.Icc ((13:ℝ)/125) (112/125)) :
    (seedLaw k).expect (fun seed => (spectrum seed w:ℝ)) / ((k*24).choose w:ℝ) ≤
      Real.exp (((k*24:ℕ):ℝ)*(Spin.Majorant.refined.toFun ((w:ℝ)/(k*24:ℕ)) -
        hEnt ((w:ℝ)/(k*24:ℕ))) + shellLogError (k*24)) := by
  have hb : (0:ℝ)<(k*24:ℕ) := by exact_mod_cast (by omega : 0<k*24)
  have hwb : w ≤ k*24 := by
    have hh : (w:ℝ)/(k*24:ℕ) ≤ 1 := by linarith [hw.2]
    exact_mod_cast (div_le_one hb).mp hh
  exact (div_le_iff₀ (by exact_mod_cast Nat.choose_pos hwb : (0:ℝ)<(k*24).choose w)).mpr
    (expected_spectrum_shell_refined hk hw)

end Spin.Structured.ConcreteOuter
