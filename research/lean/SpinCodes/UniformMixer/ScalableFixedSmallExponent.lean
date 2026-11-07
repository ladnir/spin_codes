import SpinCodes.UniformMixer.ScalableFixedFairSmall
import SpinCodes.Structured.ConcreteNativeFixedOneLimit
import SpinCodes.Structured.ConcreteNativeFixedTwoLimit

/-! The scalable one- and two-row tilts retain the native outer exponent gap. -/

noncomputable section

namespace Spin.UniformMixer

open Spin.Structured Spin.Structured.ConcreteNativeFamily Filter

theorem scalableThetaOne_le_two : scalableThetaOne ≤ 2 := by
  norm_num [scalableThetaOne]

theorem scalableThetaTwo_le_four : scalableThetaTwo ≤ 4 := by
  norm_num [scalableThetaTwo, scalableThetaOne]

theorem scalable_fixed_one_log_gap :
    4 * Real.log 2 / 39 + Real.log 2 / 2 + 1281 / 100000 +
      Real.log (1011 / 2000) + (11 / 100) * scalableThetaOne ≤
        -(3 / 100 : ℝ) := by
  have h := Spin.Structured.ConcreteNativeFamily.fixed_one_log_gap
  have ht := mul_le_mul_of_nonneg_left scalableThetaOne_le_two
    (by norm_num : (0 : ℝ) ≤ 11 / 100)
  nlinarith only [h, ht]

theorem scalable_fixed_two_log_gap :
    8 * Real.log 2 / 39 +
      2 * (Real.log 2 / 2 + 1281 / 100000) +
      Real.log (1011 / 4000) +
      (11 / 100) * scalableThetaTwo ≤ -(7 / 100 : ℝ) := by
  have h := Spin.Structured.ConcreteNativeFamily.fixed_two_log_gap
  have ht := mul_le_mul_of_nonneg_left scalableThetaTwo_le_four
    (by norm_num : (0 : ℝ) ≤ 11 / 100)
  nlinarith only [h, ht]

def scalableFixedOneBound (m : ℕ) : ℝ :=
  (Lsched m : ℝ) *
    Real.exp ((bsched m : ℝ) *
      (Real.log 2 / 2 + 1281 / 100000 +
        selectedRemainder m (nativeShellRemainder m))) *
    (1000 * (1011 / 2000 : ℝ) ^ bsched m) /
      Real.exp (-(scalableThetaOne / (Lsched m : ℝ))) ^
        threshold m

def scalableFixedTwoBound (m : ℕ) : ℝ :=
  (Lsched m : ℝ) ^ 2 *
    Real.exp (2 * (bsched m : ℝ) *
      (Real.log 2 / 2 + 1281 / 100000 +
        selectedRemainder m (nativeShellRemainder m))) *
    (1000 * (1011 / 4000 : ℝ) ^ bsched m) /
      Real.exp (-(scalableThetaTwo / (Lsched m : ℝ))) ^
        threshold m

theorem scalableFixedOneBound_le_native (m : ℕ) :
    scalableFixedOneBound m ≤ fixedOneBound m := by
  have hL : (0 : ℝ) < Lsched m := by
    exact_mod_cast Lsched_pos m
  have hθ : scalableThetaOne / (Lsched m : ℝ) ≤
      2 / (Lsched m : ℝ) :=
    div_le_div_of_nonneg_right scalableThetaOne_le_two hL.le
  have hz : Real.exp (-(2 / (Lsched m : ℝ))) ≤
      Real.exp (-(scalableThetaOne / (Lsched m : ℝ))) :=
    Real.exp_le_exp.mpr (by linarith)
  have hp := pow_le_pow_left₀ (Real.exp_pos _).le hz
    (threshold m)
  unfold scalableFixedOneBound fixedOneBound
  exact div_le_div_of_nonneg_left (by positivity)
    (pow_pos (Real.exp_pos _) _) hp

theorem scalableFixedTwoBound_le_native (m : ℕ) :
    scalableFixedTwoBound m ≤ fixedTwoBound m := by
  have hL : (0 : ℝ) < Lsched m := by
    exact_mod_cast Lsched_pos m
  have hθ : scalableThetaTwo / (Lsched m : ℝ) ≤
      4 / (Lsched m : ℝ) :=
    div_le_div_of_nonneg_right scalableThetaTwo_le_four hL.le
  have hz : Real.exp (-(4 / (Lsched m : ℝ))) ≤
      Real.exp (-(scalableThetaTwo / (Lsched m : ℝ))) :=
    Real.exp_le_exp.mpr (by linarith)
  have hp := pow_le_pow_left₀ (Real.exp_pos _).le hz
    (threshold m)
  unfold scalableFixedTwoBound fixedTwoBound
  exact div_le_div_of_nonneg_left (by positivity)
    (pow_pos (Real.exp_pos _) _) hp

theorem scalableFixedOneBound_eventually :
    ∀ᶠ m : ℕ in atTop,
      scalableFixedOneBound m ≤
        1000 * Real.exp (-(1 / 50 : ℝ) * (bsched m : ℝ)) := by
  filter_upwards [fixedOneBound_eventually] with m hm
  exact (scalableFixedOneBound_le_native m).trans hm

theorem scalableFixedTwoBound_eventually :
    ∀ᶠ m : ℕ in atTop,
      scalableFixedTwoBound m ≤
        1000 * Real.exp (-(1 / 20 : ℝ) * (bsched m : ℝ)) := by
  filter_upwards [fixedTwoBound_eventually] with m hm
  exact (scalableFixedTwoBound_le_native m).trans hm

end Spin.UniformMixer
