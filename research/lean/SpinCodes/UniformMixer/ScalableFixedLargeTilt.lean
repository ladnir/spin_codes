import SpinCodes.UniformMixer.ScalableFixedLargeProfileCount
import SpinCodes.Structured.ConcreteNativeFixedLargeLimit

/-! Charge the small Scalable tilt excess in the fixed-row exponent. -/

noncomputable section

namespace Spin.UniformMixer

open Spin.Structured Spin.Structured.ConcreteNativeFamily
  Spin.Structured.Placement

theorem tilt_div_identity (A θ θOld L : ℝ) (d : ℕ) :
    A / Real.exp (-θ / L) ^ d =
      (A / Real.exp (-θOld / L) ^ d) *
        Real.exp (((θ - θOld) / L) * (d : ℝ)) := by
  have he (t : ℝ) :
      A / Real.exp (-t / L) ^ d =
        A * Real.exp ((t / L) * (d : ℝ)) := by
    rw [← Real.exp_nat_mul]
    calc
      A / Real.exp ((d : ℝ) * (-t / L)) =
          A * Real.exp (-((d : ℝ) * (-t / L))) := by
        simp only [Real.exp_neg, div_eq_mul_inv]
      _ = _ := by congr 1; ring_nf
  rw [he θ, he θOld, mul_assoc, ← Real.exp_add]
  congr 1
  ring_nf

def scalableTiltCorrection (Q m : ℕ) : ℝ :=
  Real.exp (((scalableLargeTheta Q - paperLargeTheta Q) /
    (Lsched m : ℝ)) * (threshold m : ℝ))

theorem scalableTiltCorrection_le (Q m : ℕ) :
    scalableTiltCorrection Q m ≤
      Real.exp ((Q : ℝ) * (bsched m : ℝ) / 2000000) := by
  have hL : (0 : ℝ) < Lsched m := by
    exact_mod_cast Lsched_pos m
  have hQ : (0 : ℝ) ≤ (Q : ℝ) := Nat.cast_nonneg Q
  have hb : (0 : ℝ) ≤ (bsched m : ℝ) := Nat.cast_nonneg _
  have hθ := scalableLargeTheta_excess Q
  have hd := threshold_le m
  have ht :
      ((scalableLargeTheta Q - paperLargeTheta Q) /
        (Lsched m : ℝ)) * (threshold m : ℝ) ≤
        (Q : ℝ) * (bsched m : ℝ) / 2000000 := by
    have hdiv : 0 ≤ (scalableLargeTheta Q -
        paperLargeTheta Q) / (Lsched m : ℝ) :=
      div_nonneg hθ.1 hL.le
    have h1 := mul_le_mul_of_nonneg_left hd hdiv
    have h2 := div_le_div_of_nonneg_right hθ.2 hL.le
    have h3 := mul_le_mul_of_nonneg_right h2
      (show (0 : ℝ) ≤ 11 / 100 *
        (Lsched m : ℝ) * (bsched m : ℝ) by positivity)
    have hc :
        ((Q : ℝ) / 250000 / (Lsched m : ℝ)) *
          (11 / 100 * (Lsched m : ℝ) * (bsched m : ℝ)) =
        (11 / 25000000 : ℝ) * (Q : ℝ) * (bsched m : ℝ) := by
      field_simp
      ring
    rw [hc] at h3
    nlinarith only [h1, h3, mul_nonneg hQ hb]
  exact Real.exp_le_exp.mpr ht

end Spin.UniformMixer
