import SpinCodes.Turbo.PhysicalProfileNumericRpowDefs22
import SpinCodes.Turbo.NumericRpow
import SpinCodes.Numeric.Fixed

namespace Spin.Turbo
open Spin.Numeric Spin.Numeric.Fix

theorem integerRpowUpperCheck22_sound (p q : Nat) (an ad un ud : Int)
    (h : integerRpowUpperCheck22 p q an ad un ud = true) :
    ((an : ℝ)/(ad : ℝ))^((p : ℝ)/(q : ℝ)) ≤ (un : ℝ)/(ud : ℝ) := by
  have hc : 0 < q ∧ 0 ≤ an ∧ 0 < ad ∧ 0 ≤ un ∧ 0 < ud ∧
      an^p * ud^q ≤ un^q * ad^p := of_decide_eq_true h
  have ha : (0 : ℝ) ≤ (an : ℝ)/(ad : ℝ) :=
    div_nonneg (by exact_mod_cast hc.2.1) (by exact_mod_cast hc.2.2.1.le)
  have hu : (0 : ℝ) ≤ (un : ℝ)/(ud : ℝ) :=
    div_nonneg (by exact_mod_cast hc.2.2.2.1) (by exact_mod_cast hc.2.2.2.2.1.le)
  apply rpow_nat_div_le_of_pow_le ha hu p q hc.1
  rw [div_pow, div_pow]
  apply (div_le_div_iff₀ (pow_pos (by exact_mod_cast hc.2.2.1) p)
    (pow_pos (by exact_mod_cast hc.2.2.2.2.1) q)).mpr
  exact_mod_cast hc.2.2.2.2.2

theorem profileFixRpowUpperCheck22_sound {box : Fix} {x : ℝ}
    (hx : Mem box x) (hn : 0 ≤ x) (p q : Nat) (num den : Int)
    (h : profileFixRpowUpperCheck22 box p q num den = true) :
    x^((p : ℝ)/(q : ℝ)) ≤ (num : ℝ)/(den : ℝ) := by
  have hc : 0 < q ∧ 0 ≤ box.hi ∧ 0 < Spin.Numeric.scale ∧ 0 ≤ num ∧ 0 < den ∧
      box.hi^p * den^q ≤ num^q * Spin.Numeric.scale^p := of_decide_eq_true h
  exact (Real.rpow_le_rpow hn hx.2 (div_nonneg (Nat.cast_nonneg p) (Nat.cast_nonneg q))).trans
    (integerRpowUpperCheck22_sound p q box.hi Spin.Numeric.scale num den h)

#print axioms integerRpowUpperCheck22_sound
#print axioms profileFixRpowUpperCheck22_sound
end Spin.Turbo
