import SpinCodes.Structured.VarianceBound

/-! Integer certificates for upper bounds derived from the variance
inequality. Checking the quadratic at `cap + 1` avoids formalizing a square
root algorithm and certifies the same rounded upper bound. -/

namespace Spin

theorem variance_cap_of_certificate {M a S t cap : ℕ}
    (hv : M * t ^ 2 + a ^ 2 ≤ 2 * a * t + (M - 1) * S)
    (hm : a ≤ M * (cap + 1))
    (hc : 2 * a * (cap + 1) + (M - 1) * S < M * (cap + 1) ^ 2 + a ^ 2) :
    t ≤ cap := by
  by_contra ht
  have ht' : (cap : ℤ) + 1 ≤ t := by exact_mod_cast (show cap + 1 ≤ t by omega)
  have hm' : (a : ℤ) ≤ M * (cap + 1) := by exact_mod_cast hm
  have hv' : (M : ℤ) * t ^ 2 + a ^ 2 ≤ 2 * a * t + ((M - 1 : ℕ) : ℤ) * S := by
    exact_mod_cast hv
  have hc' : 2 * (a : ℤ) * (cap + 1) + ((M - 1 : ℕ) : ℤ) * S <
      M * (cap + 1) ^ 2 + a ^ 2 := by exact_mod_cast hc
  have hMt : (M : ℤ) * (cap + 1) ≤ M * t := mul_le_mul_of_nonneg_left ht' (by positivity)
  have hsecond : 0 ≤ (M : ℤ) * (t + (cap + 1)) - 2 * a := by nlinarith
  have hprod := mul_nonneg (sub_nonneg.mpr ht') hsecond
  nlinarith

end Spin
