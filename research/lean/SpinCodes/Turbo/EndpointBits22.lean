import Mathlib

/-! Converting an exact integer certificate `T / 2^s` into a bit-security
bound `2^(-(m/n))`. The hypothesis is a single integer inequality
`T^n * 2^m ≤ 2^(s*n)`, checked by the kernel at the endpoint. -/
noncomputable section
namespace Spin.Turbo

theorem int_div_two_pow_le_two_rpow_neg (T : ℤ) (s m n : ℕ) (hn : n ≠ 0)
    (hT : 0 ≤ T) (h : T ^ n * 2 ^ m ≤ 2 ^ (s * n)) :
    (T : ℝ) / 2 ^ s ≤ (2 : ℝ) ^ (-((m : ℝ) / n)) := by
  have hx : 0 ≤ (T : ℝ) / 2 ^ s := div_nonneg (by exact_mod_cast hT) (by positivity)
  have hy : 0 ≤ (2 : ℝ) ^ (-((m : ℝ) / n)) := by positivity
  rw [← pow_le_pow_iff_left₀ hx hy hn]
  have hnR : (n : ℝ) ≠ 0 := by exact_mod_cast hn
  have hy' : ((2 : ℝ) ^ (-((m : ℝ) / n))) ^ n = ((2 : ℝ) ^ m)⁻¹ := by
    rw [← Real.rpow_mul_natCast (by norm_num : (0 : ℝ) ≤ 2)]
    rw [show -((m : ℝ) / n) * n = -(m : ℝ) by field_simp]
    rw [Real.rpow_neg (by norm_num), Real.rpow_natCast]
  rw [hy', div_pow, ← pow_mul]
  have h' : (T : ℝ) ^ n * 2 ^ m ≤ 2 ^ (s * n) := by exact_mod_cast h
  rw [div_le_iff₀ (by positivity), inv_mul_eq_div, le_div_iff₀ (by positivity)]
  exact h'

/-- Same conclusion from a natural-number check (GMP `Nat.pow` in the kernel,
so large `n` such as `10^5` stays cheap). -/
theorem nat_div_two_pow_le_two_rpow_neg (T : ℤ) (s m n : ℕ) (hn : n ≠ 0)
    (hT : 0 ≤ T) (h : T.toNat ^ n * 2 ^ m ≤ 2 ^ (s * n)) :
    (T : ℝ) / 2 ^ s ≤ (2 : ℝ) ^ (-((m : ℝ) / n)) := by
  apply int_div_two_pow_le_two_rpow_neg T s m n hn hT
  have h' : ((T.toNat ^ n * 2 ^ m : ℕ) : ℤ) ≤ ((2 ^ (s * n) : ℕ) : ℤ) := by exact_mod_cast h
  push_cast at h'
  rwa [Int.toNat_of_nonneg hT] at h'

#print axioms int_div_two_pow_le_two_rpow_neg
#print axioms nat_div_two_pow_le_two_rpow_neg
end Spin.Turbo
