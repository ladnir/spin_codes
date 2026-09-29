import SpinCodes.Numeric.BAExponent
import SpinCodes.Structured.RouteDomination

noncomputable section
namespace Spin.Structured.ConcreteOuter
open Finset

/-- The unnormalized binary entropy of a finite weight layer, including n=0. -/
def layerEntropy (n k : ℕ) : ℝ := (n : ℝ) * Spin.Numeric.hEnt ((k : ℝ) / n)

@[simp] theorem layerEntropy_zero (n : ℕ) : layerEntropy n 0 = 0 := by simp [layerEntropy]
@[simp] theorem layerEntropy_self (n : ℕ) : layerEntropy n n = 0 := by
  by_cases hn : n = 0
  · simp [hn, layerEntropy]
  · simp [layerEntropy, hn]

/-- The probability of one word at its empirical Bernoulli parameter is exp(-entropy). -/
theorem type_power_eq_exp {n k : ℕ} (hk : k ≤ n) :
    ((k : ℝ) / n)^k * (1 - (k : ℝ) / n)^(n-k) = Real.exp (-layerEntropy n k) := by
  by_cases hk0 : k = 0
  · subst k
    simp
  by_cases hkn : k = n
  · subst k
    have hn : (n : ℝ) ≠ 0 := by exact_mod_cast hk0
    simp [hn]
  have hn : 0 < n := by omega
  have hnR : (0 : ℝ) < n := by exact_mod_cast hn
  have hkR : (0 : ℝ) < k := by exact_mod_cast Nat.pos_of_ne_zero hk0
  have hp : (0 : ℝ) < (k : ℝ) / n := by positivity
  have hp1 : (k : ℝ) / n < 1 := (div_lt_one hnR).mpr (by exact_mod_cast (show k < n by omega))
  have hc : 0 < 1 - (k : ℝ) / n := by linarith
  have ht : 0 < ((k : ℝ) / n)^k * (1 - (k : ℝ) / n)^(n-k) := by positivity
  have hh : Real.log (((k : ℝ) / n)^k * (1 - (k : ℝ) / n)^(n-k)) = -layerEntropy n k := by
    rw [Real.log_mul (by positivity) (by positivity), Real.log_pow, Real.log_pow]
    unfold layerEntropy Spin.Numeric.hEnt
    rw [Nat.cast_sub hk]
    field_simp
    ring
  rw [← hh, Real.exp_log ht]

/-- Uniform upper type bound, valid at both boundary weights and zero length. -/
theorem choose_le_exp_entropy {n k : ℕ} (hk : k ≤ n) :
    (n.choose k : ℝ) ≤ Real.exp (layerEntropy n k) := by
  by_cases hn : n = 0
  · subst n
    have : k = 0 := by omega
    subst k
    simp
  have hnR : (0 : ℝ) < n := by exact_mod_cast Nat.pos_of_ne_zero hn
  have h0 : (0 : ℝ) ≤ (k : ℝ) / n := by positivity
  have h1 : (k : ℝ) / n ≤ 1 := (div_le_one hnR).mpr (by exact_mod_cast hk)
  let P := binPMF n ((k : ℝ) / n) h0 h1
  have hm : P.p ⟨k, by omega⟩ ≤ 1 := by
    calc P.p ⟨k, by omega⟩ ≤ ∑ j, P.p j := Finset.single_le_sum (fun j _ => P.nonneg j) (mem_univ _)
      _ = 1 := P.total
  change (n.choose k : ℝ) * ((k : ℝ) / n)^k * (1 - (k : ℝ) / n)^(n-k) ≤ 1 at hm
  rw [mul_assoc, type_power_eq_exp hk, Real.exp_neg, ← div_eq_mul_inv] at hm
  exact (div_le_one (Real.exp_pos _)).mp hm

/-- Uniform lower type bound with the explicit n+1 factor; no Stirling remainder is hidden. -/
theorem exp_entropy_div_le_choose {n k : ℕ} (hk : k ≤ n) :
    Real.exp (layerEntropy n k) / ((n : ℝ) + 1) ≤ (n.choose k : ℝ) := by
  by_cases hn : n = 0
  · subst n
    have : k = 0 := by omega
    subst k
    simp
  have hnN : 0 < n := Nat.pos_of_ne_zero hn
  have hnR : (0 : ℝ) < n := by exact_mod_cast hnN
  have h0 : (0 : ℝ) ≤ (k : ℝ) / n := by positivity
  have h1 : (k : ℝ) / n ≤ 1 := (div_le_one hnR).mpr (by exact_mod_cast hk)
  have hm := types_lower_bound_at hnN hk h0 h1
  rw [mul_assoc, type_power_eq_exp hk, Real.exp_neg, ← div_eq_mul_inv] at hm
  have hh := (le_div_iff₀ (Real.exp_pos _)).mp hm
  simpa only [one_div_mul_eq_div] using hh

end Spin.Structured.ConcreteOuter
