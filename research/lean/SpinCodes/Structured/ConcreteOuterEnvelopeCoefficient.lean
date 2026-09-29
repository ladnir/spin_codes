import SpinCodes.Structured.ConcreteOuterEntropy
import SpinCodes.Structured.BAGolay

noncomputable section
namespace Spin.Structured.ConcreteOuter
open Finset Spin.Numeric

/-- A nonnegative coefficient term is at most the polynomial's positive evaluation. -/
theorem coeff_mul_pow_le_eval (P : Polynomial ℕ) (a : ℕ) {u : ℝ} (hu : 0 ≤ u) :
    (P.coeff a : ℝ) * u^a ≤ P.eval₂ (Nat.castRingHom ℝ) u := by
  rw [Polynomial.eval₂_eq_sum, Polynomial.sum]
  change (P.coeff a : ℝ) * u^a ≤ ∑ i ∈ P.support, (P.coeff i : ℝ) * u^i
  by_cases ha : a ∈ P.support
  · exact Finset.single_le_sum (f := fun i => (P.coeff i : ℝ) * u^i)
      (fun i _ => mul_nonneg (Nat.cast_nonneg _) (pow_nonneg hu _)) ha
  · rw [Polynomial.notMem_support_iff.mp ha, Nat.cast_zero, zero_mul]
    exact Finset.sum_nonneg (fun i _ => by positivity)

def golayPolynomial : Polynomial ℕ := enumerator Golay.W

def golayCoefficient (k a : ℕ) : ℕ := (golayPolynomial ^ k).coeff a

theorem golayPolynomial_eval (u : ℝ) :
    golayPolynomial.eval₂ (Nat.castRingHom ℝ) u = GolayG u := by
  rw [GolayG_eq_sum]
  simp only [golayPolynomial, enumerator, Polynomial.eval₂_finsetSum,
    Polynomial.eval₂_pow, Polynomial.eval₂_X]

theorem golayCoefficient_mul_pow_le (k a : ℕ) {u : ℝ} (hu : 0 ≤ u) :
    (golayCoefficient k a : ℝ) * u^a ≤ GolayG u ^ k := by
  have hh := coeff_mul_pow_le_eval (golayPolynomial^k) a hu
  simpa only [golayCoefficient, Polynomial.eval₂_pow, golayPolynomial_eval] using hh

theorem golayCoefficient_log_le_objective {k a : ℕ} (hk : 0 < k)
    (hc : 0 < golayCoefficient k a) {u : ℝ} (hu : 0 < u) :
    Real.log (golayCoefficient k a : ℝ) / ((k * 24 : ℕ) : ℝ) ≤
      gObj u ((a : ℝ) / (k * 24 : ℕ)) := by
  have hcR : (0 : ℝ) < golayCoefficient k a := by exact_mod_cast hc
  have hbR : (0 : ℝ) < (k * 24 : ℕ) := by exact_mod_cast (by omega : 0 < k * 24)
  have hG := GolayG_pos hu
  have hh := Real.log_le_log (by positivity : (0 : ℝ) < (golayCoefficient k a : ℝ) * u^a)
    (golayCoefficient_mul_pow_le k a hu.le)
  rw [Real.log_mul hcR.ne' (by positivity), Real.log_pow, Real.log_pow] at hh
  apply (div_le_iff₀ hbR).mpr
  have he : gObj u ((a : ℝ) / (k * 24 : ℕ)) * ((k * 24 : ℕ) : ℝ) =
      (k : ℝ) * Real.log (GolayG u) - (a : ℝ) * Real.log u := by
    unfold gObj
    push_cast
    field_simp
  rw [he]
  linarith

/-- Exact finite coefficient bound by the paper's infimum exponent, with no asymptotic loss. -/
theorem golayCoefficient_le_exp_gBA {k : ℕ} (hk : 0 < k) (a : ℕ) :
    (golayCoefficient k a : ℝ) ≤
      Real.exp (((k * 24 : ℕ) : ℝ) * gBA ((a : ℝ) / (k * 24 : ℕ))) := by
  by_cases hc : golayCoefficient k a = 0
  · rw [hc, Nat.cast_zero]
    exact (Real.exp_pos _).le
  have hcp : 0 < golayCoefficient k a := Nat.pos_of_ne_zero hc
  have hcR : (0 : ℝ) < golayCoefficient k a := by exact_mod_cast hcp
  have hbR : (0 : ℝ) < (k * 24 : ℕ) := by exact_mod_cast (by omega : 0 < k * 24)
  have hh : Real.log (golayCoefficient k a : ℝ) / ((k * 24 : ℕ) : ℝ) ≤
      gBA ((a : ℝ) / (k * 24 : ℕ)) := by
    unfold gBA
    refine le_csInf ?_ ?_
    · exact ⟨gObj 1 ((a : ℝ) / (k * 24 : ℕ)), ⟨1, by norm_num, rfl⟩⟩
    · rintro _ ⟨u, hu, rfl⟩
      exact golayCoefficient_log_le_objective hk hcp hu
  rw [← Real.exp_log hcR]
  apply Real.exp_le_exp.mpr
  have he := (div_le_iff₀ hbR).mp hh
  simpa only [mul_comm] using he

end Spin.Structured.ConcreteOuter


