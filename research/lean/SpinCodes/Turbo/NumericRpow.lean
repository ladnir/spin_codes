import Mathlib.Analysis.SpecialFunctions.Pow.Real
import Mathlib.Data.Matrix.Mul

/-!
# Exact rational certificates for fractional-power upper endpoints

The old numeric checkers reduce rational/integer claims in the kernel and
then apply a separately proved real soundness theorem. This checker follows
that pattern: it verifies a^p ≤ u^q for q>0. It trusts neither floating
fractional powers nor an exported endpoint's provenance. Its input bound
must separately be connected to the actual physical or paired matrix entry.
-/
namespace Spin.Turbo

/-- Raising both nonnegative sides to a positive integer denominator removes
fractional powers from the proposed numerical inequality. -/
theorem rpow_nat_div_le_of_pow_le {a u : ℝ} (ha : 0 ≤ a) (hu : 0 ≤ u)
    (p q : ℕ) (hq : 0 < q) (hpow : a ^ p ≤ u ^ q) :
    a ^ ((p : ℝ) / (q : ℝ)) ≤ u := by
  have hqr : (0 : ℝ) < q := by exact_mod_cast hq
  apply (Real.rpow_le_rpow_iff (Real.rpow_nonneg ha _) hu hqr).mp
  rw [← Real.rpow_mul ha]
  rw [div_mul_cancel₀ _ hqr.ne', Real.rpow_natCast, Real.rpow_natCast]
  exact hpow

/-- Proposed rational base and upper endpoint, checked with exact arithmetic. -/
def rationalRpowUpperCheck (p q : ℕ) (a u : ℚ) : Bool :=
  decide (0 < q ∧ 0 ≤ a ∧ 0 ≤ u ∧ a ^ p ≤ u ^ q)

theorem rationalRpowUpperCheck_sound (p q : ℕ) (a u : ℚ)
    (h : rationalRpowUpperCheck p q a u = true) :
    (a : ℝ) ^ ((p : ℝ) / (q : ℝ)) ≤ (u : ℝ) := by
  have hc : 0 < q ∧ 0 ≤ a ∧ 0 ≤ u ∧ a ^ p ≤ u ^ q := of_decide_eq_true h
  apply rpow_nat_div_le_of_pow_le (by exact_mod_cast hc.2.1)
    (by exact_mod_cast hc.2.2.1) p q hc.1
  exact_mod_cast hc.2.2.2

/-- An independently bounded actual entry can be substituted below the
rational base. This is the interface used after chronological multiplication. -/
theorem rationalRpowUpperCheck_sound_of_le (p q : ℕ) (a u : ℚ)
    (h : rationalRpowUpperCheck p q a u = true)
    {x : ℝ} (hx : 0 ≤ x) (hxa : x ≤ (a : ℝ)) :
    x ^ ((p : ℝ) / (q : ℝ)) ≤ (u : ℝ) := by
  exact (Real.rpow_le_rpow hx hxa (div_nonneg (Nat.cast_nonneg p)
    (Nat.cast_nonneg q))).trans (rationalRpowUpperCheck_sound p q a u h)

/-- The chronological rational product is computed by the checker, rather
than accepted as an externally asserted intermediate matrix. -/
def rationalPairRpowUpperCheck (n p q : ℕ)
    (A B U : Matrix (Fin n) (Fin n) ℚ) : Bool :=
  decide (∀ i j, rationalRpowUpperCheck p q ((A * B) i j) (U i j) = true)

/-- Every real entry of the actual ordered rational pair has its checked
fractional upper endpoint. Physical-row soundness remains a separate premise. -/
theorem rationalPairRpowUpperCheck_sound (n p q : ℕ)
    (A B U : Matrix (Fin n) (Fin n) ℚ)
    (h : rationalPairRpowUpperCheck n p q A B U = true) (i j : Fin n) :
    (∑ k : Fin n, (A i k : ℝ) * (B k j : ℝ)) ^
      ((p : ℝ) / (q : ℝ)) ≤ (U i j : ℝ) := by
  have hc : ∀ i j, rationalRpowUpperCheck p q ((A * B) i j) (U i j) = true :=
    of_decide_eq_true h
  have hv := rationalRpowUpperCheck_sound p q ((A * B) i j) (U i j) (hc i j)
  simpa only [Matrix.mul_apply, Rat.cast_sum, Rat.cast_mul] using hv

#print axioms rationalPairRpowUpperCheck_sound
#print axioms rationalRpowUpperCheck_sound
#print axioms rationalRpowUpperCheck_sound_of_le
end Spin.Turbo
