import SpinCodes.Structured.ConcreteOuterEntropy
import SpinCodes.Structured.Composition

noncomputable section
namespace Spin.Structured.ConcreteOuter
open Spin.Numeric

/-- Exact finite accumulator exponent retaining both parity-dependent binomial arguments. -/
def transitionEntropy (b a c : ℕ) : ℝ :=
  if a = 0 then 0 else
    layerEntropy (c - 1) ((a + 1) / 2 - 1) +
      layerEntropy (b - c) (a - (a + 1) / 2) - layerEntropy b a

theorem choose_product_div_le_entropy {n₁ k₁ n₂ k₂ b a : ℕ}
    (h₁ : k₁ ≤ n₁) (h₂ : k₂ ≤ n₂) (ha : a ≤ b) :
    ((n₁.choose k₁ : ℝ) * (n₂.choose k₂ : ℝ)) / (b.choose a : ℝ) ≤
      ((b : ℝ) + 1) * Real.exp (layerEntropy n₁ k₁ + layerEntropy n₂ k₂ - layerEntropy b a) := by
  have hb : (0 : ℝ) < (b : ℝ) + 1 := by positivity
  have hc : (0 : ℝ) < (b.choose a : ℝ) := by exact_mod_cast Nat.choose_pos ha
  have he : (0 : ℝ) < Real.exp (layerEntropy b a) / ((b : ℝ) + 1) := by positivity
  calc ((n₁.choose k₁ : ℝ) * (n₂.choose k₂ : ℝ)) / (b.choose a : ℝ) ≤
      (Real.exp (layerEntropy n₁ k₁) * Real.exp (layerEntropy n₂ k₂)) / (b.choose a : ℝ) := by
        exact div_le_div_of_nonneg_right
          (mul_le_mul (choose_le_exp_entropy h₁) (choose_le_exp_entropy h₂)
            (by positivity) (Real.exp_pos _).le) hc.le
    _ ≤ (Real.exp (layerEntropy n₁ k₁) * Real.exp (layerEntropy n₂ k₂)) /
        (Real.exp (layerEntropy b a) / ((b : ℝ) + 1)) :=
      div_le_div_of_nonneg_left (by positivity) he (exp_entropy_div_le_choose ha)
    _ = _ := by
      rw [Real.exp_sub, Real.exp_add]
      field_simp

/-- Every accumulator transition has an explicit finite entropy envelope.
Zero probabilities, boundary layers, and odd input weights are all included. -/
theorem transition_le_exp_entropy {b a : ℕ} (ha : a ≤ b) (c : ℕ) :
    Pt b a c ≤ ((b : ℝ) + 1) * Real.exp (transitionEntropy b a c) := by
  by_cases ha0 : a = 0
  · subst a
    simp only [Pt, accT_zero_left, transitionEntropy, if_pos rfl, Real.exp_zero,
      mul_one, Nat.choose_zero_right, Nat.cast_one, div_one]
    split_ifs <;> simp <;> positivity
  by_cases ht : accT b a c = 0
  · simp only [Pt, ht, Nat.cast_zero, zero_div]
    positivity
  have hc0 : c ≠ 0 := by
    intro h
    subst c
    exact ht (accT_zero_right b ha0)
  have hcb : c ≤ b := by
    by_contra h
    exact ht (accT_of_lt b a c ha0 (by omega))
  have he : accT b a c = (c - 1).choose ((a + 1) / 2 - 1) *
      (b - c).choose (a - (a + 1) / 2) := by
    simp only [accT, if_neg ha0, if_neg (show ¬(c = 0 ∨ b < c) by omega)]
  have h₁ : (a + 1) / 2 - 1 ≤ c - 1 := by
    by_contra h
    apply ht
    rw [he, Nat.choose_eq_zero_of_lt (by omega), zero_mul]
  have h₂ : a - (a + 1) / 2 ≤ b - c := by
    by_contra h
    apply ht
    rw [he, Nat.choose_eq_zero_of_lt (n := b - c) (k := a - (a + 1) / 2) (by omega), mul_zero]
  rw [Pt, he, Nat.cast_mul, transitionEntropy, if_neg ha0]
  exact choose_product_div_le_entropy h₁ h₂ ha

end Spin.Structured.ConcreteOuter

