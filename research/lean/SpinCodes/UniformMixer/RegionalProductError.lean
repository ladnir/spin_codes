import SpinCodes.UniformMixer.RegionalRounding

/-!
# Product comparison from rounded and exact bit gaps

For a sorted support, the rounded product uses gaps in 64-bit encoder
steps. The ideal product uses spacings between the original bit sites.
The row-error calculus bounds their difference by the summed gap error.
-/

noncomputable section

namespace Spin.UniformMixer

open Spin.Structured.FiniteKernel (RowError)

def roundedSupportProduct (γ : ℝ) (R : ℕ) :
    List (Fin (64 * R)) → ℕ → Matrix Bool Bool ℝ
  | [], prev => timeRegionGap γ (((R - prev : ℕ) : ℝ) / R)
  | x :: xs, prev =>
      timeRegionGap γ (((x.val / 64 - prev : ℕ) : ℝ) / R) *
        coarseP * roundedSupportProduct γ R xs (x.val / 64 + 1)

def idealSupportProduct (γ : ℝ) (R : ℕ) :
    List (Fin (64 * R)) → ℕ → Matrix Bool Bool ℝ
  | [], prevBit =>
      timeRegionGap γ (1 - (prevBit : ℝ) / (64 * R))
  | x :: xs, prevBit =>
      timeRegionGap γ (((x.val : ℝ) - prevBit) / (64 * R)) *
        coarseP * idealSupportProduct γ R xs x.val

theorem supportGapError_nonneg (R : ℕ)
    (xs : List (Fin (64 * R))) (prev prevBit : ℕ) :
    0 ≤ supportGapError R xs prev prevBit := by
  induction xs generalizing prev prevBit with
  | nil => exact abs_nonneg _
  | cons x tail ih =>
      exact add_nonneg (abs_nonneg _) (ih _ _)

theorem roundedSupportProduct_substochastic {R prev : ℕ}
    {xs : List (Fin (64 * R))} (h : ValidSupport R prev xs)
    {γ : ℝ} (hγ : 0 ≤ γ) :
    Spin.Structured.FiniteKernel.Substochastic
      (roundedSupportProduct γ R xs prev) := by
  induction h with
  | nil _ =>
      exact timeRegionGap_substochastic hγ (by positivity)
  | cons _ _ _ ih =>
      exact ((timeRegionGap_substochastic hγ (by positivity)).mul
        coarseP_substochastic).mul ih

theorem supportProducts_rowError {R prev : ℕ}
    {xs : List (Fin (64 * R))} (h : ValidSupport R prev xs)
    (hR : 0 < R) {γ : ℝ} (hγ : 0 ≤ γ) :
    ∀ prevBit : ℕ,
      ((prev = 0 ∧ prevBit = 0) ∨ prev = prevBit / 64 + 1) →
      prevBit < 64 * R →
      RowError (roundedSupportProduct γ R xs prev)
        (idealSupportProduct γ R xs prevBit)
        (γ * supportGapError R xs prev prevBit) := by
  induction h with
  | nil hprev =>
      intro prevBit hrelation hbit
      have hRp : (0 : ℝ) < R := by exact_mod_cast hR
      have hden : (0 : ℝ) < 64 * R := by positivity
      have hbitReal : (prevBit : ℝ) ≤ 64 * R := by
        exact_mod_cast (Nat.le_of_lt hbit)
      have hideal : 0 ≤ 1 - (prevBit : ℝ) / (64 * R) := by
        have hdiv : (prevBit : ℝ) / (64 * R) ≤ 1 :=
          (div_le_iff₀ hden).mpr (by simpa using hbitReal)
        linarith
      simpa [roundedSupportProduct, idealSupportProduct,
        supportGapError] using
        (timeRegionGap_rowError hγ (by positivity :
          0 ≤ (((R - _ : ℕ) : ℝ) / R)) hideal)
  | @cons prev x xs hgap hinterior htail ih =>
      intro prevBit hrelation hbit
      have hRp : (0 : ℝ) < R := by exact_mod_cast hR
      have hden : (0 : ℝ) < 64 * R := by positivity
      have hbit_le : prevBit ≤ x.val := by
        rcases hrelation with ⟨hp, hb⟩ | hp
        · omega
        · have hprevmod := Nat.mod_lt prevBit (by omega : 0 < 64)
          have hprevdiv := Nat.div_add_mod prevBit 64
          have hxdiv := Nat.div_add_mod x.val 64
          omega
      have hideal : 0 ≤ ((x.val : ℝ) - prevBit) / (64 * R) := by
        apply div_nonneg _ hden.le
        exact sub_nonneg.mpr (by exact_mod_cast hbit_le)
      have hhead := timeRegionGap_rowError hγ
        (by positivity :
          0 ≤ (((x.val / 64 - prev : ℕ) : ℝ) / R)) hideal
      have hfirst := hhead.mul_right coarseP_substochastic
      have htailErr := ih x.val (Or.inr rfl) x.isLt
      have htailNonneg : 0 ≤ γ *
          supportGapError R xs (x.val / 64 + 1) x.val :=
        mul_nonneg hγ (supportGapError_nonneg R xs _ _)
      have hprod := hfirst.mul htailErr
        (roundedSupportProduct_substochastic htail hγ)
        ((timeRegionGap_substochastic hγ hideal).mul
          coarseP_substochastic)
        htailNonneg
      have heq :
          γ * |(((x.val / 64 - prev : ℕ) : ℝ) / R) -
            ((x.val : ℝ) - prevBit) / (64 * R)| +
          γ * supportGapError R xs (x.val / 64 + 1) x.val =
          γ * supportGapError R (x :: xs) prev prevBit := by
        simp only [supportGapError]
        ring
      rw [heq] at hprod
      exact hprod

end Spin.UniformMixer
