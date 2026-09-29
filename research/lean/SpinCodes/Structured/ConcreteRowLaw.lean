import SpinCodes.Structured.ConcreteShuffle

namespace Spin.Structured.Routing

open Finset

noncomputable def rowDensity {n : ℕ} (S : Finset (Fin n)) : ℝ := (S.card : ℝ) / n

lemma rowDensity_nonneg {n : ℕ} (S : Finset (Fin n)) : 0 ≤ rowDensity S := by
  unfold rowDensity
  positivity

lemma rowDensity_le_one {n : ℕ} (S : Finset (Fin n)) : rowDensity S ≤ 1 := by
  unfold rowDensity
  exact div_le_one_of_le₀ (by exact_mod_cast card_le_width S) (by positivity)

noncomputable def rowBernoulli {n : ℕ} (S : Finset (Fin n)) : FinPMF (Finset (Fin n)) :=
  poissonBinom (fun _ => rowDensity S) (fun _ => rowDensity_nonneg S)
    (fun _ => rowDensity_le_one S)

noncomputable def rowCost {n : ℕ} (S : Finset (Fin n)) : ℝ :=
  if S = ∅ then 1 else (n : ℝ) + 1

lemma rowCost_nonneg {n : ℕ} (S : Finset (Fin n)) : 0 ≤ rowCost S := by
  unfold rowCost
  split_ifs <;> positivity

lemma shuffleLaw_dominates_rowBernoulli {n : ℕ} (hn : 0 < n)
    (S : Finset (Fin n)) : Dominates (shuffleLaw S) (rowBernoulli S) (rowCost S) := by
  classical
  intro T
  by_cases he : S = ∅
  · subst S
    by_cases ht : T = ∅
    · subst T
      simp [shuffleLaw_apply, rowBernoulli, poissonBinom_const_apply,
        rowDensity, rowCost]
    · have hc : T.card ≠ 0 := Finset.card_ne_zero.mpr (Finset.nonempty_iff_ne_empty.mpr ht)
      rw [shuffleLaw_apply, Finset.card_empty]
      rw [if_neg hc]
      exact mul_nonneg (rowCost_nonneg _) ((rowBernoulli _).nonneg T)
  · by_cases ht : T.card = S.card
    · rw [shuffleLaw_apply, if_pos ht]
      have hm := types_lower_bound_at hn (card_le_width S)
        (rowDensity_nonneg S) (rowDensity_le_one S)
      have hn1 : (0 : ℝ) < (n : ℝ) + 1 := by positivity
      have hc := choose_card_pos S
      have hmul := (div_le_iff₀ hn1).mp hm
      rw [div_le_iff₀ hc]
      unfold rowBernoulli
      rw [poissonBinom_const_apply, ht]
      simpa only [rowCost, if_neg he,
        rowDensity, mul_assoc, mul_comm, mul_left_comm] using hmul
    · rw [shuffleLaw_apply, if_neg ht]
      exact mul_nonneg (rowCost_nonneg _) (((rowBernoulli _).nonneg T))

end Spin.Structured.Routing



