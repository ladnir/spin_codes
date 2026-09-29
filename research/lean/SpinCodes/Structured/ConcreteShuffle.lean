import SpinCodes.Structured.Domination
import SpinCodes.Structured.Interleaver

namespace Spin.Structured.Routing

open Finset

lemma card_le_width {n : ℕ} (S : Finset (Fin n)) : S.card ≤ n := by
  simpa using Finset.card_le_card (Finset.subset_univ S)

lemma choose_card_pos {n : ℕ} (S : Finset (Fin n)) :
    (0 : ℝ) < (n.choose S.card : ℝ) := by
  exact_mod_cast Nat.choose_pos (card_le_width S)

lemma card_layer_filter (n k : ℕ) :
    (univ.filter (fun T : Finset (Fin n) => T.card = k)).card = n.choose k := by
  have h : univ.filter (fun T : Finset (Fin n) => T.card = k) =
      Finset.powersetCard k univ := by
    ext T
    simp [Finset.mem_powersetCard]
  rw [h, Finset.card_powersetCard, Finset.card_univ, Fintype.card_fin]

/-- The uniform probability law on the weight layer of the supplied support. -/
noncomputable def shuffleLaw {n : ℕ} (S : Finset (Fin n)) :
    FinPMF (Finset (Fin n)) where
  p := fun T => if T.card = S.card then 1 / (n.choose S.card : ℝ) else 0
  nonneg := fun T => by split_ifs <;> positivity
  total := by
    rw [← Finset.sum_filter, Finset.sum_const, nsmul_eq_mul, card_layer_filter]
    exact mul_one_div_cancel (choose_card_pos S).ne'

@[simp] lemma shuffleLaw_apply {n : ℕ} (S T : Finset (Fin n)) :
    (shuffleLaw S).p T = if T.card = S.card then 1 / (n.choose S.card : ℝ) else 0 := rfl

lemma shuffleLaw_fiber_uniform {n : ℕ} (S T U : Finset (Fin n))
    (h : T.card = U.card) : (shuffleLaw S).p T = (shuffleLaw S).p U := by
  simp only [shuffleLaw_apply, h]

lemma shuffleLaw_prob_card {n : ℕ} (S : Finset (Fin n)) (k : ℕ) :
    (shuffleLaw S).prob (fun T => T.card = k) = if S.card = k then 1 else 0 := by
  classical
  by_cases h : S.card = k
  · subst k
    rw [FinPMF.prob]
    rw [Finset.sum_congr rfl (fun T hT =>
      show (shuffleLaw S).p T = 1 / (n.choose S.card : ℝ) from by
        rw [shuffleLaw_apply, if_pos (Finset.mem_filter.mp hT).2])]
    rw [Finset.sum_const, nsmul_eq_mul, card_layer_filter]
    simp only [if_pos rfl]
    exact mul_one_div_cancel (choose_card_pos S).ne'
  · rw [FinPMF.prob, if_neg h]
    apply Finset.sum_eq_zero
    intro T hT
    have ht := (Finset.mem_filter.mp hT).2
    simp [shuffleLaw_apply, ht, Ne.symm h]

end Spin.Structured.Routing

