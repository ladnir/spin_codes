import SpinCodes.Structured.ConcretePlacementLaw

/-! Collision and short-gap union bounds under the actual fixed-cardinality region shuffle. -/
noncomputable section
namespace Spin.Structured.Placement
open Finset Routing
attribute [local instance] Classical.propDecidable

theorem shuffle_prob_hits {n : Nat} (S B : Finset (Fin n)) :
    (shuffleLaw S).prob (fun T => ∃ x ∈ B, x ∈ T) ≤ (B.card : ℝ) * S.card / n := by
  calc _ ≤ ∑ x ∈ B, (shuffleLaw S).prob (fun T => x ∈ T) :=
      Spin.FinPMF.prob_biUnion_le _ B _
    _ = _ := by
      simp only [shuffle_prob_mem, Finset.sum_const, nsmul_eq_mul]
      ring

theorem shuffle_prob_pairs {n : Nat} (S : Finset (Fin n)) (B : Finset (Fin n × Fin n))
    (hB : ∀ p ∈ B, p.1 ≠ p.2) :
    (shuffleLaw S).prob (fun T => ∃ p ∈ B, p.1 ∈ T ∧ p.2 ∈ T) ≤
      (B.card : ℝ) * (((S.card : ℝ) * (S.card - 1)) / ((n : ℝ) * (n - 1))) := by
  calc _ ≤ ∑ p ∈ B, (shuffleLaw S).prob (fun T => p.1 ∈ T ∧ p.2 ∈ T) :=
      Spin.FinPMF.prob_biUnion_le _ B _
    _ = _ := by
      rw [Finset.sum_congr rfl (fun p hp => shuffle_prob_pair S p.1 p.2 (hB p hp))]
      simp only [Finset.sum_const, nsmul_eq_mul]

def closePairs (n D : Nat) : Finset (Fin n × Fin n) :=
  univ.filter (fun p => p.1 ≠ p.2 ∧ Nat.dist p.1.val p.2.val ≤ D)

theorem close_degree_bound {n : Nat} (D : Nat) (x : Fin n) :
    (univ.filter (fun y : Fin n => x ≠ y ∧ Nat.dist x.val y.val ≤ D)).card ≤ 2 * D + 1 := by
  let B := univ.filter (fun y : Fin n => x ≠ y ∧ Nat.dist x.val y.val ≤ D)
  have hs : B.image Fin.val ⊆ Finset.Icc (x.val - D) (x.val + D) := by
    intro y hy
    obtain ⟨i, hi, rfl⟩ := Finset.mem_image.mp hy
    have hd := (Finset.mem_filter.mp hi).2.2
    simp only [Finset.mem_Icc]
    unfold Nat.dist at hd
    omega
  calc B.card = (B.image Fin.val).card := (Finset.card_image_of_injective _ Fin.val_injective).symm
    _ ≤ (Finset.Icc (x.val - D) (x.val + D)).card := Finset.card_le_card hs
    _ ≤ 2 * D + 1 := by rw [Nat.card_Icc]; omega

theorem closePairs_card_bound (n D : Nat) : (closePairs n D).card ≤ n * (2 * D + 1) := by
  unfold closePairs
  rw [Finset.card_filter, Fintype.sum_prod_type]
  calc (∑ x : Fin n, ∑ y : Fin n, if x ≠ y ∧ Nat.dist x.val y.val ≤ D then 1 else 0) ≤
      ∑ _x : Fin n, (2 * D + 1) := by
        apply Finset.sum_le_sum
        intro x _
        simpa only [Finset.card_filter] using close_degree_bound D x
    _ = _ := by simp only [Finset.sum_const, Finset.card_univ, Fintype.card_fin, nsmul_eq_mul, Nat.cast_id]

theorem card_times_pred_nonneg {n : Nat} (S : Finset (Fin n)) :
    0 ≤ (S.card : ℝ) * (S.card - 1) := by
  by_cases h : S.card = 0
  · simp only [h, Nat.cast_zero, zero_mul, le_refl]
  · have hn : 1 ≤ S.card := by omega
    have hr : (1 : ℝ) ≤ S.card := by exact_mod_cast hn
    exact mul_nonneg (Nat.cast_nonneg _) (by linarith)

theorem shuffle_prob_close {n : Nat} (hn : 2 ≤ n) (S : Finset (Fin n)) (D : Nat) :
    (shuffleLaw S).prob (fun T => ∃ x ∈ T, ∃ y ∈ T, x ≠ y ∧ Nat.dist x.val y.val ≤ D) ≤
      (2 * (D : ℝ) + 1) * (S.card : ℝ) ^ 2 / ((n : ℝ) - 1) := by
  have hnR : (2 : ℝ) ≤ n := by exact_mod_cast hn
  have hn0 : (n : ℝ) ≠ 0 := by linarith
  have hn1 : (n : ℝ) - 1 > 0 := by linarith
  have hB (p : Fin n × Fin n) (hp : p ∈ closePairs n D) : p.1 ≠ p.2 :=
    (Finset.mem_filter.mp hp).2.1
  have h := shuffle_prob_pairs S (closePairs n D) hB
  have he : (shuffleLaw S).prob
      (fun T => ∃ x ∈ T, ∃ y ∈ T, x ≠ y ∧ Nat.dist x.val y.val ≤ D) =
      (shuffleLaw S).prob (fun T => ∃ p ∈ closePairs n D, p.1 ∈ T ∧ p.2 ∈ T) := by
    apply Spin.FinPMF.prob_congr
    intro T
    simp only [closePairs, Finset.mem_filter, Finset.mem_univ, true_and]
    constructor
    · rintro ⟨x, hx, y, hy, hxy, hd⟩
      exact ⟨(x, y), ⟨hxy, hd⟩, hx, hy⟩
    · rintro ⟨⟨x, y⟩, ⟨hxy, hd⟩, hx, hy⟩
      exact ⟨x, hx, y, hy, hxy, hd⟩
  rw [he]
  have hcard : ((closePairs n D).card : ℝ) ≤ (n : ℝ) * (2 * (D : ℝ) + 1) := by
    exact_mod_cast closePairs_card_bound n D
  have hp : 0 ≤ ((S.card : ℝ) * (S.card - 1)) / ((n : ℝ) * (n - 1)) :=
    div_nonneg (card_times_pred_nonneg S) (mul_nonneg (Nat.cast_nonneg _) hn1.le)
  calc _ ≤ ((closePairs n D).card : ℝ) *
      (((S.card : ℝ) * (S.card - 1)) / ((n : ℝ) * (n - 1))) := h
    _ ≤ ((n : ℝ) * (2 * (D : ℝ) + 1)) *
      (((S.card : ℝ) * (S.card - 1)) / ((n : ℝ) * (n - 1))) := mul_le_mul_of_nonneg_right hcard hp
    _ = (2 * (D : ℝ) + 1) * ((S.card : ℝ) * (S.card - 1)) / ((n : ℝ) - 1) := by
      field_simp
    _ ≤ (2 * (D : ℝ) + 1) * (S.card : ℝ) ^ 2 / ((n : ℝ) - 1) := by
      apply div_le_div_of_nonneg_right _ hn1.le
      apply mul_le_mul_of_nonneg_left _ (by positivity)
      nlinarith [Nat.cast_nonneg (α := ℝ) S.card]

end Spin.Structured.Placement
