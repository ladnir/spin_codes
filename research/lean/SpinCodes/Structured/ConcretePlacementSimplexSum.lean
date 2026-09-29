import SpinCodes.Structured.ConcretePlacementSimplexLaw

/-! Exact good-placement sum and quantitative omitted mass for the actual shuffle law. -/
noncomputable section
namespace Spin.Structured.Placement
open Finset Routing
attribute [local instance] Classical.propDecidable

theorem shuffle_good_blocks_sum {R a : Nat} (S : Finset (Fin (128 * R)))
    (hS : S.card = a) (H : Nat) (f : Finset (Fin (128 * R)) → ℝ) :
    (shuffleLaw S).expect (fun T => if ¬ badBlockPlacement H T then f T else 0) =
      ((128 : ℝ)^a / ((128 * R).choose a : ℝ)) *
        ∑ B : BlockSubset R a, if ¬ badBlockIndices (orderedBlocks B) H then
          (Spin.piPMF (fun _ : Fin a => Spin.FinPMF.uniform (Fin 128))).expect
            (fun coords => f (singletonSupport (orderedBlocks B) coords)) else 0 := by
  have he : (shuffleLaw S).expect (fun T => if ¬ badBlockPlacement H T then f T else 0) =
      (shuffleLaw S).expect (fun T => if (occupiedBlocks T).card = a then
        (if ¬ badBlockPlacement H T then f T else 0) else 0) := by
    unfold Spin.FinPMF.expect
    apply sum_congr rfl
    intro T _
    by_cases ht : T.card = a
    · by_cases hg : ¬ badBlockPlacement H T
      · have hc : (occupiedBlocks T).card = a := (good_occupiedBlocks_card T H hg).trans ht
        simp only [hc, hg, ite_true]
      · simp [hg]
    · simp [shuffleLaw_apply, hS, ht]
  rw [he, shuffle_distinct_blocks_sum S hS]
  congr 1
  apply sum_congr rfl
  intro B _
  simp only [badBlockPlacement_singletonSupport]
  by_cases hg : ¬ badBlockIndices (orderedBlocks B) H
  · simp [hg]
  · simp only [hg, ite_false, Spin.FinPMF.expect_const]

theorem expect_omit_event {Ω : Type*} [Fintype Ω] (P : Spin.FinPMF Ω)
    (E : Ω → Prop) (f : Ω → ℝ) (hf : ∀ x, 0 ≤ f x ∧ f x ≤ 1) :
    |P.expect f - P.expect (fun x => if ¬ E x then f x else 0)| ≤ P.prob E := by
  have he : P.expect f - P.expect (fun x => if ¬ E x then f x else 0) =
      P.expect (fun x => if E x then f x else 0) := by
    unfold Spin.FinPMF.expect
    rw [← sum_sub_distrib]
    apply sum_congr rfl
    intro x _
    by_cases hx : E x <;> simp [hx]
  rw [he, abs_of_nonneg]
  · rw [Spin.FinPMF.prob_eq_expect_indicator]
    apply Spin.FinPMF.expect_mono
    intro x
    by_cases hx : E x <;> simp [hx, (hf x).2]
  · apply Spin.FinPMF.expect_nonneg
    intro x
    split_ifs <;> first | exact (hf x).1 | exact le_rfl

theorem shuffle_good_blocks_sum_error {R a : Nat} (hR : 1 ≤ R)
    (S : Finset (Fin (128 * R))) (hS : S.card = a) (H : Nat)
    (f : Finset (Fin (128 * R)) → ℝ) (hf : ∀ T, 0 ≤ f T ∧ f T ≤ 1) :
    |(shuffleLaw S).expect f -
      ((128 : ℝ)^a / ((128 * R).choose a : ℝ)) *
        ∑ B : BlockSubset R a, if ¬ badBlockIndices (orderedBlocks B) H then
          (Spin.piPMF (fun _ : Fin a => Spin.FinPMF.uniform (Fin 128))).expect
            (fun coords => f (singletonSupport (orderedBlocks B) coords)) else 0| ≤
      256 * (H : ℝ) * a / (128 * R) +
        (256 * ((H : ℝ) + 1) + 1) * (a : ℝ)^2 / (128 * (R : ℝ) - 1) := by
  rw [← shuffle_good_blocks_sum S hS H f]
  have h := (expect_omit_event (shuffleLaw S) (badBlockPlacement H) f hf).trans
    (shuffle_prob_badBlockPlacement hR S H)
  simpa only [hS] using h

theorem card_blockSubset (R a : Nat) : Fintype.card (BlockSubset R a) = R.choose a := by
  rw [Fintype.card_subtype]
  exact card_layer_filter R a

theorem shuffle_prob_distinct_blocks {R a : Nat} (S : Finset (Fin (128 * R))) (hS : S.card = a) :
    (shuffleLaw S).prob (fun T => (occupiedBlocks T).card = a) =
      (128 : ℝ)^a * (R.choose a : ℝ) / ((128 * R).choose a : ℝ) := by
  rw [Spin.FinPMF.prob_eq_expect_indicator, shuffle_distinct_blocks_sum S hS]
  simp only [Spin.FinPMF.expect_const, sum_const, card_univ, card_blockSubset, nsmul_eq_mul, mul_one]
  ring

end Spin.Structured.Placement

