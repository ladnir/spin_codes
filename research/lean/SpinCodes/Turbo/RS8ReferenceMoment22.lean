import SpinCodes.Turbo.RS8UniformReference22
import SpinCodes.Turbo.RemainingRoute22

/-!
# Transfer the RS8 reference bound through the actual route and inner

The completed Structured proof first bounds every complete outer word,
then multiplies by the nonnegative conditional failure probability. The
same step is valid here because the sampled regional permutations and the
continuous inner state are kept together in `turboRemainingFailure22`.
-/

noncomputable section

namespace Spin.Turbo

open Finset

theorem finite_weighted_sum_le_factor
    {α : Type*} [Fintype α]
    (a b f : α → ℝ) (c : ℝ)
    (hab : ∀ x, a x ≤ c * b x) (hf : ∀ x, 0 ≤ f x) :
    (∑ x, a x * f x) ≤ c * ∑ x, b x * f x := by
  rw [Finset.mul_sum]
  apply Finset.sum_le_sum
  intro x _
  calc
    a x * f x ≤ (c * b x) * f x := mul_le_mul_of_nonneg_right (hab x) (hf x)
    _ = c * (b x * f x) := by ring

theorem turboShuffledCountMoment22_le_reference
    (R : ℕ) (points : Fin 8 → Field16)
    (hpoints : Function.Injective points)
    (S : Finset (Fin (R * 16))) (d : ℕ) :
    (∑ words : TurboOuterWord22 R,
      turboShuffledTupleCount22 R points S words *
        turboRemainingFailure22 words d) ≤
    turboOldUniformFactor22 ^ S.card *
      ∑ words : TurboOuterWord22 R,
        (∏ i : Fin (R * 16),
          turboReferenceRowDensity22 points (i ∈ S) (words i)) *
          turboRemainingFailure22 words d := by
  simp only [turboShuffledTupleCount22_eq_prod]
  apply finite_weighted_sum_le_factor
  · intro words
    exact turboPositionProduct22_le_reference R points hpoints S words
  · intro words
    exact (turboRemainingLaw22 R).prob_nonneg _

/-- The comparison is on the original shared-outer setup experiment,
not a fresh setup independently sampled for each message. -/
theorem turboSetup22_activeSet_moment_le_reference
    (R : ℕ) (points : Fin 8 → Field16)
    (hpoints : Function.Injective points)
    (S : Finset (Fin (R * 16))) (d : ℕ) :
    (turboOuterLaw22 R).expect (fun outerSeed =>
      ∑ msg ∈ turboActiveMessages22 S,
        (turboSetup22 R points).qd d
          (turboOuterWord22 points outerSeed msg)) ≤
    turboOldUniformFactor22 ^ S.card *
      ∑ words : TurboOuterWord22 R,
        (∏ i : Fin (R * 16),
          turboReferenceRowDensity22 points (i ∈ S) (words i)) *
          turboRemainingFailure22 words d := by
  rw [turboSetup22_activeSet_expect_eq_shuffled_count]
  exact turboShuffledCountMoment22_le_reference R points hpoints S d

#print axioms turboShuffledCountMoment22_le_reference
#print axioms turboSetup22_activeSet_moment_le_reference

end Spin.Turbo
