import SpinCodes.Turbo.RS8AtomEnvelope
import SpinCodes.Turbo.RowShuffledArray22

/-!
# RS8 outer comparison with a uniform-value group reference

The old Structured proof uses `ConcreteOuterCountingProduct.tupleCount_le_fairRows`
to multiply a local pointwise count bound, followed by
`Dominates.expect_le` to transfer a nonnegative route-and-inner test.
Here each active Turbo group uses the uniform complete 256-bit packet-word
law as its reference. Inactive groups retain their exact zero-message law.
-/

noncomputable section

namespace Spin.Turbo

open Finset

def turboOldUniformFactor22 : ℝ :=
  (2 : ℝ) ^ 256 / (65535 : ℝ) ^ 8

def turboReferenceRowDensity22 (points : Fin 8 → Field16)
    (active : Prop) [Decidable active] (word : PacketWord64) : ℝ :=
  if active then 1 / (2 : ℝ) ^ 256
  else (shuffledPacketLaw (sampledRS8RowLaw points 0)).p word

theorem turboOldUniformFactor22_pos : 0 < turboOldUniformFactor22 := by
  unfold turboOldUniformFactor22
  positivity

theorem turboOldUniformFactor22_uniform_atom :
    turboOldUniformFactor22 * (1 / (2 : ℝ) ^ 256) =
      1 / (65535 : ℝ) ^ 8 := by
  unfold turboOldUniformFactor22
  have hc : (2 : ℝ) ^ 256 ≠ 0 := by positivity
  rw [div_mul_eq_mul_div, mul_one_div_cancel hc]

theorem turboShuffledPositionCount22_active_le_reference
    (points : Fin 8 → Field16) (hpoints : Function.Injective points)
    (word : PacketWord64) :
    turboShuffledPositionCount22 points True word ≤
      turboOldUniformFactor22 * (1 / (2 : ℝ) ^ 256) := by
  change rs8ExpectedShuffledWordCount points word ≤ _
  rw [turboOldUniformFactor22_uniform_atom]
  exact rs8ExpectedShuffledWordCount_le_old_uniform_atom points hpoints word

theorem turboShuffledPositionCount22_nonneg
    (points : Fin 8 → Field16) (active : Prop)
    [Decidable active] (word : PacketWord64) :
    0 ≤ turboShuffledPositionCount22 points active word := by
  by_cases h : active
  · rw [turboShuffledPositionCount22, if_pos h]
    unfold rs8ExpectedShuffledWordCount
    apply Finset.sum_nonneg
    intro rows _
    exact (shuffledPacketLaw (sampledRS8RowLaw points rows)).nonneg word
  · rw [turboShuffledPositionCount22, if_neg h]
    exact (shuffledPacketLaw (sampledRS8RowLaw points 0)).nonneg word

theorem turboShuffledPositionCount22_le_reference
    (points : Fin 8 → Field16) (hpoints : Function.Injective points)
    (active : Prop) [Decidable active] (word : PacketWord64) :
    turboShuffledPositionCount22 points active word ≤
      (if active then turboOldUniformFactor22 else 1) *
        turboReferenceRowDensity22 points active word := by
  by_cases h : active
  · simpa only [h, ite_true, turboReferenceRowDensity22] using
      turboShuffledPositionCount22_active_le_reference points hpoints word
  · simp only [h, ite_false, one_mul, le_refl,
      turboShuffledPositionCount22, turboReferenceRowDensity22]

/-- Generic finite product comparison used by the old Structured outer
counting proof. Keeping the arithmetic separate avoids kernel expansion
of the Turbo word space in the final specialization. -/
theorem product_le_active_factor
    {ι : Type*} [Fintype ι] [DecidableEq ι]
    (S : Finset ι) (c : ℝ) (a b : ι → ℝ)
    (ha : ∀ i, 0 ≤ a i)
    (hab : ∀ i, a i ≤ (if i ∈ S then c else 1) * b i) :
    (∏ i, a i) ≤ c ^ S.card * ∏ i, b i := by
  calc
    (∏ i, a i) ≤ ∏ i, (if i ∈ S then c else 1) * b i :=
      Finset.prod_le_prod₀ (fun i _ => ha i) (fun i _ => hab i)
    _ = (∏ i, if i ∈ S then c else 1) * ∏ i, b i := by
      rw [Finset.prod_mul_distrib]
    _ = c ^ S.card * ∏ i, b i := by
      rw [Finset.prod_ite]
      simp

theorem turboPositionProduct22_le_reference
    (R : ℕ) (points : Fin 8 → Field16)
    (hpoints : Function.Injective points)
    (S : Finset (Fin (R * 16)))
    (words : TurboOuterWord22 R) :
    (∏ i : Fin (R * 16),
      turboShuffledPositionCount22 points (i ∈ S) (words i)) ≤
      turboOldUniformFactor22 ^ S.card *
        ∏ i : Fin (R * 16),
          turboReferenceRowDensity22 points (i ∈ S) (words i) := by
  apply product_le_active_factor S turboOldUniformFactor22
  · intro i
    exact turboShuffledPositionCount22_nonneg points (i ∈ S) (words i)
  · intro i
    exact turboShuffledPositionCount22_le_reference points hpoints (i ∈ S) (words i)

#print axioms turboPositionProduct22_le_reference

end Spin.Turbo
