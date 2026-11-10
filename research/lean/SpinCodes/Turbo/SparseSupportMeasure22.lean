import SpinCodes.Turbo.SparseSupportArray22
import SpinCodes.Turbo.RemainingRoute22

/-!
# Exact RS8 support-class expected measure

Audit: `RS8ShuffledPlacement.rs8_expected_shuffled_word_count` gives each
complete-word atom Cv/(choose(64,v)*15^v), and `RowShuffledArray22` factors
these expected atoms between groups. We disintegrate that existing measure
by exact supports; no probability law on independently resampled messages
is introduced. This is the exact-value counterpart of the old
`ConcreteOuterCountingProduct.tupleCount_eq_prod` bridge.
-/
noncomputable section
namespace Spin.Turbo
open Finset

theorem regionalPacketSupport22_eq_empty_iff {n : ℕ} (word : Fin n → Packet4) :
    regionalPacketSupport22 word = ∅ ↔ word = 0 := by
  constructor
  · intro hs
    funext i
    by_contra hn
    have hi : i ∈ regionalPacketSupport22 word :=
      Finset.mem_filter.mpr ⟨Finset.mem_univ _, hn⟩
    simpa [hs] using hi
  · rintro rfl
    simp [regionalPacketSupport22]

/-- Expected support-class multiplicity, with exact nonzero labels
normalized separately by regionalNonzeroGiven22. -/
def sparseSupportCount22 (points : Fin 8 → Field16)
    (active : Prop) [Decidable active] (T : Finset (Fin 64)) : ℝ :=
  if active then rs8ExpectedPacketCount points T.card / (Nat.choose 64 T.card : ℝ)
  else if T = ∅ then 1 else 0

theorem turboShuffledPositionCount22_eq_support_measure
    (points : Fin 8 → Field16) (active : Prop) [Decidable active]
    (word : PacketWord64) :
    turboShuffledPositionCount22 points active word =
      ∑ T : Finset (Fin 64), sparseSupportCount22 points active T *
        (regionalNonzeroGiven22 T).p word := by
  simp_rw [regionalNonzeroGiven22_apply]
  simp only [mul_ite, mul_zero, Finset.sum_ite_eq, Finset.mem_univ, ite_true]
  by_cases ha : active
  · simp only [sparseSupportCount22, turboShuffledPositionCount22, ha, ite_true]
    rw [rs8_expected_shuffled_word_count]
    change rs8ExpectedPacketCount points (regionalPacketSupport22 word).card /
      ((Nat.choose 64 (regionalPacketSupport22 word).card : ℝ) * 15 ^ (regionalPacketSupport22 word).card) = _
    simp only [div_eq_mul_inv, mul_inv_rev, mul_one]
    ring
  · simp only [sparseSupportCount22, turboShuffledPositionCount22, ha, ite_false]
    rw [shuffledRS8RowLaw_zero_apply]
    by_cases hs : regionalPacketSupport22 word = ∅
    · have hw : word = 0 := (regionalPacketSupport22_eq_empty_iff word).mp hs
      rw [hs]
      simp [hw]
    · have hw : word ≠ 0 := mt (regionalPacketSupport22_eq_empty_iff word).mpr hs
      simp only [hs, hw, ite_false, zero_mul]

/-- Exact disintegration of the original active-set expected packet array
measure into complete-valued exact-support product laws. -/
theorem turboShuffledTupleCount22_eq_support_measure
    (R : ℕ) (points : Fin 8 → Field16) (S : Finset (Fin (R * 16)))
    (words : TurboOuterWord22 R) :
    turboShuffledTupleCount22 R points S words =
      ∑ supports : Fin (R * 16) → Finset (Fin 64),
        (∏ i, sparseSupportCount22 points (i ∈ S) (supports i)) *
          (sparseSupportArrayLaw22 R supports).p words := by
  rw [turboShuffledTupleCount22_eq_prod]
  simp_rw [turboShuffledPositionCount22_eq_support_measure]
  rw [Spin.prod_sum_eq_sum_prod]
  apply Finset.sum_congr rfl
  intro supports _
  rw [Finset.prod_mul_distrib]
  rfl

#print axioms turboShuffledPositionCount22_eq_support_measure
#print axioms turboShuffledTupleCount22_eq_support_measure
end Spin.Turbo

