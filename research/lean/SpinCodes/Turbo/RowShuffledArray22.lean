import SpinCodes.Turbo.GroupArrayLaw22
import SpinCodes.Turbo.RS8ShuffledPlacement

/-!
# The sampled row-shuffle law of a Turbo outer array

The first component of the packet route independently permutes the 64
packet positions of every RS8 group. This file keeps that permutation
in the actual product setup and connects the resulting complete-word
array to the already proved one-group shuffled placement law.

The old counterpart is
`Structured/ConcreteRoutePermutation.uniform_rows`; its support-only
law is replaced by a complete four-bit packet-value law here.
-/

noncomputable section

namespace Spin.Turbo

open Spin.Structured.ConcreteRoute

def piPairEquiv {ι : Type*} [Fintype ι]
    {α β : ι → Type*} :
    ((∀ i, α i) × (∀ i, β i)) ≃ (∀ i, α i × β i) where
  toFun xy i := (xy.1 i, xy.2 i)
  invFun f := (fun i => (f i).1, fun i => (f i).2)
  left_inv := by intro xy; cases xy; rfl
  right_inv := by intro f; funext i; exact Prod.ext rfl rfl

/-- Independent product arrays can be paired coordinatewise. -/
theorem piPMF_prod_pair
    {ι : Type*} [Fintype ι] [DecidableEq ι]
    {α β : ι → Type*}
    [∀ i, Fintype (α i)] [∀ i, Fintype (β i)]
    (P : ∀ i, FinPMF (α i))
    (Q : ∀ i, FinPMF (β i)) :
    ((piPMF P).prod (piPMF Q)).map piPairEquiv =
      piPMF (fun i => (P i).prod (Q i)) := by
  classical
  apply FinPMF.ext
  funext y
  rw [FinPMF.map_equiv_apply]
  change (piPMF P).p (fun i => (y i).1) *
    (piPMF Q).p (fun i => (y i).2) =
      (piPMF (fun i => (P i).prod (Q i))).p y
  simp only [piPMF_apply]
  rw [← Finset.prod_mul_distrib]
  rfl

abbrev TurboRowPerms22 (R : ℕ) :=
  Fin (R * 16) → Equiv.Perm (Fin 64)

def turboRowPermLaw22 (R : ℕ) : FinPMF (TurboRowPerms22 R) :=
  piPMF (fun _ : Fin (R * 16) =>
    FinPMF.uniform (Equiv.Perm (Fin 64)))

def turboRowShuffledWord22 (points : Fin 8 → Field16)
    (msg : TurboMessage22 R)
    (seed : TurboOuterSeed22 R × TurboRowPerms22 R) :
    TurboOuterWord22 R :=
  fun i => shufflePacketWord (seed.2 i)
    (sampledRS8RowWord points (msg i) (seed.1 i))

def turboRowShuffledLaw22 (R : ℕ) (points : Fin 8 → Field16)
    (msg : TurboMessage22 R) : FinPMF (TurboOuterWord22 R) :=
  ((turboOuterLaw22 R).prod (turboRowPermLaw22 R)).map
    (turboRowShuffledWord22 points msg)

/-- For each fixed message, the actual outer multipliers and row
permutations induce independent complete row-word coordinates. -/
theorem turboRowShuffledLaw22_product
    (R : ℕ) (points : Fin 8 → Field16)
    (msg : TurboMessage22 R) :
    turboRowShuffledLaw22 R points msg =
      piPMF (fun i : Fin (R * 16) =>
        shuffledPacketLaw (sampledRS8RowLaw points (msg i))) := by
  let P : Fin (R * 16) → FinPMF
      (Fin 8 → {c : Field32 // c ≠ 0}) :=
    fun _ => piPMF (fun _ : Fin 8 =>
      FinPMF.uniform {c : Field32 // c ≠ 0})
  let Q : Fin (R * 16) → FinPMF (Equiv.Perm (Fin 64)) :=
    fun _ => FinPMF.uniform (Equiv.Perm (Fin 64))
  let f : ∀ i : Fin (R * 16),
      (Fin 8 → {c : Field32 // c ≠ 0}) × Equiv.Perm (Fin 64) →
        PacketWord64 :=
    fun i p => shufflePacketWord p.2
      (sampledRS8RowWord points (msg i) p.1)
  have hlocal (i : Fin (R * 16)) :
      ((P i).prod (Q i)).map (f i) =
        shuffledPacketLaw (sampledRS8RowLaw points (msg i)) := by
    rw [shuffledPacketLaw, sampledRS8RowLaw,
      sampledRS8PacketColumnsLaw]
    apply FinPMF.eq_of_expect_eq
    intro g
    simp only [FinPMF.expect_map, FinPMF.expect_prod]
    rfl
  calc
    turboRowShuffledLaw22 R points msg =
        (((piPMF P).prod (piPMF Q)).map piPairEquiv).map
          (fun pairs i => f i (pairs i)) := by
      rw [FinPMF.map_comp]
      rfl
    _ = (piPMF (fun i => (P i).prod (Q i))).map
          (fun pairs i => f i (pairs i)) := by
      rw [piPMF_prod_pair]
    _ = piPMF (fun i => ((P i).prod (Q i)).map (f i)) :=
      Spin.piPMF_map _ _
    _ = _ := by simp only [hlocal]

/-- Expected number of active-set messages with a specified complete
packet array after the actual independent row shuffles. -/
def turboShuffledTupleCount22 (R : ℕ) (points : Fin 8 → Field16)
    (S : Finset (Fin (R * 16)))
    (words : TurboOuterWord22 R) : ℝ :=
  ∑ msg ∈ turboActiveMessages22 S,
    (turboRowShuffledLaw22 R points msg).p words

def turboShuffledPositionCount22 (points : Fin 8 → Field16)
    (active : Prop) [Decidable active] (word : PacketWord64) : ℝ :=
  if active then rs8ExpectedShuffledWordCount points word
  else (shuffledPacketLaw (sampledRS8RowLaw points 0)).p word

theorem turboShuffledPositionCount22_eq_sum
    (points : Fin 8 → Field16)
    (active : Prop) [Decidable active] (word : PacketWord64) :
    turboShuffledPositionCount22 points active word =
      ∑ msg : RS8Group Field16,
        if (msg ≠ 0 ↔ active) then
          (shuffledPacketLaw (sampledRS8RowLaw points msg)).p word
        else 0 := by
  classical
  by_cases h : active
  · simp [turboShuffledPositionCount22, h,
      rs8ExpectedShuffledWordCount, Spin.nonzeroMsgs,
      Finset.sum_filter]
  · simp [turboShuffledPositionCount22, h]

/-- The shuffled message-count law still factors by group. The regional
permutations and inner state are not included in this product. -/
theorem turboShuffledTupleCount22_eq_prod
    (R : ℕ) (points : Fin 8 → Field16)
    (S : Finset (Fin (R * 16)))
    (words : TurboOuterWord22 R) :
    turboShuffledTupleCount22 R points S words =
      ∏ i : Fin (R * 16),
        turboShuffledPositionCount22 points (i ∈ S) (words i) := by
  classical
  simp only [turboShuffledPositionCount22_eq_sum,
    Spin.prod_sum_eq_sum_prod, turboShuffledTupleCount22,
    turboRowShuffledLaw22_product, Spin.piPMF_apply,
    turboActiveMessages22, Finset.sum_filter]
  apply Finset.sum_congr rfl
  intro msg _
  by_cases h : ∀ i, msg i ≠ 0 ↔ i ∈ S
  · have hs : turboActivePositions22 msg = S :=
      (Finset.mem_filter.mp ((turboActiveMessages22_iff S msg).mpr h)).2
    rw [if_pos hs]
    apply Finset.prod_congr rfl
    intro i _
    simp only [if_pos (h i)]
  · have hs : turboActivePositions22 msg ≠ S := by
      intro he
      exact h ((turboActiveMessages22_iff S msg).mp
        (Finset.mem_filter.mpr ⟨Finset.mem_univ _, he⟩))
    rw [if_neg hs]
    obtain ⟨i, hi⟩ := not_forall.mp h
    symm
    exact Finset.prod_eq_zero (Finset.mem_univ i) (if_neg hi)

/-- At an active group, the exact local factor is the paper's `C_v`
divided by the number of complete words of packet weight `v`. -/
theorem turboShuffledPositionCount22_active
    (points : Fin 8 → Field16) (word : PacketWord64) :
    turboShuffledPositionCount22 points True word =
      rs8ExpectedPacketCount points (packetWordSupport word).card /
        ((Nat.choose 64 (packetWordSupport word).card : ℝ) *
          (15 : ℝ) ^ (packetWordSupport word).card) := by
  simpa [turboShuffledPositionCount22] using
    rs8_expected_shuffled_word_count points word

#print axioms piPMF_prod_pair
#print axioms turboRowShuffledLaw22_product
#print axioms turboShuffledTupleCount22_eq_prod
#print axioms turboShuffledPositionCount22_active

end Spin.Turbo
