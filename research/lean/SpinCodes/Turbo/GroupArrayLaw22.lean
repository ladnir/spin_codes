import SpinCodes.Turbo.Occupation22
import SpinCodes.FiniteProductLaw

/-!
# Joint outer-word law for one fixed Turbo message

The outer setup has independent RS8 group coordinates. For a fixed full
message, its complete array of four-bit packet values has the product of
the actual sampled single-group laws. This concerns one message at a
time; different messages still use the same sampled setup.

The old proof's counterpart is
`Structured/ConcreteOuterCountingProduct.tupleCount_eq_prod`.
-/

noncomputable section

namespace Spin.Turbo

open Finset

/-- Product law of the complete outer packet rows for a fixed message. -/
theorem turboOuterWord22_product_law
    (R : ℕ) (points : Fin 8 → Field16)
    (msg : TurboMessage22 R) :
    (turboOuterLaw22 R).map
      (fun seed => turboOuterWord22 points seed msg) =
    piPMF (fun i : Fin (R * 16) => sampledRS8RowLaw points (msg i)) := by
  let P : Fin (R * 16) → FinPMF
      (Fin 8 → {c : Field32 // c ≠ 0}) :=
    fun _ => piPMF (fun _ : Fin 8 =>
      FinPMF.uniform {c : Field32 // c ≠ 0})
  have h := Spin.piPMF_map P
    (fun i coeff => sampledRS8RowWord points (msg i) coeff)
  convert h using 1
  · rfl
  · apply congrArg
    funext i
    rw [sampledRS8RowLaw, sampledRS8PacketColumnsLaw,
      FinPMF.map_comp]
    rfl

/-- The expected number of messages with active set `S` whose complete
outer packet array equals `words`. It is a first moment over one shared
outer setup, not a probability law on independent messages. -/
def turboTupleCount22 (R : ℕ) (points : Fin 8 → Field16)
    (S : Finset (Fin (R * 16)))
    (words : TurboOuterWord22 R) : ℝ :=
  ∑ msg ∈ turboActiveMessages22 S,
    ((turboOuterLaw22 R).map
      (fun seed => turboOuterWord22 points seed msg)).p words

/-- Local expected packet-word count at an active group; an inactive
group has the zero-message law. -/
def turboPositionCount22 (points : Fin 8 → Field16)
    (active : Prop) [Decidable active] (word : PacketWord64) : ℝ :=
  if active then
    ∑ msg ∈ Spin.nonzeroMsgs (RS8Group Field16),
      (sampledRS8RowLaw points msg).p word
  else (sampledRS8RowLaw points 0).p word

theorem turboPositionCount22_eq_sum (points : Fin 8 → Field16)
    (active : Prop) [Decidable active] (word : PacketWord64) :
    turboPositionCount22 points active word =
      ∑ msg : RS8Group Field16,
        if (msg ≠ 0 ↔ active) then
          (sampledRS8RowLaw points msg).p word else 0 := by
  classical
  by_cases h : active
  · simp [turboPositionCount22, h, Spin.nonzeroMsgs,
      Finset.sum_filter]
  · simp [turboPositionCount22, h]

theorem turboActiveMessages22_iff {R : ℕ}
    (S : Finset (Fin (R * 16))) (msg : TurboMessage22 R) :
    msg ∈ turboActiveMessages22 S ↔
      ∀ i, msg i ≠ 0 ↔ i ∈ S := by
  simp only [turboActiveMessages22, Finset.mem_filter, Finset.mem_univ,
    true_and]
  constructor
  · intro h i
    have hi := Finset.ext_iff.mp h i
    simpa only [turboActivePositions22, Finset.mem_filter,
      Finset.mem_univ, true_and] using hi
  · intro h
    ext i
    simpa only [turboActivePositions22, Finset.mem_filter,
      Finset.mem_univ, true_and] using h i

/-- Exact factorization of the active-set expected outer array count.
It uses independence of setup coordinates *between groups* only after
fixing one message, then interchanges the finite sum over messages with
the product. It does not assert independent codewords for messages
sharing the same setup. -/
theorem turboTupleCount22_eq_prod
    (R : ℕ) (points : Fin 8 → Field16)
    (S : Finset (Fin (R * 16)))
    (words : TurboOuterWord22 R) :
    turboTupleCount22 R points S words =
      ∏ i : Fin (R * 16),
        turboPositionCount22 points (i ∈ S) (words i) := by
  classical
  simp only [turboPositionCount22_eq_sum,
    Spin.prod_sum_eq_sum_prod, turboTupleCount22,
    turboOuterWord22_product_law, Spin.piPMF_apply,
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

#print axioms turboOuterWord22_product_law
#print axioms turboTupleCount22_eq_prod

end Spin.Turbo
