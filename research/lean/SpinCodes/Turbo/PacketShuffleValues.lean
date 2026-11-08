import SpinCodes.Turbo.PacketRoute
import SpinCodes.Turbo.PacketSupportCount

/-!
# Uniform packet shuffling with values retained

The completed Structured proof counts the permutation fiber of a support.
For complete packet words, a target word and a permutation determine the
source labels uniquely. The fiber equivalence below connects the new
value-carrying law to that old support count.
-/

noncomputable section

namespace Spin.Turbo

open Spin.Structured.Routing Finset

abbrev PacketWord64 := Fin 64 → Packet4

def packetWordSupport (x : PacketWord64) : Finset (Fin 64) :=
  Finset.univ.filter (fun i => x i ≠ 0)

/-- The row-shuffle convention used by `routePackets`: the output at `i`
reads the source at `σ i`. -/
def shufflePacketWord (σ : Equiv.Perm (Fin 64))
    (x : PacketWord64) : PacketWord64 := x ∘ σ

theorem shufflePacketWord_support (σ : Equiv.Perm (Fin 64))
    (x : PacketWord64) :
    packetWordSupport (shufflePacketWord σ x) =
      shuffleSupport σ (packetWordSupport x) := by
  ext i
  change i ∈ Finset.univ.filter (fun j => x (σ j) ≠ 0) ↔
    i ∈ shuffleSupport σ (Finset.univ.filter (fun j => x j ≠ 0))
  rw [mem_shuffleSupport]
  simp only [Finset.mem_filter, Finset.mem_univ, true_and]

abbrev LabeledSupport (S : Finset (Fin 64)) :=
  {x : PacketWord64 // packetWordSupport x = S}

noncomputable instance (S : Finset (Fin 64)) : Fintype (LabeledSupport S) :=
  Fintype.ofFinite _

/-- A word with support exactly `S` is a choice of one nonzero four-bit
label at each member of `S`. -/
def packetWordSupportEquiv (S : Finset (Fin 64)) :
    LabeledSupport S ≃ (S → {p : Packet4 // p ≠ 0}) := by
  classical
  refine ⟨(fun x i => ⟨x.val i.val, ?_⟩),
    (fun f => ⟨(fun i => if h : i ∈ S then (f ⟨i, h⟩).val else 0), ?_⟩),
    ?_, ?_⟩
  · have hi : i.val ∈ packetWordSupport x.val := by
      rw [x.property]
      exact i.property
    exact (Finset.mem_filter.mp hi).2
  · ext i
    simp only [packetWordSupport, Finset.mem_filter,
      Finset.mem_univ, true_and]
    by_cases hi : i ∈ S
    · simp [hi, (f ⟨i, hi⟩).property]
    · simp [hi]
  · intro x
    apply Subtype.ext
    funext i
    by_cases hi : i ∈ S
    · simp [hi]
    · have hz : x.val i = 0 := by
        by_contra hn
        have hm : i ∈ packetWordSupport x.val :=
          Finset.mem_filter.mpr ⟨Finset.mem_univ _, hn⟩
        rw [x.property] at hm
        exact hi hm
      simp [hi, hz]
  · intro f
    funext i
    apply Subtype.ext
    simp

theorem packetWordSupport_card (S : Finset (Fin 64)) :
    Fintype.card (LabeledSupport S) = 15 ^ S.card := by
  rw [Fintype.card_congr (packetWordSupportEquiv S)]
  simp

noncomputable instance (S : Finset (Fin 64)) :
    Nonempty (LabeledSupport S) :=
  ⟨(packetWordSupportEquiv S).symm (fun _ => ⟨1, one_ne_zero⟩)⟩

abbrev SupportPermFiber (S T : Finset (Fin 64)) :=
  {σ : Equiv.Perm (Fin 64) // shuffleSupport σ S = T}

abbrev PacketShuffleFiber (S : Finset (Fin 64)) (y : PacketWord64) :=
  {p : LabeledSupport S × Equiv.Perm (Fin 64) //
    shufflePacketWord p.2 p.1.val = y}

/-- For each target packet word, the permutation determines its unique
source labeling. The probability calculation can therefore reuse the old
support-permutation fiber count without discarding packet values. -/
def packetShuffleFiberEquiv (S : Finset (Fin 64)) (y : PacketWord64) :
    PacketShuffleFiber S y ≃ SupportPermFiber S (packetWordSupport y) where
  toFun p := ⟨p.val.2, by
    have hs := shufflePacketWord_support p.val.2 p.val.1.val
    rw [p.val.1.property, p.property] at hs
    exact hs.symm⟩
  invFun σ := ⟨
    ((⟨shufflePacketWord σ.val.symm y, by
      calc
        packetWordSupport (shufflePacketWord σ.val.symm y) =
            shuffleSupport σ.val.symm (packetWordSupport y) :=
          shufflePacketWord_support _ _
        _ = shuffleSupport σ.val.symm (shuffleSupport σ.val S) :=
          congrArg _ σ.property.symm
        _ = S := by ext i; simp⟩ : LabeledSupport S), σ.val),
    by
      funext i
      simp [shufflePacketWord]⟩
  left_inv p := by
    apply Subtype.ext
    apply Prod.ext
    · apply Subtype.ext
      funext i
      have h := congrFun p.property (p.val.2.symm i)
      simpa [shufflePacketWord] using h.symm
    · rfl
  right_inv σ := by
    apply Subtype.ext
    rfl

/-- Sample active packet labels independently and uniformly, then apply an
independent uniform permutation to the complete packet word. -/
def packetShuffleLaw (S : Finset (Fin 64)) : FinPMF PacketWord64 :=
  ((FinPMF.uniform (LabeledSupport S)).prod
    (FinPMF.uniform (Equiv.Perm (Fin 64)))).map
      (fun p => shufflePacketWord p.2 p.1.val)

private theorem uniform_prod_packet (S : Finset (Fin 64)) :
    (FinPMF.uniform (LabeledSupport S)).prod
      (FinPMF.uniform (Equiv.Perm (Fin 64))) =
        FinPMF.uniform (LabeledSupport S × Equiv.Perm (Fin 64)) := by
  ext p
  change (1 / (Fintype.card (LabeledSupport S) : ℝ)) *
      (1 / (Fintype.card (Equiv.Perm (Fin 64)) : ℝ)) =
    1 / (Fintype.card (LabeledSupport S × Equiv.Perm (Fin 64)) : ℝ)
  rw [Fintype.card_prod]
  push_cast
  ring

/-- Complete packet words are uniform within their packet-weight layer.
This is the old support-shuffle fiber count together with the independent
nonzero labels; it does not assume a fresh coefficient per message. -/
theorem packetShuffleLaw_apply (S : Finset (Fin 64))
    (y : PacketWord64) :
    (packetShuffleLaw S).p y =
      if (packetWordSupport y).card = S.card then
        1 / ((Nat.choose 64 S.card : ℝ) * (15 : ℝ) ^ S.card)
      else 0 := by
  classical
  have hp : (packetShuffleLaw S).p y =
      (packetShuffleLaw S).prob (fun z => z = y) := by
    simp [FinPMF.prob, Finset.sum_filter]
  rw [hp, packetShuffleLaw, FinPMF.map_prob, uniform_prod_packet,
    FinPMF.uniform_prob]
  have hcard :
      (Finset.univ.filter (fun p : LabeledSupport S × Equiv.Perm (Fin 64) =>
        shufflePacketWord p.2 p.1.val = y)).card =
        Fintype.card (PacketShuffleFiber S y) := by
    rw [Fintype.card_subtype]
  rw [hcard, Fintype.card_congr (packetShuffleFiberEquiv S y),
    Fintype.card_subtype, Fintype.card_prod, packetWordSupport_card]
  by_cases ht : (packetWordSupport y).card = S.card
  · rw [if_pos ht]
    have htotal := shuffle_fiber_total S (packetWordSupport y) ht
    have htotalR :
        (Fintype.card (Equiv.Perm (Fin 64)) : ℝ) =
          ((Finset.univ.filter (fun σ : Equiv.Perm (Fin 64) =>
            shuffleSupport σ S = packetWordSupport y)).card : ℝ) *
              (Nat.choose 64 S.card : ℝ) := by exact_mod_cast htotal
    have hc : (0 : ℝ) < (Nat.choose 64 S.card : ℝ) := choose_card_pos S
    have hp : (0 : ℝ) < (Fintype.card (Equiv.Perm (Fin 64)) : ℝ) := by
      exact_mod_cast Fintype.card_pos
    have hl : (0 : ℝ) < (15 : ℝ) ^ S.card := by positivity
    push_cast
    field_simp
    nlinarith
  · rw [if_neg ht]
    have he : Finset.univ.filter (fun σ : Equiv.Perm (Fin 64) =>
        shuffleSupport σ S = packetWordSupport y) = ∅ := by
      apply Finset.filter_eq_empty_iff.mpr
      intro σ _ hs
      exact ht (by rw [← hs, shuffleSupport_card])
    rw [he]
    simp only [Finset.card_empty, Nat.cast_zero, zero_div]

#print axioms shufflePacketWord_support
#print axioms packetShuffleFiberEquiv
#print axioms packetWordSupport_card
#print axioms packetShuffleLaw_apply

end Spin.Turbo
