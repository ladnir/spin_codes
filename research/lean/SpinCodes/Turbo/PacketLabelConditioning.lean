import SpinCodes.Turbo.ColumnAdjointField

/-!
# Packet labels conditional on their support

The sampled field adjoint produces a uniform nonzero eight-packet column.
This file counts the exact support fiber and derives the conditional law
of its nonzero four-bit labels. It is a value-level bridge; the old
Structured proof counted only supports under a uniform permutation.
-/

noncomputable section

namespace Spin.Turbo

open Finset

abbrev NonzeroPacketColumn := {p : PacketColumn // p ≠ 0}

abbrev NonzeroColumnSupport (s : Finset (Fin 8)) :=
  {p : NonzeroPacketColumn // packetColumnSupport p.val = s}

noncomputable instance (s : Finset (Fin 8)) :
    Fintype (NonzeroColumnSupport s) := Fintype.ofFinite _

/-- For nonempty `s`, the nonzero-column support fiber is exactly the
unrestricted support fiber already counted by `packetColumnSupport_card`. -/
def nonzeroColumnSupportEquiv (s : Finset (Fin 8)) (hs : s.Nonempty) :
    NonzeroColumnSupport s ≃
      {p : PacketColumn // packetColumnSupport p = s} where
  toFun p := ⟨p.val.val, p.property⟩
  invFun p := ⟨⟨p.val, by
    intro hz
    have hzero : packetColumnSupport (0 : PacketColumn) = ∅ := by
      simp [packetColumnSupport]
    have he : s = ∅ := by
      rw [← p.property, hz]
      exact hzero
    exact hs.ne_empty he⟩, p.property⟩
  left_inv p := by
    apply Subtype.ext
    apply Subtype.ext
    rfl
  right_inv p := by
    apply Subtype.ext
    rfl

theorem nonzeroColumnSupport_card (s : Finset (Fin 8))
    (hs : s.Nonempty) :
    Fintype.card (NonzeroColumnSupport s) = 15 ^ s.card := by
  rw [Fintype.card_congr (nonzeroColumnSupportEquiv s hs)]
  exact packetColumnSupport_card s

/-- The exact probability of one nonempty packet support in a uniform
nonzero 32-bit column. -/
theorem uniform_nonzero_column_support_prob (s : Finset (Fin 8))
    (hs : s.Nonempty) :
    (FinPMF.uniform NonzeroPacketColumn).prob
      (fun p => packetColumnSupport p.val = s) =
      (15 : ℝ) ^ s.card / (2 ^ 32 - 1 : ℝ) := by
  classical
  rw [FinPMF.uniform_prob]
  have hcard :
      (Finset.univ.filter (fun p : NonzeroPacketColumn =>
        packetColumnSupport p.val = s)).card =
        Fintype.card (NonzeroColumnSupport s) := by
    rw [Fintype.card_subtype]
  rw [hcard, nonzeroColumnSupport_card s hs, packetColumn_nonzero_card]
  norm_num

theorem uniform_nonzero_column_support_pos (s : Finset (Fin 8))
    (hs : s.Nonempty) :
    0 < (FinPMF.uniform NonzeroPacketColumn).prob
      (fun p => packetColumnSupport p.val = s) := by
  rw [uniform_nonzero_column_support_prob s hs]
  positivity

/-- Given its nonempty support, each compatible nonzero column has the
same probability `15^(-|s|)`. This is the uniform-label premise required
by the value-carrying packet shuffle. -/
theorem uniform_nonzero_column_conditional_apply
    (s : Finset (Fin 8)) (hs : s.Nonempty)
    (p : NonzeroPacketColumn) :
    ((FinPMF.uniform NonzeroPacketColumn).condition
      (fun q => packetColumnSupport q.val = s)
      (uniform_nonzero_column_support_pos s hs)).p p =
        if packetColumnSupport p.val = s then
          1 / (15 : ℝ) ^ s.card else 0 := by
  classical
  rw [FinPMF.condition_p_apply,
    uniform_nonzero_column_support_prob s hs]
  by_cases hp : packetColumnSupport p.val = s
  · simp only [if_pos hp]
    change (1 / (Fintype.card NonzeroPacketColumn : ℝ)) /
      ((15 : ℝ) ^ s.card / (2 ^ 32 - 1)) =
        1 / (15 : ℝ) ^ s.card
    rw [packetColumn_nonzero_card]
    norm_num
    field_simp
  · simp only [if_neg hp]

#print axioms nonzeroColumnSupport_card
#print axioms uniform_nonzero_column_support_prob
#print axioms uniform_nonzero_column_conditional_apply

end Spin.Turbo
