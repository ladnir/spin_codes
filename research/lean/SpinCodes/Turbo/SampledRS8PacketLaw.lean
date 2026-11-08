import SpinCodes.Turbo.RS8FieldMoment
import SpinCodes.Turbo.PacketLabelConditioning

/-!
# The complete sampled packet law of one RS8 group

The eight multipliers are independent coordinates of one sampled setup.
Their output on a fixed RS8 group message is a product law over complete
eight-packet columns. Zero outer columns remain zero, while each nonzero
column has the uniform nonzero value law proved for the field adjoint.
-/

noncomputable section

namespace Spin.Turbo

open Finset

def sampledPacketColumn (x : Bits32)
    (c : {c : Field32 // c ≠ 0}) : PacketColumn :=
  bits32PacketEquiv (coordinateAdjoint field32Coordinates c.val x)

def sampledRS8PacketColumns (points : Fin 8 → Field16)
    (rows : RS8Group Field16)
    (coeff : Fin 8 → {c : Field32 // c ≠ 0}) :
    Fin 8 → PacketColumn :=
  fun j => sampledPacketColumn (rs8Bits points rows j) (coeff j)

def sampledPacketColumnLaw (x : Bits32) : FinPMF PacketColumn :=
  (FinPMF.uniform {c : Field32 // c ≠ 0}).map (sampledPacketColumn x)

def sampledRS8PacketColumnsLaw (points : Fin 8 → Field16)
    (rows : RS8Group Field16) : FinPMF (Fin 8 → PacketColumn) :=
  (piPMF (fun _ : Fin 8 =>
    FinPMF.uniform {c : Field32 // c ≠ 0})).map
      (sampledRS8PacketColumns points rows)

/-- The complete group law factors by column. This reuses the old
`FiniteProductLaw.piPMF_map` rather than resampling for each message. -/
theorem sampledRS8PacketColumnsLaw_product
    (points : Fin 8 → Field16) (rows : RS8Group Field16) :
    sampledRS8PacketColumnsLaw points rows =
      piPMF (fun j : Fin 8 =>
        sampledPacketColumnLaw (rs8Bits points rows j)) := by
  exact piPMF_map (fun _ : Fin 8 =>
    FinPMF.uniform {c : Field32 // c ≠ 0})
      (fun j c => sampledPacketColumn (rs8Bits points rows j) c)

/-- An inactive RS8 column contributes the zero packet column for every
sampled coefficient, hence with probability one. -/
theorem sampledPacketColumn_zero
    (c : {c : Field32 // c ≠ 0}) :
    sampledPacketColumn 0 c = 0 := by
  funext i j
  simp [sampledPacketColumn, coordinateAdjoint, bits32PacketEquiv]

/-- For a nonzero RS8 column, the complete packet value is uniform over
all nonzero columns. -/
theorem sampledPacketColumnLaw_nonzero
    (x : Bits32) (hx : x ≠ 0) :
    sampledPacketColumnLaw x =
      (FinPMF.uniform NonzeroPacketColumn).map Subtype.val := by
  let e := (nonzeroAdjointEquiv field32Coordinates.toAddEquiv binaryPair
    binaryPair_nondegenerate binaryPair_zero_left
    (coordinateAdjoint field32Coordinates)
    (coordinateAdjoint_pairing field32Coordinates) x hx).trans
      nonzeroBits32PacketEquiv
  have he (c : {c : Field32 // c ≠ 0}) :
      (e c).val = sampledPacketColumn x c := rfl
  calc
    sampledPacketColumnLaw x =
        ((FinPMF.uniform {c : Field32 // c ≠ 0}).map e).map
          Subtype.val := by
      rw [FinPMF.map_comp]
      rfl
    _ = _ := by rw [field32_uniform_packet_column x hx]

/-- The pushforward of a uniform nonzero packet column has one atom of
mass `1/(2^32-1)` at each nonzero value and no mass at zero. -/
theorem uniformNonzeroPacket_value_apply (y : PacketColumn) :
    ((FinPMF.uniform NonzeroPacketColumn).map Subtype.val).p y =
      if y = 0 then 0 else 1 / (2 ^ 32 - 1 : ℝ) := by
  classical
  rw [FinPMF.map_p]
  by_cases hy : y = 0
  · subst y
    simp only [ite_true]
    apply Finset.sum_eq_zero
    intro p _
    simp [p.property]
  · rw [if_neg hy]
    have h (p : NonzeroPacketColumn) :
        (p.val = y) ↔ p = ⟨y, hy⟩ := by
      constructor
      · intro hp; exact Subtype.ext hp
      · intro hp; exact congrArg Subtype.val hp
    simp only [h]
    simp only [Finset.sum_ite_eq', Finset.mem_univ, if_pos]
    change 1 / (Fintype.card NonzeroPacketColumn : ℝ) = _
    rw [packetColumn_nonzero_card]
    norm_num

/-- Point probabilities for a randomized active RS column. -/
theorem sampledPacketColumnLaw_nonzero_apply
    (x : Bits32) (hx : x ≠ 0) (y : PacketColumn) :
    (sampledPacketColumnLaw x).p y =
      if y = 0 then 0 else 1 / (2 ^ 32 - 1 : ℝ) := by
  rw [sampledPacketColumnLaw_nonzero x hx]
  exact uniformNonzeroPacket_value_apply y

theorem sampledPacketColumnLaw_zero_apply (y : PacketColumn) :
    (sampledPacketColumnLaw 0).p y =
      if y = 0 then 1 else 0 := by
  classical
  rw [sampledPacketColumnLaw, FinPMF.map_p]
  simp only [sampledPacketColumn_zero]
  by_cases hy : y = 0
  · subst y
    simp only [ite_true]
    exact (FinPMF.uniform {c : Field32 // c ≠ 0}).total
  · have hz : (0 : PacketColumn) ≠ y := Ne.symm hy
    simp [hy, hz]

/-- Exact point probabilities for a fixed RS8 group message under its
single sampled eight-coordinate multiplier tuple. Inactive RS columns
must output zero; each active column takes any one nonzero 32-bit packet
value with equal probability. -/
theorem sampledRS8PacketColumnsLaw_apply
    (points : Fin 8 → Field16) (rows : RS8Group Field16)
    (y : Fin 8 → PacketColumn) :
    (sampledRS8PacketColumnsLaw points rows).p y =
      ∏ j : Fin 8,
        if j ∈ rs8Support points rows then
          (if y j = 0 then 0 else 1 / (2 ^ 32 - 1 : ℝ))
        else (if y j = 0 then 1 else 0) := by
  classical
  rw [sampledRS8PacketColumnsLaw_product, piPMF_apply]
  apply Finset.prod_congr rfl
  intro j _
  by_cases hj : j ∈ rs8Support points rows
  · rw [if_pos hj]
    exact sampledPacketColumnLaw_nonzero_apply
      (rs8Bits points rows j)
      ((rs8Bits_ne_zero_iff points rows j).mpr hj) (y j)
  · rw [if_neg hj]
    have hz : rs8Bits points rows j = 0 := by
      by_contra hn
      exact hj ((rs8Bits_ne_zero_iff points rows j).mp hn)
    rw [hz]
    exact sampledPacketColumnLaw_zero_apply (y j)

theorem packetColumnSupport_eq_empty_iff (x : PacketColumn) :
    packetColumnSupport x = ∅ ↔ x = 0 := by
  constructor
  · intro h
    funext i
    by_contra hn
    have hi : i ∈ packetColumnSupport x :=
      Finset.mem_filter.mpr ⟨Finset.mem_univ _, hn⟩
    rw [h] at hi
    simp at hi
  · intro h
    subst x
    simp [packetColumnSupport]

/-- Under the actual sampled eight-column law, all complete output words
with the same packet support in every column have identical probability.
This is the group-level label symmetry needed before its row shuffle. -/
theorem sampledRS8PacketColumnsLaw_fiber_uniform
    (points : Fin 8 → Field16) (rows : RS8Group Field16)
    (y z : Fin 8 → PacketColumn)
    (hs : ∀ j, packetColumnSupport (y j) = packetColumnSupport (z j)) :
    (sampledRS8PacketColumnsLaw points rows).p y =
      (sampledRS8PacketColumnsLaw points rows).p z := by
  rw [sampledRS8PacketColumnsLaw_apply,
    sampledRS8PacketColumnsLaw_apply]
  apply Finset.prod_congr rfl
  intro j _
  have hz : y j = 0 ↔ z j = 0 := by
    rw [← packetColumnSupport_eq_empty_iff,
      ← packetColumnSupport_eq_empty_iff, hs j]
  by_cases hy : y j = 0
  · have hz' : z j = 0 := hz.mp hy
    simp [hy, hz']
  · have hz' : z j ≠ 0 := fun h => hy (hz.mpr h)
    simp [hy, hz']

#print axioms sampledRS8PacketColumnsLaw_product
#print axioms sampledPacketColumn_zero
#print axioms sampledPacketColumnLaw_nonzero
#print axioms uniformNonzeroPacket_value_apply
#print axioms sampledPacketColumnLaw_nonzero_apply
#print axioms sampledPacketColumnLaw_zero_apply
#print axioms sampledRS8PacketColumnsLaw_apply
#print axioms sampledRS8PacketColumnsLaw_fiber_uniform

end Spin.Turbo
