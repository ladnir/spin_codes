import SpinCodes.Turbo.PhysicalSmallBirth22Inputs
import SpinCodes.Turbo.PhysicalSmallBirth22NumericData

namespace Spin.Turbo
open Spin.UniformMixer
open Spin.Structured

def smallBirthOnePacked22 (a : SmallBirthOne22) : Nat :=
  (a.2.val + 1) <<< (4 * a.1.val)

def smallBirthTwoPacked22 (a : SmallBirthTwo22) : Nat :=
  smallBirthOnePacked22 (a.1.val.1, a.2.1) ^^^
    smallBirthOnePacked22 (a.1.val.2, a.2.2)

theorem smallBirthOnePacked22_input (a : SmallBirthOne22) :
    smallBirthPackedInput22 (smallBirthOnePacked22 a) =
      packetStepInput22 (smallBirthOneWord22 a) := by
  funext j
  obtain ⟨⟨slot, bit⟩, rfl⟩ := (finProdFinEquiv : Fin 16 × Fin 4 ≃ Fin 64).surjective j
  simp only [packetStepInput22, Equiv.symm_apply_apply, smallBirthOneWord22]
  change bitScalar (((a.2.val+1) <<< (4*a.1.val)).testBit (bit.val+4*slot.val)) = _
  rw [Nat.testBit_shiftLeft]
  by_cases hs : slot = a.1
  · subst slot
    have hle : 4*a.1.val ≤ bit.val+4*a.1.val := by omega
    have he : bit.val+4*a.1.val-4*a.1.val = bit.val := by omega
    simp [hle, he, smallBirthLabel22, profilePackedBits]
  · rw [if_neg hs]
    by_cases hlt : slot.val < a.1.val
    · have hlt' : ¬bit.val+4*slot.val ≥ 4*a.1.val := by omega
      simp [hlt', bitScalar]
    · have hidx : 4 ≤ bit.val+4*slot.val-4*a.1.val := by
        have hne : slot.val ≠ a.1.val := fun h => hs (Fin.ext h)
        omega
      have hv : a.2.val+1 < 2^(bit.val+4*slot.val-4*a.1.val) := by
        have hp : 2^4 ≤ 2^(bit.val+4*slot.val-4*a.1.val) :=
          Nat.pow_le_pow_right (by decide) hidx
        have ha := a.2.isLt
        norm_num at hp
        omega
      rw [Nat.testBit_lt_two_pow hv]
      simp [bitScalar]
theorem smallBirthTwoWord22_add (a : SmallBirthTwo22) :
    smallBirthTwoWord22 a = smallBirthOneWord22 (a.1.val.1, a.2.1) +
      smallBirthOneWord22 (a.1.val.2, a.2.2) := by
  funext i
  simp only [smallBirthTwoWord22, smallBirthOneWord22, Pi.add_apply]
  have hne : a.1.val.1 ≠ a.1.val.2 := ne_of_lt a.1.property
  split_ifs with ha hb hb <;> simp_all

theorem smallBirthTwoPacked22_input (a : SmallBirthTwo22) :
    smallBirthPackedInput22 (smallBirthTwoPacked22 a) =
      packetStepInput22 (smallBirthTwoWord22 a) := by
  rw [smallBirthTwoPacked22, smallBirthPackedInput22_xor,
    smallBirthOnePacked22_input, smallBirthOnePacked22_input, smallBirthTwoWord22_add]
  rfl

theorem smallBirthClassWeight22_actual (r : State22) :
    smallBirthClassWeight22 (coordinateWeight coordMap r) = physicalBirthIndex22 r := by
  by_cases hr : r = 0
  · subst r
    simp [smallBirthClassWeight22, physicalBirthIndex22, coordinateWeight]
  · have hw : coordinateWeight coordMap r ≠ 0 := fun h => hr (productionWeight_zero_only r h)
    simp only [smallBirthClassWeight22, physicalBirthIndex22, hr, hw, ite_false]

def smallBirthPacketWeight22 (p : Packet4) : Nat :=
  ∑ i : Fin 4, if p i = 1 then 1 else 0

theorem smallBirthPacketWeight22_zero : smallBirthPacketWeight22 0 = 0 := by
  simp [smallBirthPacketWeight22]

theorem smallBirthStepWeight22 (x : PacketStep16) :
    stepWeight (packetStepInput22 x) 0 = ∑ slot, smallBirthPacketWeight22 (x slot) := by
  unfold stepWeight packetStepInput22 smallBirthPacketWeight22
  simp only [map_zero, zero_add]
  rw [← Equiv.sum_comp (finProdFinEquiv : Fin 16 × Fin 4 ≃ Fin 64), Fintype.sum_prod_type]
  simp only [Equiv.symm_apply_apply]

theorem smallBirthOneWord22_weight (a : SmallBirthOne22) :
    stepWeight (packetStepInput22 (smallBirthOneWord22 a)) 0 =
      smallBirthPacketWeight22 (smallBirthLabel22 a.2) := by
  rw [smallBirthStepWeight22]
  simp [smallBirthOneWord22, apply_ite, smallBirthPacketWeight22_zero]

theorem smallBirthTwoWord22_weight (a : SmallBirthTwo22) :
    stepWeight (packetStepInput22 (smallBirthTwoWord22 a)) 0 =
      smallBirthPacketWeight22 (smallBirthLabel22 a.2.1) +
        smallBirthPacketWeight22 (smallBirthLabel22 a.2.2) := by
  rw [smallBirthStepWeight22]
  have hne : a.1.val.1 ≠ a.1.val.2 := ne_of_lt a.1.property
  have ht (i : Fin 16) : smallBirthPacketWeight22 (smallBirthTwoWord22 a i) =
      (if i = a.1.val.1 then smallBirthPacketWeight22 (smallBirthLabel22 a.2.1) else 0) +
      (if i = a.1.val.2 then smallBirthPacketWeight22 (smallBirthLabel22 a.2.2) else 0) := by
    by_cases h1 : i = a.1.val.1
    · subst i
      simp [smallBirthTwoWord22, hne]
    · by_cases h2 : i = a.1.val.2
      · simp [smallBirthTwoWord22, h1, h2, hne.symm]
      · simp [smallBirthTwoWord22, h1, h2, smallBirthPacketWeight22_zero]
  simp_rw [ht]
  rw [Finset.sum_add_distrib]
  simp

set_option maxRecDepth 10000 in
theorem smallBirthLabel22_weight : ∀ v : Fin 15,
    smallBirthPacketWeight22 (smallBirthLabel22 v) = MapSpectrum.fastWeight 4 (v.val + 1) := by
  decide +kernel
end Spin.Turbo





