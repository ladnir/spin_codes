import SpinCodes.Turbo.PhysicalSmallBirth22Packed

noncomputable section
namespace Spin.Turbo
open Spin.UniformMixer
open Finset

def smallBirthLabel22 (v : Fin 15) : Packet4 := profilePackedBits 4 (v.val + 1)

set_option maxRecDepth 10000 in
theorem smallBirthLabel22_nonzero : ∀ v, smallBirthLabel22 v ≠ 0 := by decide +kernel

set_option maxRecDepth 10000 in
theorem smallBirthLabel22_injective : Function.Injective smallBirthLabel22 := by decide +kernel

abbrev SmallBirthOne22 := Fin 16 × Fin 15
abbrev SmallBirthSupportTwo22 := {ab : Fin 16 × Fin 16 // ab.1 < ab.2}
abbrev SmallBirthTwo22 := SmallBirthSupportTwo22 × Fin 15 × Fin 15

def smallBirthOneWord22 (a : SmallBirthOne22) : PacketStep16 :=
  fun i => if i = a.1 then smallBirthLabel22 a.2 else 0

def smallBirthTwoWord22 (a : SmallBirthTwo22) : PacketStep16 :=
  fun i => if i = a.1.val.1 then smallBirthLabel22 a.2.1 else
    if i = a.1.val.2 then smallBirthLabel22 a.2.2 else 0

theorem smallBirthOneWord22_support (a : SmallBirthOne22) :
    packetStepSupport22 (smallBirthOneWord22 a) = {a.1} := by
  ext i
  simp only [packetStepSupport22, mem_filter, mem_univ, true_and, mem_singleton]
  simp only [smallBirthOneWord22]
  split_ifs with h
  · simp [h, smallBirthLabel22_nonzero]
  · simp [h]

theorem smallBirthTwoWord22_support (a : SmallBirthTwo22) :
    packetStepSupport22 (smallBirthTwoWord22 a) = {a.1.val.1, a.1.val.2} := by
  ext i
  simp only [packetStepSupport22, mem_filter, mem_univ, true_and, mem_insert, mem_singleton]
  simp only [smallBirthTwoWord22]
  split_ifs with h h'
  · simp [h, smallBirthLabel22_nonzero]
  · simp [h', smallBirthLabel22_nonzero]
  · simp [h, h']

def smallBirthOneLayer22 (a : SmallBirthOne22) : PacketStepLayer22 1 :=
  ⟨smallBirthOneWord22 a, by rw [smallBirthOneWord22_support]; simp⟩

def smallBirthTwoLayer22 (a : SmallBirthTwo22) : PacketStepLayer22 2 :=
  ⟨smallBirthTwoWord22 a, by
    rw [smallBirthTwoWord22_support]
    simp [ne_of_lt a.1.property]⟩

theorem smallBirthOneLayer22_injective : Function.Injective smallBirthOneLayer22 := by
  intro a b h
  have hw := congrArg Subtype.val h
  change smallBirthOneWord22 a = smallBirthOneWord22 b at hw
  have hs := congrArg packetStepSupport22 hw
  rw [smallBirthOneWord22_support, smallBirthOneWord22_support] at hs
  have hab : a.1 = b.1 := singleton_injective hs
  apply Prod.ext hab
  apply smallBirthLabel22_injective
  have hv := congrFun hw a.1
  simpa only [smallBirthOneWord22, hab, ite_true] using hv

theorem smallBirthTwoLayer22_injective : Function.Injective smallBirthTwoLayer22 := by
  intro a b h
  have hw := congrArg Subtype.val h
  change smallBirthTwoWord22 a = smallBirthTwoWord22 b at hw
  have hs := congrArg packetStepSupport22 hw
  rw [smallBirthTwoWord22_support, smallBirthTwoWord22_support] at hs
  have ha : a.1.val.1 = b.1.val.1 ∨ a.1.val.1 = b.1.val.2 := by
    have hm : a.1.val.1 ∈ ({b.1.val.1, b.1.val.2} : Finset (Fin 16)) := by rw [← hs]; simp
    simpa only [mem_insert, mem_singleton] using hm
  have hb : a.1.val.2 = b.1.val.1 ∨ a.1.val.2 = b.1.val.2 := by
    have hm : a.1.val.2 ∈ ({b.1.val.1, b.1.val.2} : Finset (Fin 16)) := by rw [← hs]; simp
    simpa only [mem_insert, mem_singleton] using hm
  have hal := a.1.property
  have hbl := b.1.property
  have ha' : a.1.val.1 = b.1.val.1 := by rcases ha with h | h; exact h; rcases hb with h' | h' <;> omega
  have hb' : a.1.val.2 = b.1.val.2 := by rcases hb with h | h; omega; exact h
  apply Prod.ext
  · apply Subtype.ext
    exact Prod.ext ha' hb'
  · apply Prod.ext
    · apply smallBirthLabel22_injective
      have hv := congrFun hw a.1.val.1
      simpa only [smallBirthTwoWord22, ha', ite_true] using hv
    · apply smallBirthLabel22_injective
      have hv := congrFun hw a.1.val.2
      have hne : a.1.val.2 ≠ a.1.val.1 := ne_of_gt hal
      have hne' : a.1.val.2 ≠ b.1.val.1 := by rw [← ha']; exact hne
      dsimp only [smallBirthTwoWord22] at hv
      rw [if_neg hne, if_pos rfl, if_neg hne', if_pos hb'] at hv
      exact hv

set_option maxRecDepth 100000 in
theorem smallBirthSupportTwo22_card : Fintype.card SmallBirthSupportTwo22 = 120 := by decide +kernel

def smallBirthOneEquiv22 : SmallBirthOne22 ≃ PacketStepLayer22 1 :=
  Equiv.ofBijective smallBirthOneLayer22 ((Fintype.bijective_iff_injective_and_card _).mpr
    ⟨smallBirthOneLayer22_injective, by rw [packetStepLayer_card22]; simp⟩)

def smallBirthTwoEquiv22 : SmallBirthTwo22 ≃ PacketStepLayer22 2 :=
  Equiv.ofBijective smallBirthTwoLayer22 ((Fintype.bijective_iff_injective_and_card _).mpr
    ⟨smallBirthTwoLayer22_injective, by
      rw [packetStepLayer_card22]
      simp only [Fintype.card_prod, Fintype.card_fin, smallBirthSupportTwo22_card]
      norm_num [Nat.choose]⟩)

end Spin.Turbo



