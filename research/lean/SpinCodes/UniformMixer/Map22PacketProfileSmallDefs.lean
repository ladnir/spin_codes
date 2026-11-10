import SpinCodes.UniformMixer.Map22PacketProfileDefs

/-! Small quadratic truth tables used by canonical packet restrictions. -/
noncomputable section
namespace Spin.UniformMixer

abbrev ProfileQuad3 := Fin 7 → ZMod 2
abbrev ProfileQuad2 := Fin 4 → ZMod 2

def profileQuad3Rows : List Nat := [0xff,0xaa,0xcc,0xf0,0x88,0xa0,0xc0]
def profileQuad2Rows : List Nat := [0xf,0xa,0xc,0x8]

def profileQuad3Eval (q : ProfileQuad3) (v : Fin 8) : ZMod 2 :=
  ∑ i : Fin 7, bitScalar ((profileQuad3Rows.getD i.val 0).testBit v.val) * q i

def profileQuad2Eval (q : ProfileQuad2) (v : Fin 4) : ZMod 2 :=
  ∑ i : Fin 4, bitScalar ((profileQuad2Rows.getD i.val 0).testBit v.val) * q i

def profileQuad3Weight (q : ProfileQuad3) : Nat :=
  ∑ v : Fin 8, if profileQuad3Eval q v = 1 then 1 else 0

def profileQuad2Weight (q : ProfileQuad2) : Nat :=
  ∑ v : Fin 4, if profileQuad2Eval q v = 1 then 1 else 0

def profileQuad4Restrict3 (q : ProfileQuad4) : ProfileQuad3 :=
  ![q 0,q 1,q 2,q 3,q 5,q 6,q 8]

def profileQuad4Restrict2 (q : ProfileQuad4) : ProfileQuad2 :=
  ![q 0,q 1,q 2,q 5]

/-- Actual packed coefficient coordinates, for all small binary dimensions. -/
def profilePackedBits (n q : Nat) : Fin n → ZMod 2 :=
  fun i => bitScalar (q.testBit i.val)

theorem profilePackedBits_injective (n : Nat) :
    Function.Injective (fun q : Fin (2^n) => profilePackedBits n q.val) := by
  intro q r h
  apply Fin.ext
  apply Nat.eq_of_testBit_eq
  intro i
  by_cases hi : i < n
  · exact bitScalar_injective (congrFun h ⟨i,hi⟩)
  · have hpow : 2^n ≤ 2^i := Nat.pow_le_pow_right (by decide) (by omega)
    rw [Nat.testBit_lt_two_pow (q.isLt.trans_le hpow),
      Nat.testBit_lt_two_pow (r.isLt.trans_le hpow)]

def profilePackedBitsEquiv (n : Nat) : Fin (2^n) ≃ (Fin n → ZMod 2) :=
  Equiv.ofBijective (fun q => profilePackedBits n q.val)
    ((Fintype.bijective_iff_injective_and_card _).mpr
      ⟨profilePackedBits_injective n, by simp [ZMod.card]⟩)

end Spin.UniformMixer
