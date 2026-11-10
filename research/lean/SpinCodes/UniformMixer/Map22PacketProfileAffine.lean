import SpinCodes.UniformMixer.Map22PacketProfileSmallDefs

noncomputable section
namespace Spin.UniformMixer
open Spin.Structured

def profileAffinePackedEquiv : Fin 32 ≃ ProfileAffine4 := profilePackedBitsEquiv 5

def profileAffinePairEquiv : Fin 1024 ≃ ProfileAffine4 × ProfileAffine4 :=
  ((finProdFinEquiv : Fin 32 × Fin 32 ≃ Fin 1024).symm.trans (Equiv.prodComm _ _)).trans
    (Equiv.prodCongr profileAffinePackedEquiv profileAffinePackedEquiv)

theorem profileAffinePacked_eval (p : Fin 32) (v : Fin 16) :
    profileAffine4Eval (profileAffinePackedEquiv p) v =
      bitScalar ((PackedMap.eval (profileQuad4Rows.take 5) p.val).testBit v.val) := by
  rw [packed_eval_coordinate]
  rfl

theorem profileAffinePair_fst (p : Fin 1024) :
    (profileAffinePairEquiv p).1 = profileAffinePackedEquiv ⟨p.val % 32, Nat.mod_lt _ (by decide)⟩ := rfl

theorem profileAffinePair_snd (p : Fin 1024) :
    (profileAffinePairEquiv p).2 = profileAffinePackedEquiv ⟨p.val / 32, by omega⟩ := rfl

theorem bitScalar_eq_zero (b : Bool) : bitScalar b = 0 ↔ b = false := by
  cases b <;> simp [bitScalar]

#print axioms profileAffinePacked_eval
end Spin.UniformMixer
