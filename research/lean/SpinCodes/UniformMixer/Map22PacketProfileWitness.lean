import SpinCodes.UniformMixer.Map22PacketProfileSmallDefs

/-!
# Small checked changes of coordinates for packet-profile zero sets

A witness is a permutation of only sixteen packet positions and eleven
basis images in the four-variable quadratic space. Its checker proves that
the permutation preserves the whole quadratic space, not just a sampled word.
-/
noncomputable section
namespace Spin.UniformMixer
open Spin.Structured

structure ProfileZeroWitness where
  kind : Fin 3
  zeroMask : Nat
  permCode : Nat
  invCode : Nat
  coeffs : List Nat

def profilePerm16 (code : Nat) (v : Fin 16) : Fin 16 :=
  ⟨(code >>> (4*v.val)) % 16, Nat.mod_lt _ (by decide)⟩

namespace ProfileZeroWitness

def dimension (w : ProfileZeroWitness) : Nat := w.kind.val + 2

def check (w : ProfileZeroWitness) : Bool := decide (
  w.coeffs.length = 11 ∧
  (∀ v : Fin 16, profilePerm16 w.invCode (profilePerm16 w.permCode v) = v) ∧
  (∀ v : Fin 16, profilePerm16 w.permCode (profilePerm16 w.invCode v) = v) ∧
  (∀ v : Fin 16, w.zeroMask.testBit (profilePerm16 w.permCode v).val =
    decide (v.val < 2^w.dimension)) ∧
  (∀ j : Fin 11, ∀ v : Fin 16,
    (PackedMap.eval profileQuad4Rows (w.coeffs.getD j.val 0)).testBit v.val =
      (profileQuad4Rows.getD j.val 0).testBit (profilePerm16 w.permCode v).val))

def pointEquiv (w : ProfileZeroWitness) (hw : w.check = true) : Fin 16 ≃ Fin 16 where
  toFun := profilePerm16 w.permCode
  invFun := profilePerm16 w.invCode
  left_inv := (of_decide_eq_true hw).2.1
  right_inv := (of_decide_eq_true hw).2.2.1

def transform (w : ProfileZeroWitness) (q : ProfileQuad4) : ProfileQuad4 :=
  fun i => ∑ j : Fin 11, bitScalar ((w.coeffs.getD j.val 0).testBit i.val) * q j

theorem transform_eval (w : ProfileZeroWitness) (hw : w.check = true)
    (q : ProfileQuad4) (v : Fin 16) :
    profileQuad4Eval (w.transform q) v = profileQuad4Eval q (w.pointEquiv hw v) := by
  have hb := (of_decide_eq_true hw).2.2.2.2
  unfold profileQuad4Eval transform
  simp only [Finset.mul_sum]
  rw [Finset.sum_comm]
  apply Finset.sum_congr rfl
  intro j _
  calc
    (∑ i : Fin 11, bitScalar ((profileQuad4Rows.getD i.val 0).testBit v.val) *
        (bitScalar ((w.coeffs.getD j.val 0).testBit i.val) * q j)) =
      (∑ i : Fin 11, bitScalar ((profileQuad4Rows.getD i.val 0).testBit v.val) *
        bitScalar ((w.coeffs.getD j.val 0).testBit i.val)) * q j := by
      rw [Finset.sum_mul]
      apply Finset.sum_congr rfl
      intro i _; ring
    _ = _ := by
      have hp : (∑ i : Fin 11,
          bitScalar ((profileQuad4Rows.getD i.val 0).testBit v.val) *
            bitScalar ((w.coeffs.getD j.val 0).testBit i.val)) =
          bitScalar ((PackedMap.eval profileQuad4Rows (w.coeffs.getD j.val 0)).testBit v.val) :=
        (packed_eval_coordinate profileQuad4Rows (w.coeffs.getD j.val 0) v.val).symm
      rw [hp, hb j v]
      rfl

theorem transform_injective (w : ProfileZeroWitness) (hw : w.check = true) :
    Function.Injective w.transform := by
  intro q r h
  apply profileQuad4Eval_injective
  funext v
  have he := congrArg (fun t => profileQuad4Eval t ((w.pointEquiv hw).symm v)) h
  simpa only [w.transform_eval hw, Equiv.apply_symm_apply] using he

def quadraticEquiv (w : ProfileZeroWitness) (hw : w.check = true) :
    ProfileQuad4 ≃ ProfileQuad4 :=
  Equiv.ofBijective w.transform
    ((Fintype.bijective_iff_injective_and_card _).mpr ⟨w.transform_injective hw, rfl⟩)

theorem zero_alignment (w : ProfileZeroWitness) (hw : w.check = true) (v : Fin 16) :
    w.zeroMask.testBit (w.pointEquiv hw v).val = true ↔ v.val < 2^w.dimension := by
  rw [show w.zeroMask.testBit (w.pointEquiv hw v).val =
      decide (v.val < 2^w.dimension) from (of_decide_eq_true hw).2.2.2.1 v]
  simp

#print axioms transform_eval
#print axioms quadraticEquiv
end ProfileZeroWitness
end Spin.UniformMixer
