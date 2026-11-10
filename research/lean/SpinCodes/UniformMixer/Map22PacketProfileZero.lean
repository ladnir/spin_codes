import SpinCodes.UniformMixer.Map22PacketProfileWitness
import SpinCodes.UniformMixer.Map22PacketProfileRestriction

noncomputable section
namespace Spin.UniformMixer
set_option maxHeartbeats 1000000

def profileMaskedWeight (mask : Nat) (q : ProfileQuad4) : Nat :=
  ∑ v : Fin 16, if mask.testBit v.val then (if profileQuad4Eval q v = 1 then 1 else 0) else 0

def profileCanonicalWeight (d : Nat) (q : ProfileQuad4) : Nat :=
  ∑ v : Fin 16, if v.val < 2^d then (if profileQuad4Eval q v = 1 then 1 else 0) else 0

theorem profileZeroWitness_card (w : ProfileZeroWitness) (hw : w.check = true) :
    (Finset.univ.filter (fun v : Fin 16 => w.zeroMask.testBit v.val = true)).card =
      2^w.dimension := by
  rw [Finset.card_filter, ← Equiv.sum_comp (w.pointEquiv hw)]
  have he : (∑ v : Fin 16,
      if w.zeroMask.testBit (w.pointEquiv hw v).val = true then 1 else 0) =
      ∑ v : Fin 16, if v.val < 2^w.dimension then 1 else 0 := by
    apply Finset.sum_congr rfl
    intro v _
    simp only [w.zero_alignment hw v]
  rw [he]
  exact (show ∀ k : Fin 3,
    (∑ v : Fin 16, if v.val < 2^(k.val+2) then 1 else 0) = 2^(k.val+2) by decide) w.kind

theorem profileMaskedWeight_transform (w : ProfileZeroWitness) (hw : w.check = true)
    (q : ProfileQuad4) :
    profileMaskedWeight w.zeroMask q = profileCanonicalWeight w.dimension (w.transform q) := by
  unfold profileMaskedWeight profileCanonicalWeight
  rw [← Equiv.sum_comp (w.pointEquiv hw)]
  apply Finset.sum_congr rfl
  intro v _
  rw [w.transform_eval hw]
  have hz := w.zero_alignment hw v
  by_cases h : v.val < 2^w.dimension
  · rw [if_pos h, (hz.mpr h)]
    rfl
  · have hh : w.zeroMask.testBit (w.pointEquiv hw v).val = false := by
      cases hb : w.zeroMask.testBit (w.pointEquiv hw v).val
      · rfl
      · exact False.elim (h (hz.mp hb))
    rw [if_neg h, hh]
    rfl

theorem profileCanonicalWeight_two (q : ProfileQuad4) :
    profileCanonicalWeight 2 q = profileQuad2Weight (profileQuad4Restrict2 q) := by
  simp only [profileCanonicalWeight, profileQuad2Weight, Fin.sum_univ_succ]
  simp [profileQuad4Eval, profileQuad2Eval,
    profileQuad4Restrict2, profileQuad4Rows, profileQuad2Rows,
    Fin.sum_univ_succ, bitScalar, Nat.testBit_succ]

theorem profileCanonicalWeight_three (q : ProfileQuad4) :
    profileCanonicalWeight 3 q = profileQuad3Weight (profileQuad4Restrict3 q) := by
  simp only [profileCanonicalWeight, profileQuad3Weight, Fin.sum_univ_succ]
  simp [profileQuad4Eval, profileQuad3Eval,
    profileQuad4Restrict3, profileQuad4Rows, profileQuad3Rows,
    Fin.sum_univ_succ, bitScalar, Nat.testBit_succ]

theorem profileCanonicalWeight_four (q : ProfileQuad4) :
    profileCanonicalWeight 4 q = profileQuad4Weight q := by
  apply Finset.sum_congr rfl
  intro v _
  simp only [show v.val < 2^4 from v.isLt, if_true]

theorem profileMaskedWeight_sum {M : Type*} [AddCommMonoid M]
    (w : ProfileZeroWitness) (hw : w.check = true) (f : Nat → M) :
    (∑ q : ProfileQuad4, f (profileMaskedWeight w.zeroMask q)) =
      ∑ q : ProfileQuad4, f (profileCanonicalWeight w.dimension q) := by
  simp_rw [profileMaskedWeight_transform w hw]
  change (∑ q : ProfileQuad4, (fun r => f (profileCanonicalWeight w.dimension r))
      (w.quadraticEquiv hw q)) = _
  exact Equiv.sum_comp (w.quadraticEquiv hw) (fun r => f (profileCanonicalWeight w.dimension r))

theorem profileMaskedWeight_sum_two {M : Type*} [AddCommMonoid M]
    (w : ProfileZeroWitness) (hw : w.check = true) (hk : w.kind = 0) (f : Nat → M) :
    (∑ q : ProfileQuad4, f (profileMaskedWeight w.zeroMask q)) =
      128 • ∑ q : ProfileQuad2, f (profileQuad2Weight q) := by
  rw [profileMaskedWeight_sum w hw]
  have hd : w.dimension = 2 := by simp [ProfileZeroWitness.dimension, hk]
  simp_rw [hd, profileCanonicalWeight_two]
  exact profileQuad4_restrict2_sum (fun q => f (profileQuad2Weight q))

theorem profileMaskedWeight_sum_three {M : Type*} [AddCommMonoid M]
    (w : ProfileZeroWitness) (hw : w.check = true) (hk : w.kind = 1) (f : Nat → M) :
    (∑ q : ProfileQuad4, f (profileMaskedWeight w.zeroMask q)) =
      16 • ∑ q : ProfileQuad3, f (profileQuad3Weight q) := by
  rw [profileMaskedWeight_sum w hw]
  have hd : w.dimension = 3 := by simp [ProfileZeroWitness.dimension, hk]
  simp_rw [hd, profileCanonicalWeight_three]
  exact profileQuad4_restrict3_sum (fun q => f (profileQuad3Weight q))

theorem profileMaskedWeight_sum_four {M : Type*} [AddCommMonoid M]
    (w : ProfileZeroWitness) (hw : w.check = true) (hk : w.kind = 2) (f : Nat → M) :
    (∑ q : ProfileQuad4, f (profileMaskedWeight w.zeroMask q)) =
      ∑ q : ProfileQuad4, f (profileQuad4Weight q) := by
  rw [profileMaskedWeight_sum w hw]
  have hd : w.dimension = 4 := by simp [ProfileZeroWitness.dimension, hk]
  simp_rw [hd, profileCanonicalWeight_four]

#print axioms profileMaskedWeight_transform
#print axioms profileCanonicalWeight_three
end Spin.UniformMixer
