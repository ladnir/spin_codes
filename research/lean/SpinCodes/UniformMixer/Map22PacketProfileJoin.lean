import SpinCodes.UniformMixer.Map22PacketProfileZero
import SpinCodes.UniformMixer.Map22PacketProfileHistogram
import SpinCodes.UniformMixer.Map22PacketProfileShapes
import SpinCodes.UniformMixer.Map22PacketProfileAffine

noncomputable section
namespace Spin.UniformMixer
set_option maxHeartbeats 1000000

theorem profileSplit6_odd_sum {M : Type*} [AddCommMonoid M]
    (b c : ProfileAffine4) (f : (Fin 5 → Nat) → M) :
    (∑ q : ProfileQuad4, f (profileOfNibbles (profileSplit6Eval (1,b,c,q)))) =
      ∑ w : Fin 17, profileQuad4WeightCounts.getD w.val 0 • f (profileOddShape w.val) := by
  simp_rw [profileSplit6_odd_shape]
  calc
    _ = ∑ q : ProfileQuad4, f (profileOddShape (profileQuad4Weight q)) := by
      change (∑ q : ProfileQuad4, (fun r => f (profileOddShape (profileQuad4Weight r)))
        (profileQuad4Shift (profileAffineProduct b c) q)) = _
      exact Equiv.sum_comp (profileQuad4Shift (profileAffineProduct b c))
        (fun r => f (profileOddShape (profileQuad4Weight r)))
    _ = _ := profileQuad4_weight_sum (fun r => f (profileOddShape r))

theorem profileSplit6_even_mask_shape (b c : ProfileAffine4) (q : ProfileQuad4)
    (mask : Nat) (hm : ∀ v : Fin 16, mask.testBit v.val = true ↔
      profileAffine4Eval b v = 0 ∧ profileAffine4Eval c v = 0) :
    profileOfNibbles (profileSplit6Eval (0,b,c,q)) =
      profileEvenShape
        (Finset.univ.filter (fun v : Fin 16 => mask.testBit v.val = true)).card
        (profileMaskedWeight mask q) := by
  have hs := profileOfNibbles_even_shape (profileSplit6Eval (0,b,c,q))
    (profileQuad4Eval q) (Finset.univ.filter (fun v : Fin 16 => mask.testBit v.val = true))
    (by intro v; simpa only [Finset.mem_filter, Finset.mem_univ, true_and, hm v]
      using profileEvenNibbleWeight b c q v)
  have he : ((Finset.univ.filter (fun v : Fin 16 => mask.testBit v.val = true)).filter
      (fun v => profileQuad4Eval q v = 1)).card = profileMaskedWeight mask q := by
    rw [Finset.card_filter]
    unfold profileMaskedWeight
    rw [Finset.sum_filter]
  exact hs.trans (congrArg (profileEvenShape _) he)

theorem profileSplit6_even_witness_shape (b c : ProfileAffine4) (q : ProfileQuad4)
    (w : ProfileZeroWitness) (hw : w.check = true)
    (hm : ∀ v : Fin 16, w.zeroMask.testBit v.val = true ↔
      profileAffine4Eval b v = 0 ∧ profileAffine4Eval c v = 0) :
    profileOfNibbles (profileSplit6Eval (0,b,c,q)) =
      profileEvenShape (2^w.dimension) (profileMaskedWeight w.zeroMask q) := by
  rw [profileSplit6_even_mask_shape b c q w.zeroMask hm, profileZeroWitness_card w hw]

theorem profileSplit6_even_witness_sum_two {M : Type*} [AddCommMonoid M]
    (b c : ProfileAffine4) (w : ProfileZeroWitness) (hw : w.check = true)
    (hk : w.kind = 0) (f : (Fin 5 → Nat) → M)
    (hm : ∀ v : Fin 16, w.zeroMask.testBit v.val = true ↔
      profileAffine4Eval b v = 0 ∧ profileAffine4Eval c v = 0) :
    (∑ q : ProfileQuad4, f (profileOfNibbles (profileSplit6Eval (0,b,c,q)))) =
      128 • ∑ a : Fin 5, profileQuad2WeightCounts.getD a.val 0 • f (profileEvenShape 4 a.val) := by
  simp_rw [profileSplit6_even_witness_shape b c _ w hw hm]
  have hd : 2^w.dimension = 4 := by simp [ProfileZeroWitness.dimension, hk]
  rw [hd, profileMaskedWeight_sum_two w hw hk (fun r => f (profileEvenShape 4 r)),
    profileQuad2_weight_sum (fun r => f (profileEvenShape 4 r))]

theorem profileSplit6_even_witness_sum_three {M : Type*} [AddCommMonoid M]
    (b c : ProfileAffine4) (w : ProfileZeroWitness) (hw : w.check = true)
    (hk : w.kind = 1) (f : (Fin 5 → Nat) → M)
    (hm : ∀ v : Fin 16, w.zeroMask.testBit v.val = true ↔
      profileAffine4Eval b v = 0 ∧ profileAffine4Eval c v = 0) :
    (∑ q : ProfileQuad4, f (profileOfNibbles (profileSplit6Eval (0,b,c,q)))) =
      16 • ∑ a : Fin 9, profileQuad3WeightCounts.getD a.val 0 • f (profileEvenShape 8 a.val) := by
  simp_rw [profileSplit6_even_witness_shape b c _ w hw hm]
  have hd : 2^w.dimension = 8 := by simp [ProfileZeroWitness.dimension, hk]
  rw [hd, profileMaskedWeight_sum_three w hw hk (fun r => f (profileEvenShape 8 r)),
    profileQuad3_weight_sum (fun r => f (profileEvenShape 8 r))]

theorem profileSplit6_even_witness_sum_four {M : Type*} [AddCommMonoid M]
    (b c : ProfileAffine4) (w : ProfileZeroWitness) (hw : w.check = true)
    (hk : w.kind = 2) (f : (Fin 5 → Nat) → M)
    (hm : ∀ v : Fin 16, w.zeroMask.testBit v.val = true ↔
      profileAffine4Eval b v = 0 ∧ profileAffine4Eval c v = 0) :
    (∑ q : ProfileQuad4, f (profileOfNibbles (profileSplit6Eval (0,b,c,q)))) =
      ∑ a : Fin 17, profileQuad4WeightCounts.getD a.val 0 • f (profileEvenShape 16 a.val) := by
  simp_rw [profileSplit6_even_witness_shape b c _ w hw hm]
  have hd : 2^w.dimension = 16 := by simp [ProfileZeroWitness.dimension, hk]
  rw [hd, profileMaskedWeight_sum_four w hw hk (fun r => f (profileEvenShape 16 r)),
    profileQuad4_weight_sum (fun r => f (profileEvenShape 16 r))]

#print axioms profileSplit6_odd_sum
#print axioms profileSplit6_even_witness_sum_three

end Spin.UniformMixer
