import SpinCodes.UniformMixer.Map22PacketProfileClasses
import SpinCodes.UniformMixer.Map22PacketProfileJoin
import SpinCodes.UniformMixer.Map22PacketProfileBridge

noncomputable section
namespace Spin.UniformMixer
set_option maxHeartbeats 1000000

def profileEvenClassSum {M : Type*} [AddCommMonoid M]
    (f : (Fin 5 → Nat) → M) (k : Nat) : M :=
  if k = 0 then 2048 • f (profileEvenShape 0 0) else
  if k = 1 then 128 • ∑ a : Fin 5, profileQuad2WeightCounts.getD a.val 0 • f (profileEvenShape 4 a.val) else
  if k = 2 then 16 • ∑ a : Fin 9, profileQuad3WeightCounts.getD a.val 0 • f (profileEvenShape 8 a.val) else
    ∑ a : Fin 17, profileQuad4WeightCounts.getD a.val 0 • f (profileEvenShape 16 a.val)

theorem profileSplit6_even_empty_shape (b c : ProfileAffine4) (q : ProfileQuad4)
    (h : ∀ v, ¬ (profileAffine4Eval b v = 0 ∧ profileAffine4Eval c v = 0)) :
    profileOfNibbles (profileSplit6Eval (0,b,c,q)) = profileEvenShape 0 0 := by
  have hs := profileOfNibbles_even_shape (profileSplit6Eval (0,b,c,q))
    (profileQuad4Eval q) ∅ (by intro v; simpa [h v] using profileEvenNibbleWeight b c q v)
  simpa only [Finset.card_empty, Finset.filter_empty] using hs

theorem profileSplit6_even_pair_sum {M : Type*} [AddCommMonoid M]
    (p : Fin 1024) (f : (Fin 5 → Nat) → M) :
    (∑ q : ProfileQuad4, f (profileOfNibbles (profileSplit6Eval
      (0,(profileAffinePairEquiv p).1,(profileAffinePairEquiv p).2,q)))) =
      profileEvenClassSum f (profileZeroClass p) := by
  by_cases hp : profileZeroIndex p = 0
  · simp_rw [profileSplit6_even_empty_shape _ _ _ (profileClassZero_empty p hp)]
    simp only [profileZeroClass, hp, ite_true, profileEvenClassSum, Finset.sum_const,
      Finset.card_univ, Fintype.card_fun, Fintype.card_fin, ZMod.card]
    rfl
  · have hw := profileClassWitness_checked p hp
    have hm := profileClassWitness_alignment p hp
    have hk : (profileClassWitness p).kind = 0 ∨ (profileClassWitness p).kind = 1 ∨
        (profileClassWitness p).kind = 2 := by omega
    rcases hk with hk | hk | hk
    · simpa only [profileZeroClass, hp, ite_false, hk, Fin.val_zero, zero_add,
        profileEvenClassSum, Nat.one_ne_zero, ite_true] using
        profileSplit6_even_witness_sum_two _ _ (profileClassWitness p) hw hk f hm
    · simpa [profileZeroClass, hp, hk, profileEvenClassSum] using
        profileSplit6_even_witness_sum_three _ _ (profileClassWitness p) hw hk f hm
    · simpa [profileZeroClass, hp, hk, profileEvenClassSum] using
        profileSplit6_even_witness_sum_four _ _ (profileClassWitness p) hw hk f hm

theorem profileSplit6_even_total {M : Type*} [AddCommMonoid M]
    (f : (Fin 5 → Nat) → M) :
    (∑ bc : ProfileAffine4 × ProfileAffine4, ∑ q : ProfileQuad4,
      f (profileOfNibbles (profileSplit6Eval (0,bc.1,bc.2,q)))) =
      ∑ k : Fin 4, ([93,840,90,1] : List Nat).getD k.val 0 • profileEvenClassSum f k.val := by
  rw [← Equiv.sum_comp profileAffinePairEquiv]
  simp_rw [profileSplit6_even_pair_sum]
  exact profileZeroClass_sum (profileEvenClassSum f)

theorem profileSplit6_odd_total {M : Type*} [AddCommMonoid M]
    (f : (Fin 5 → Nat) → M) :
    (∑ bc : ProfileAffine4 × ProfileAffine4, ∑ q : ProfileQuad4,
      f (profileOfNibbles (profileSplit6Eval (1,bc.1,bc.2,q)))) =
      1024 • ∑ w : Fin 17, profileQuad4WeightCounts.getD w.val 0 • f (profileOddShape w.val) := by
  simp_rw [profileSplit6_odd_sum]
  simp only [Finset.sum_const, Finset.card_univ, Fintype.card_prod,
    Fintype.card_fun, Fintype.card_fin, ZMod.card]
  rfl

theorem profileSplit6_total {M : Type*} [AddCommMonoid M]
    (f : (Fin 5 → Nat) → M) :
    (∑ t : ProfileSplit6, f (profileOfNibbles (profileSplit6Eval t))) =
      (∑ k : Fin 4, ([93,840,90,1] : List Nat).getD k.val 0 • profileEvenClassSum f k.val) +
      1024 • ∑ w : Fin 17, profileQuad4WeightCounts.getD w.val 0 • f (profileOddShape w.val) := by
  have htwo (g : ZMod 2 → M) : (∑ a : ZMod 2, g a) = g 0 + g 1 := by
    change (∑ a : Fin 2, g a) = _
    exact Fin.sum_univ_two _
  simp only [Fintype.sum_prod_type]
  rw [htwo]
  have he := profileSplit6_even_total f
  have ho := profileSplit6_odd_total f
  simp only [Fintype.sum_prod_type] at he ho
  rw [he, ho]

/-- The complete profile observable identity for the actual sampled 22-bit
production state. All finite tables enter through kernel-checked histograms
or checked basis-preserver/affine-pair certificates. -/
theorem productionPacketProfile_census {M : Type*} [AddCommMonoid M]
    (f : (Fin 5 → Nat) → M) :
    (∑ q : State22, f (packetProfile64 (fun j => coordMap j q))) =
      (∑ k : Fin 4, ([93,840,90,1] : List Nat).getD k.val 0 • profileEvenClassSum f k.val) +
      1024 • ∑ w : Fin 17, profileQuad4WeightCounts.getD w.val 0 • f (profileOddShape w.val) := by
  simp_rw [productionProfileSplit_profile]
  calc
    _ = ∑ t : ProfileSplit6, f (profileOfNibbles (profileSplit6Eval t)) :=
      Equiv.sum_comp productionProfileSplitEquiv (fun t => f (profileOfNibbles (profileSplit6Eval t)))
    _ = _ := profileSplit6_total f

#print axioms productionPacketProfile_census

end Spin.UniformMixer
