import SpinCodes.UniformMixer.Map22PacketProfileClassifier
import SpinCodes.UniformMixer.Map22PacketProfileAffine
import SpinCodes.UniformMixer.Map22PacketProfileHistogram

noncomputable section
namespace Spin.UniformMixer
set_option maxHeartbeats 1000000
set_option maxRecDepth 10000

theorem profileZeroClass_ofFn : List.ofFn profileZeroClass = profileZeroKinds := by
  rw [List.ofFn_mul (m := 32) (n := 32)]
  have hb : ∀ b : Fin 32,
      List.ofFn (fun r : Fin 32 => profileZeroClass ⟨b.val*32+r.val, by omega⟩) =
        profileZeroKindChunks.getD b.val [] := by
    intro b
    fin_cases b
    · exact profileZeroKinds_block0
    · exact profileZeroKinds_block1
    · exact profileZeroKinds_block2
    · exact profileZeroKinds_block3
    · exact profileZeroKinds_block4
    · exact profileZeroKinds_block5
    · exact profileZeroKinds_block6
    · exact profileZeroKinds_block7
    · exact profileZeroKinds_block8
    · exact profileZeroKinds_block9
    · exact profileZeroKinds_block10
    · exact profileZeroKinds_block11
    · exact profileZeroKinds_block12
    · exact profileZeroKinds_block13
    · exact profileZeroKinds_block14
    · exact profileZeroKinds_block15
    · exact profileZeroKinds_block16
    · exact profileZeroKinds_block17
    · exact profileZeroKinds_block18
    · exact profileZeroKinds_block19
    · exact profileZeroKinds_block20
    · exact profileZeroKinds_block21
    · exact profileZeroKinds_block22
    · exact profileZeroKinds_block23
    · exact profileZeroKinds_block24
    · exact profileZeroKinds_block25
    · exact profileZeroKinds_block26
    · exact profileZeroKinds_block27
    · exact profileZeroKinds_block28
    · exact profileZeroKinds_block29
    · exact profileZeroKinds_block30
    · exact profileZeroKinds_block31
  simp_rw [hb]
  exact (show (List.ofFn (fun b : Fin 32 => profileZeroKindChunks.getD b.val [])).flatten =
    profileZeroKinds by decide)

theorem profile_count_ofFn {n : Nat} (f : Fin n → Nat) (k : Nat) :
    (List.ofFn f).count k = ∑ i : Fin n, if f i = k then 1 else 0 := by
  induction n with
  | zero => simp
  | succ n ih =>
    rw [List.ofFn_succ, List.count_cons, Fin.sum_univ_succ, ih]
    by_cases h : f 0 = k
    · simp only [h, beq_self_eq_true, ite_true]
      omega
    · have hh : (f 0 == k) = false := beq_eq_false_iff_ne.mpr h
      simp only [h, hh, Bool.false_eq_true, ite_false, zero_add, add_zero]

theorem profileZeroClass_count (k : Fin 4) :
    (∑ p : Fin 1024, if profileZeroClass p = k.val then 1 else 0) =
      ([93,840,90,1] : List Nat).getD k.val 0 := by
  rw [← profile_count_ofFn, profileZeroClass_ofFn]
  exact profileZeroKinds_counts k

theorem profileZeroClass_bound (p : Fin 1024) : profileZeroClass p ≤ 3 := by
  unfold profileZeroClass
  split_ifs
  · omega
  · have := (profileClassWitness p).kind.isLt
    omega

theorem profileZeroClass_sum {M : Type*} [AddCommMonoid M] (f : Nat → M) :
    (∑ p : Fin 1024, f (profileZeroClass p)) =
      ∑ k : Fin 4, ([93,840,90,1] : List Nat).getD k.val 0 • f k.val := by
  rw [profile_sum_by_weight 3 profileZeroClass profileZeroClass_bound]
  simp_rw [profileZeroClass_count]

theorem profileAffinePair_zero_iff (p : Fin 1024) (v : Fin 16) :
    (profileAffine4Eval (profileAffinePairEquiv p).1 v = 0 ∧
      profileAffine4Eval (profileAffinePairEquiv p).2 v = 0) ↔
      (!(profilePackedAffine (p.val%32)).testBit v.val &&
        !(profilePackedAffine (p.val/32)).testBit v.val) = true := by
  rw [profileAffinePair_fst, profileAffinePair_snd, profileAffinePacked_eval,
    profileAffinePacked_eval]
  simp only [profilePackedAffine, bitScalar_eq_zero]
  cases (Spin.Structured.PackedMap.eval (profileQuad4Rows.take 5) (p.val%32)).testBit v.val <;>
    cases (Spin.Structured.PackedMap.eval (profileQuad4Rows.take 5) (p.val/32)).testBit v.val <;> decide

theorem profileClassWitness_checked (p : Fin 1024) (hp : profileZeroIndex p ≠ 0) :
    (profileClassWitness p).check = true := by
  have hh := profileZeroClassifier_checked p
  simp only [profileZeroClassifierValid, hp, ite_false] at hh
  apply profileZeroBank_checked
  unfold profileClassWitness
  rw [List.getD_eq_getElem _ _ hh.1]
  exact List.getElem_mem _

theorem profileClassWitness_alignment (p : Fin 1024) (hp : profileZeroIndex p ≠ 0)
    (v : Fin 16) :
    (profileClassWitness p).zeroMask.testBit v.val = true ↔
      profileAffine4Eval (profileAffinePairEquiv p).1 v = 0 ∧
      profileAffine4Eval (profileAffinePairEquiv p).2 v = 0 := by
  have hh := profileZeroClassifier_checked p
  simp only [profileZeroClassifierValid, hp, ite_false] at hh
  rw [hh.2 v, profileAffinePair_zero_iff]
  cases (profilePackedAffine (p.val%32)).testBit v.val <;>
    cases (profilePackedAffine (p.val/32)).testBit v.val <;> decide

theorem profileClassZero_empty (p : Fin 1024) (hp : profileZeroIndex p = 0)
    (v : Fin 16) :
    ¬ (profileAffine4Eval (profileAffinePairEquiv p).1 v = 0 ∧
      profileAffine4Eval (profileAffinePairEquiv p).2 v = 0) := by
  have hh := profileZeroClassifier_checked p
  simp only [profileZeroClassifierValid, hp, ite_true] at hh
  rw [profileAffinePair_zero_iff]
  have h := hh v
  cases ha : (profilePackedAffine (p.val%32)).testBit v.val <;>
    cases hb : (profilePackedAffine (p.val/32)).testBit v.val <;> simp_all

#print axioms profileZeroClass_sum
end Spin.UniformMixer
