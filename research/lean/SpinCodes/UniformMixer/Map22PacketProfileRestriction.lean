import SpinCodes.UniformMixer.Map22PacketProfileSmallDefs

/-!
# Uniform fibers of canonical quadratic restrictions

Coordinate projection onto a canonical subcube keeps seven or four of
the eleven coefficients. The other coefficients remain free. Explicit
equivalences establish the fiber multiplicities without enumerating the
production states or assuming a uniform restriction distribution.
-/
noncomputable section
namespace Spin.UniformMixer

def profileQuad4Restriction3Equiv :
    ProfileQuad4 ≃ ProfileQuad3 × (Fin 4 → ZMod 2) where
  toFun q := (profileQuad4Restrict3 q, ![q 4,q 7,q 9,q 10])
  invFun ab := ![ab.1 0,ab.1 1,ab.1 2,ab.1 3,ab.2 0,
    ab.1 4,ab.1 5,ab.2 1,ab.1 6,ab.2 2,ab.2 3]
  left_inv q := by funext i; fin_cases i <;> rfl
  right_inv ab := by
    apply Prod.ext
    · funext i; fin_cases i <;> rfl
    · funext i; fin_cases i <;> rfl

def profileQuad4Restriction2Equiv :
    ProfileQuad4 ≃ ProfileQuad2 × (Fin 7 → ZMod 2) where
  toFun q := (profileQuad4Restrict2 q, ![q 3,q 4,q 6,q 7,q 8,q 9,q 10])
  invFun ab := ![ab.1 0,ab.1 1,ab.1 2,ab.2 0,ab.2 1,
    ab.1 3,ab.2 2,ab.2 3,ab.2 4,ab.2 5,ab.2 6]
  left_inv q := by funext i; fin_cases i <;> rfl
  right_inv ab := by
    apply Prod.ext
    · funext i; fin_cases i <;> rfl
    · funext i; fin_cases i <;> rfl

theorem profileQuad4Restriction3Equiv_fst (q : ProfileQuad4) :
    (profileQuad4Restriction3Equiv q).1 = profileQuad4Restrict3 q := rfl

theorem profileQuad4Restriction2Equiv_fst (q : ProfileQuad4) :
    (profileQuad4Restriction2Equiv q).1 = profileQuad4Restrict2 q := rfl

theorem profileQuad4Restrict3_eval (q : ProfileQuad4) (v : Fin 8) :
    profileQuad4Eval q (Fin.castLE (by decide : 8 ≤ 16) v) =
      profileQuad3Eval (profileQuad4Restrict3 q) v := by
  fin_cases v <;>
    simp [profileQuad4Eval, profileQuad3Eval, profileQuad4Restrict3,
      profileQuad4Rows, profileQuad3Rows, Fin.sum_univ_succ, bitScalar, Nat.testBit_succ]

theorem profileQuad4Restrict2_eval (q : ProfileQuad4) (v : Fin 4) :
    profileQuad4Eval q (Fin.castLE (by decide : 4 ≤ 16) v) =
      profileQuad2Eval (profileQuad4Restrict2 q) v := by
  fin_cases v <;>
    simp [profileQuad4Eval, profileQuad2Eval, profileQuad4Restrict2,
      profileQuad4Rows, profileQuad2Rows, Fin.sum_univ_succ, bitScalar, Nat.testBit_succ]

theorem profileQuad4_restrict3_sum {M : Type*} [AddCommMonoid M]
    (f : ProfileQuad3 → M) :
    (∑ q : ProfileQuad4, f (profileQuad4Restrict3 q)) =
      16 • ∑ a : ProfileQuad3, f a := by
  rw [← Equiv.sum_comp profileQuad4Restriction3Equiv.symm, Fintype.sum_prod_type]
  have he (a : ProfileQuad3) (b : Fin 4 → ZMod 2) :
      profileQuad4Restrict3 (profileQuad4Restriction3Equiv.symm (a,b)) = a := by
    exact congrArg Prod.fst (profileQuad4Restriction3Equiv.apply_symm_apply (a,b))
  simp only [he, Finset.sum_const, Finset.card_univ, Fintype.card_fun,
    Fintype.card_fin, ZMod.card]
  norm_num
  rw [Finset.sum_nsmul]

theorem profileQuad4_restrict2_sum {M : Type*} [AddCommMonoid M]
    (f : ProfileQuad2 → M) :
    (∑ q : ProfileQuad4, f (profileQuad4Restrict2 q)) =
      128 • ∑ a : ProfileQuad2, f a := by
  rw [← Equiv.sum_comp profileQuad4Restriction2Equiv.symm, Fintype.sum_prod_type]
  have he (a : ProfileQuad2) (b : Fin 7 → ZMod 2) :
      profileQuad4Restrict2 (profileQuad4Restriction2Equiv.symm (a,b)) = a := by
    exact congrArg Prod.fst (profileQuad4Restriction2Equiv.apply_symm_apply (a,b))
  simp only [he, Finset.sum_const, Finset.card_univ, Fintype.card_fun,
    Fintype.card_fin, ZMod.card]
  norm_num
  rw [Finset.sum_nsmul]

theorem profileQuad4_restrict3_count (P : ProfileQuad3 → Prop) [DecidablePred P] :
    (Finset.univ.filter (fun q : ProfileQuad4 => P (profileQuad4Restrict3 q))).card =
      16 * (Finset.univ.filter P).card := by
  simpa only [Finset.card_filter, nsmul_eq_mul, Nat.cast_ofNat] using
    profileQuad4_restrict3_sum (fun a => if P a then (1 : ℕ) else 0)

theorem profileQuad4_restrict2_count (P : ProfileQuad2 → Prop) [DecidablePred P] :
    (Finset.univ.filter (fun q : ProfileQuad4 => P (profileQuad4Restrict2 q))).card =
      128 * (Finset.univ.filter P).card := by
  simpa only [Finset.card_filter, nsmul_eq_mul, Nat.cast_ofNat] using
    profileQuad4_restrict2_sum (fun a => if P a then (1 : ℕ) else 0)

#print axioms profileQuad4_restrict3_count
#print axioms profileQuad4_restrict2_count
#print axioms profileQuad4Restrict3_eval
end Spin.UniformMixer
