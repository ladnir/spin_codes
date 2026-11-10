import SpinCodes.UniformMixer.Map22PacketProfileSmallDefs

/-! Symbolic affine-product closure and the two packet-profile cases. -/
noncomputable section
namespace Spin.UniformMixer

def profileAffineProduct (b c : ProfileAffine4) : ProfileQuad4 := ![
  b 0*c 0,
  b 0*c 1+b 1*c 0+b 1*c 1,
  b 0*c 2+b 2*c 0+b 2*c 2,
  b 0*c 3+b 3*c 0+b 3*c 3,
  b 0*c 4+b 4*c 0+b 4*c 4,
  b 1*c 2+b 2*c 1,
  b 1*c 3+b 3*c 1,
  b 1*c 4+b 4*c 1,
  b 2*c 3+b 3*c 2,
  b 2*c 4+b 4*c 2,
  b 3*c 4+b 4*c 3]

theorem profileQuad4Eval_add (q r : ProfileQuad4) (v : Fin 16) :
    profileQuad4Eval (q+r) v = profileQuad4Eval q v + profileQuad4Eval r v := by
  simp [profileQuad4Eval, Pi.add_apply, mul_add, Finset.sum_add_distrib]

theorem profileAffineProduct_eval (b c : ProfileAffine4) (v : Fin 16) :
    profileQuad4Eval (profileAffineProduct b c) v = profileAffine4Eval b v * profileAffine4Eval c v := by
  fin_cases v <;>
    simp [profileQuad4Eval, profileAffine4Eval, profileQuad4Rows, profileAffineProduct,
      Fin.sum_univ_succ, bitScalar, Nat.testBit_succ] <;> ring

def profileQuad4Shift (p : ProfileQuad4) : ProfileQuad4 ≃ ProfileQuad4 where
  toFun q := q+p
  invFun q := q-p
  left_inv q := by simp
  right_inv q := by simp

/-- With low quadratic coefficient one, all packets have weight one or three.
The governing high polynomial is one translated four-variable quadratic. -/
theorem profileOddNibbleWeight (b c : ProfileAffine4) (q : ProfileQuad4) (v : Fin 16) :
    profileNibbleWordWeight (profileSplit6Eval (1,b,c,q)) v =
      if profileQuad4Eval (q+profileAffineProduct b c) v = 0 then 1 else 3 := by
  change profileNibbleWeight 1 (profileAffine4Eval b v) (profileAffine4Eval c v) (profileQuad4Eval q v) = _
  rw [profileNibbleWeight_table]
  simp only [one_ne_zero, ite_false]
  rw [profileQuad4Eval_add, profileAffineProduct_eval]

/-- With low quadratic coefficient zero, the only non-balanced packets lie
in the common zero set of the two affine high-variable coefficients. -/
theorem profileEvenNibbleWeight (b c : ProfileAffine4) (q : ProfileQuad4) (v : Fin 16) :
    profileNibbleWordWeight (profileSplit6Eval (0,b,c,q)) v =
      if profileAffine4Eval b v = 0 ∧ profileAffine4Eval c v = 0 then
        (if profileQuad4Eval q v = 0 then 0 else 4) else 2 := by
  change profileNibbleWeight 0 (profileAffine4Eval b v) (profileAffine4Eval c v) (profileQuad4Eval q v) = _
  rw [profileNibbleWeight_table]
  simp only [ite_true]

#print axioms profileAffineProduct_eval
#print axioms profileOddNibbleWeight
#print axioms profileEvenNibbleWeight
end Spin.UniformMixer
