import SpinCodes.UniformMixer.Map22PackedBridge

/-!
# Four-variable quadratics and packet profiles

The small basis is used to authenticate the packet-profile census structurally.
No enumeration of the 22-bit production-state space is introduced here.
-/
noncomputable section
namespace Spin.UniformMixer

abbrev ProfileQuad4 := Fin 11 → ZMod 2
abbrev ProfileAffine4 := Fin 5 → ZMod 2
abbrev ProfileSplit6 := ZMod 2 × ProfileAffine4 × ProfileAffine4 × ProfileQuad4

def profileQuad4Rows : List Nat :=
  [0xffff, 0xaaaa, 0xcccc, 0xf0f0, 0xff00,
   0x8888, 0xa0a0, 0xaa00, 0xc0c0, 0xcc00, 0xf000]

def profileQuad4Monomial : Fin 11 → Nat := ![0,1,2,4,8,3,5,9,6,10,12]

def profileQuad4Eval (q : ProfileQuad4) (v : Fin 16) : ZMod 2 :=
  ∑ i : Fin 11, bitScalar ((profileQuad4Rows.getD i.val 0).testBit v.val) * q i

def profileAffine4Eval (a : ProfileAffine4) (v : Fin 16) : ZMod 2 :=
  ∑ i : Fin 5, bitScalar ((profileQuad4Rows.getD i.val 0).testBit v.val) * a i

def profileQuad4Recover (f : Fin 16 → ZMod 2) (i : Fin 11) : ZMod 2 :=
  ∑ v : Fin 16, (if Nat.land v.val (profileQuad4Monomial i) = v.val then 1 else 0) * f v

theorem profileQuad4_basis_inverse : ∀ i j : Fin 11,
    (∑ v : Fin 16,
      (if Nat.land v.val (profileQuad4Monomial i) = v.val then (1 : ZMod 2) else 0) *
        bitScalar ((profileQuad4Rows.getD j.val 0).testBit v.val)) =
      if j = i then 1 else 0 := by decide

theorem profileQuad4Recover_eval (q : ProfileQuad4) :
    profileQuad4Recover (profileQuad4Eval q) = q := by
  funext i
  unfold profileQuad4Recover profileQuad4Eval
  simp only [Finset.mul_sum]
  rw [Finset.sum_comm]
  have hs (j : Fin 11) :
      (∑ v : Fin 16,
        (if Nat.land v.val (profileQuad4Monomial i) = v.val then (1 : ZMod 2) else 0) *
          (bitScalar ((profileQuad4Rows.getD j.val 0).testBit v.val) * q j)) =
      (∑ v : Fin 16,
        (if Nat.land v.val (profileQuad4Monomial i) = v.val then (1 : ZMod 2) else 0) *
          bitScalar ((profileQuad4Rows.getD j.val 0).testBit v.val)) * q j := by
    rw [Finset.sum_mul]
    apply Finset.sum_congr rfl
    intro v _; ring
  simp_rw [hs, profileQuad4_basis_inverse]
  simp

theorem profileQuad4Eval_injective : Function.Injective profileQuad4Eval := by
  intro q r h
  have hh := congrArg profileQuad4Recover h
  simpa only [profileQuad4Recover_eval] using hh

def profileNibbleEval (a b c d : ZMod 2) (u : Fin 4) : ZMod 2 :=
  d + b * bitScalar (u.val.testBit 0) + c * bitScalar (u.val.testBit 1) +
    a * bitScalar (u.val.testBit 0) * bitScalar (u.val.testBit 1)

def profileNibbleWeight (a b c d : ZMod 2) : Nat :=
  ∑ u : Fin 4, if profileNibbleEval a b c d u = 1 then 1 else 0

/-- The low-two-variable truth table is the complete packet-weight reduction. -/
theorem profileNibbleWeight_table : ∀ a b c d : ZMod 2,
    profileNibbleWeight a b c d =
      if a = 0 then (if b = 0 ∧ c = 0 then (if d = 0 then 0 else 4) else 2)
      else (if d + b*c = 0 then 1 else 3) := by decide

def profileSplit6Eval (t : ProfileSplit6) (vu : Fin 16 × Fin 4) : ZMod 2 :=
  profileNibbleEval t.1 (profileAffine4Eval t.2.1 vu.1)
    (profileAffine4Eval t.2.2.1 vu.1) (profileQuad4Eval t.2.2.2 vu.1) vu.2

def profileNibbleWordWeight (f : Fin 16 × Fin 4 → ZMod 2) (v : Fin 16) : Nat :=
  ∑ u : Fin 4, if f (v,u) = 1 then 1 else 0

def profileOfNibbles (f : Fin 16 × Fin 4 → ZMod 2) (r : Fin 5) : Nat :=
  (Finset.univ.filter (fun v => profileNibbleWordWeight f v = r.val)).card

def packetProfile64 (f : Fin 64 → ZMod 2) : Fin 5 → Nat :=
  profileOfNibbles (fun vu => f (finProdFinEquiv vu))

def profileQuad4Weight (q : ProfileQuad4) : Nat :=
  ∑ v : Fin 16, if profileQuad4Eval q v = 1 then 1 else 0

#print axioms profileQuad4Recover_eval
#print axioms profileNibbleWeight_table
end Spin.UniformMixer
