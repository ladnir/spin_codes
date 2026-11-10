import SpinCodes.Turbo.PhysicalProfileNumeric22

noncomputable section
namespace Spin.Turbo
open Spin.Numeric Spin.Numeric.Fix Polynomial

def profileRealFactors22 (profile : List Nat) (bank : List (ℝ × ℝ)) : List (ℝ × ℝ) :=
  ((List.range 5).map (fun r => List.replicate (profile.getD r 0) (bank.getD r (0,0)))).flatten

def profileRealMomentBank22 (z : ℝ) : List (ℝ × ℝ) :=
  (List.range 5).map (fun r => (z^r,(1+z)^4-z^r))

def profileRealFourierBank22 (z : ℝ) : List (ℝ × ℝ) :=
  (List.range 5).map (fun r => (1,(1+z)^(4-r)*(1-z)^r-1))

def profileBankPolynomial22 (profile : List Nat) (bank : List (ℝ × ℝ)) : Polynomial ℝ :=
  ∏ r : Fin 5, (C (bank.getD r.val (0,0)).1 + X*C (bank.getD r.val (0,0)).2) ^ profile.getD r.val 0

theorem profileForall₂_getD22 {A B : Type*} {R : A → B → Prop}
    {as : List A} {bs : List B} (h : List.Forall₂ R as bs)
    (da : A) (db : B) (hd : R da db) (j : Nat) : R (as.getD j da) (bs.getD j db) := by
  induction h generalizing j with
  | nil => exact hd
  | cons hab htail ih => cases j with | zero => exact hab | succ j => exact ih j

theorem profileForall₂_replicate22 {A B : Type*} {R : A → B → Prop}
    {a : A} {b : B} (h : R a b) (n : Nat) :
    List.Forall₂ R (List.replicate n a) (List.replicate n b) := by
  induction n with
  | zero => exact .nil
  | succ n ih => exact .cons h ih

theorem profileFixFactors22_mem (profile : List Nat) {fs : List (Fix × Fix)} {rs : List (ℝ × ℝ)}
    (h : List.Forall₂ (fun f r => Mem f.1 r.1 ∧ Mem f.2 r.2) fs rs) :
    List.Forall₂ (fun f r => Mem f.1 r.1 ∧ Mem f.2 r.2)
      (profileFixFactors22 profile fs) (profileRealFactors22 profile rs) := by
  unfold profileFixFactors22 profileRealFactors22
  apply List.rel_flatten
  rw [List.forall₂_map_left_iff, List.forall₂_map_right_iff, List.forall₂_same]
  intro r _
  apply profileForall₂_replicate22
  apply profileForall₂_getD22 h (Fix.zero,Fix.zero) (0,0)
  constructor <;> simp [Mem,Fix.zero]

theorem profileRealFactors22_polynomial (profile : List Nat) (bank : List (ℝ × ℝ)) :
    ((profileRealFactors22 profile bank).map (fun p : ℝ × ℝ => C p.1 + X*C p.2)).prod =
      profileBankPolynomial22 profile bank := by
  have hg : ∀ ks : List Nat,
      (((ks.map (fun r => List.replicate (profile.getD r 0) (bank.getD r (0,0)))).flatten).map
        (fun p : ℝ × ℝ => C p.1 + X*C p.2)).prod =
        (ks.map (fun r => (C (bank.getD r (0,0)).1 + X*C (bank.getD r (0,0)).2)^profile.getD r 0)).prod := by
    intro ks
    induction ks with
    | nil => simp
    | cons k ks ih =>
      simp only [List.map_cons, List.flatten_cons, List.map_append, List.prod_append,
        List.map_replicate, List.prod_replicate, List.prod_cons, ih]
  rw [profileRealFactors22, hg]
  simp [profileBankPolynomial22, List.range_succ, Fin.prod_univ_succ]

theorem profileFixBankCoefficient22_sound (profile : List Nat)
    {fs : List (Fix × Fix)} {rs : List (ℝ × ℝ)}
    (h : List.Forall₂ (fun f r => Mem f.1 r.1 ∧ Mem f.2 r.2) fs rs) (j : Nat) :
    Mem ((profileFixCoefficients22 (profileFixFactors22 profile fs)).getD j Fix.zero)
      ((profileBankPolynomial22 profile rs).coeff j) := by
  rw [← profileRealFactors22_polynomial]
  exact profileFixCoefficients22_sound (profileFixFactors22_mem profile h) j

theorem profileFixMomentBank22_mem {box : Fix} {z : ℝ} (hz : Mem box z) :
    List.Forall₂ (fun f r => Mem f.1 r.1 ∧ Mem f.2 r.2)
      (profileFixMomentBank22 box) (profileRealMomentBank22 z) := by
  have hone : Mem (Fix.ofInt 1) (1 : ℝ) := by simpa using ofInt_mem 1
  have hz2 : Mem (Fix.mul box box) (z^2) := by simpa only [pow_two] using mul_mem hz hz
  have hz3 : Mem (Fix.mul (Fix.mul box box) box) (z^3) := by
    convert mul_mem hz2 hz using 1 <;> ring
  have hz4 : Mem (Fix.mul (Fix.mul box box) (Fix.mul box box)) (z^4) := by
    convert mul_mem hz2 hz2 using 1 <;> ring
  have ht : Mem (Fix.pow (Fix.add (Fix.ofInt 1) box) 4) ((1+z)^4) :=
    pow_mem (add_mem hone hz) 4
  change List.Forall₂ _ _ [(z^0,(1+z)^4-z^0), (z^1,(1+z)^4-z^1),
    (z^2,(1+z)^4-z^2), (z^3,(1+z)^4-z^3), (z^4,(1+z)^4-z^4)]
  exact .cons ⟨by simpa only [pow_zero] using hone, by simpa only [pow_zero] using sub_mem ht hone⟩
    (.cons ⟨by simpa only [pow_one] using hz, by simpa only [pow_one] using sub_mem ht hz⟩
      (.cons ⟨hz2,sub_mem ht hz2⟩ (.cons ⟨hz3,sub_mem ht hz3⟩ (.cons ⟨hz4,sub_mem ht hz4⟩ .nil))))

theorem profileFixFourierBank22_mem {box : Fix} {z : ℝ} (hz : Mem box z) :
    List.Forall₂ (fun f r => Mem f.1 r.1 ∧ Mem f.2 r.2)
      (profileFixFourierBank22 box) (profileRealFourierBank22 z) := by
  have hone : Mem (Fix.ofInt 1) (1 : ℝ) := by simpa using ofInt_mem 1
  unfold profileFixFourierBank22 profileRealFourierBank22
  rw [List.forall₂_map_left_iff, List.forall₂_map_right_iff, List.forall₂_same]
  intro r _
  exact ⟨hone, sub_mem (mul_mem (pow_mem (add_mem hone hz) (4-r))
    (pow_mem (sub_mem hone hz) r)) hone⟩

theorem profileChoose16Values22_checked : ∀ j : Fin 17,
    profileChoose16Values22.getD j.val 0 = Nat.choose 16 j.val := by decide

theorem profileLayerDenominator22_pos (j : Fin 17) : 0 < profileLayerDenominator22 j.val := by
  unfold profileLayerDenominator22
  apply mul_pos
  · have h : 0 < profileChoose16Values22.getD j.val 0 := by
      exact (show ∀ k : Fin 17, 0 < profileChoose16Values22.getD k.val 0 by decide) j
    exact_mod_cast h
  · positivity

theorem profileLayerDenominator22_real (j : Fin 17) :
    (profileLayerDenominator22 j.val : ℝ) = (Nat.choose 16 j.val : ℝ) * (15 : ℝ)^j.val := by
  unfold profileLayerDenominator22
  push_cast
  rw [profileChoose16Values22_checked]

def profileNumericMoment22 (profile : List Nat) (z : ℝ) (j : Fin 17) : ℝ :=
  (profileBankPolynomial22 profile (profileRealMomentBank22 z)).coeff j.val /
    ((Nat.choose 16 j.val : ℝ) * (15 : ℝ)^j.val)

def profileNumericFourier22 (profile : List Nat) (z : ℝ) (j : Fin 17) : ℝ :=
  (profileBankPolynomial22 profile (profileRealFourierBank22 z)).coeff j.val /
    ((Nat.choose 16 j.val : ℝ) * (15 : ℝ)^j.val)

theorem profileFixMoment22_sound (profile : List Nat) {box : Fix} {z : ℝ}
    (hz : Mem box z) (j : Fin 17) : Mem (profileFixMoment22 profile box j.val)
      (profileNumericMoment22 profile z j) := by
  have h := divInt_mem (profileFixBankCoefficient22_sound profile (profileFixMomentBank22_mem hz) j.val)
    (profileLayerDenominator22_pos j)
  simpa only [profileFixMoment22, profileNumericMoment22, profileLayerDenominator22_real] using h

theorem profileFixFourier22_sound (profile : List Nat) {box : Fix} {z : ℝ}
    (hz : Mem box z) (j : Fin 17) : Mem (profileFixFourier22 profile box j.val)
      (profileNumericFourier22 profile z j) := by
  have h := divInt_mem (profileFixBankCoefficient22_sound profile (profileFixFourierBank22_mem hz) j.val)
    (profileLayerDenominator22_pos j)
  simpa only [profileFixFourier22, profileNumericFourier22, profileLayerDenominator22_real] using h

theorem profileFixSecondMoment22_sound (profile : List Nat) {box : Fix} {z : ℝ}
    (hz : Mem box z) (j : Fin 17) : Mem (profileFixSecondMoment22 profile box j.val)
      (profileNumericMoment22 profile (z*z) j) :=
  profileFixMoment22_sound profile (mul_mem hz hz) j

#print axioms profileFixMoment22_sound
#print axioms profileFixFourier22_sound
end Spin.Turbo
