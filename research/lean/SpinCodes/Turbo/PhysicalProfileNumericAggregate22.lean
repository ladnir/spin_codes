import SpinCodes.Turbo.PhysicalProfileNumericFactors22
import SpinCodes.Turbo.PhysicalProfileNumericAggregateDefs22

noncomputable section
namespace Spin.Turbo
open Spin.Numeric Spin.Numeric.Fix

private theorem aggregateZeroMem22 : Mem Fix.zero 0 := by simp [Mem,Fix.zero]

theorem profileFixClassSum22_mem (entries : List (Nat × List Nat)) (i : Nat)
    (value : List Nat → Fix) (realValue : List Nat → ℝ)
    (hv : ∀ p, Mem (value p) (realValue p)) :
    Mem (profileFixClassSum22 entries i value)
      ((entries.map (fun e => if profileListClassAllows22 i e.2 then
        (e.1 : ℝ)*realValue e.2 else 0)).sum) := by
  induction entries with
  | nil => exact aggregateZeroMem22
  | cons e entries ih =>
    simp only [profileFixClassSum22, List.foldr_cons, List.map_cons, List.sum_cons]
    cases h : profileListClassAllows22 i e.2
    · simpa only [profileFixClassSum22, h, Bool.false_eq_true, ite_false, zero_add] using ih
    · simp only [h, ite_true]
      exact add_mem (mul_mem (by simpa only [Int.cast_natCast] using ofInt_mem (e.1 : Int))
        (hv e.2)) ih

theorem profileFixAllSum22_mem (entries : List (Nat × List Nat))
    (value : List Nat → Fix) (realValue : List Nat → ℝ)
    (hv : ∀ p, Mem (value p) (realValue p)) :
    Mem (profileFixAllSum22 entries value)
      ((entries.map (fun e => (e.1 : ℝ)*realValue e.2)).sum) := by
  induction entries with
  | nil => exact aggregateZeroMem22
  | cons e entries ih =>
    exact add_mem (mul_mem (by simpa only [Int.cast_natCast] using ofInt_mem (e.1 : Int))
      (hv e.2)) ih

theorem profileFixClassMean22_mem (entries : List (Nat × List Nat)) (i : Nat)
    (value : List Nat → Fix) (realValue : List Nat → ℝ)
    (hv : ∀ p, Mem (value p) (realValue p)) (hc : 0 < profileListClassCount22 entries i) :
    Mem (profileFixClassMean22 entries i value)
      (((entries.map (fun e => if profileListClassAllows22 i e.2 then
        (e.1 : ℝ)*realValue e.2 else 0)).sum) / (profileListClassCount22 entries i : ℝ)) := by
  have hh := divInt_mem (profileFixClassSum22_mem entries i value realValue hv)
    (show (0 : Int) < (profileListClassCount22 entries i : Int) by exact_mod_cast hc)
  simpa only [profileFixClassMean22, Int.cast_natCast] using hh

theorem profileFixAllMean22_mem (entries : List (Nat × List Nat))
    (value : List Nat → Fix) (realValue : List Nat → ℝ)
    (hv : ∀ p, Mem (value p) (realValue p)) :
    Mem (profileFixAllMean22 entries value)
      (((entries.map (fun e => (e.1 : ℝ)*realValue e.2)).sum) / (4194304 : ℝ)) := by
  have hh := divInt_mem (profileFixAllSum22_mem entries value realValue hv)
    (show (0 : Int) < 4194304 by decide)
  simpa only [profileFixAllMean22, Int.cast_ofNat] using hh

theorem profileFixClassUpperCheck22_sound (entries : List (Nat × List Nat)) (i : Nat)
    (value : List Nat → Fix) (realValue : List Nat → ℝ)
    (hv : ∀ p, Mem (value p) (realValue p)) {num den : Int} (hd : 0 < den)
    (h : profileFixClassUpperCheck22 entries i value num den = true)
    (e : Nat × List Nat) (he : e ∈ entries) (ha : profileListClassAllows22 i e.2 = true) :
    realValue e.2 ≤ (num : ℝ)/(den : ℝ) := by
  have hh := List.all_eq_true.mp h e he
  simp only [ha, ite_true] at hh
  exact profileFixUpperCheck22_sound (hv e.2) hd hh

theorem profileFixRawMoment22_mem (profile : List Nat) {box : Fix} {z : ℝ}
    (hz : Mem box z) (j : Nat) : Mem (profileFixRawMoment22 profile box j)
      ((profileBankPolynomial22 profile (profileRealMomentBank22 z)).coeff j) :=
  profileFixBankCoefficient22_sound profile (profileFixMomentBank22_mem hz) j

theorem profileFixRawFourier22_mem (profile : List Nat) {box : Fix} {z : ℝ}
    (hz : Mem box z) (j : Nat) : Mem (profileFixRawFourier22 profile box j)
      ((profileBankPolynomial22 profile (profileRealFourierBank22 z)).coeff j) :=
  profileFixBankCoefficient22_sound profile (profileFixFourierBank22_mem hz) j

theorem profileFixClassMomentMean22_mem (entries : List (Nat × List Nat)) (i : Nat)
    {box : Fix} {z : ℝ} (hz : Mem box z) (j : Fin 17)
    (hc : 0 < profileListClassCount22 entries i) :
    Mem (profileFixClassMomentMean22 entries i box j.val)
      (((entries.map (fun e => if profileListClassAllows22 i e.2 then (e.1 : ℝ) *
        (profileBankPolynomial22 e.2 (profileRealMomentBank22 z)).coeff j.val else 0)).sum) /
          ((profileListClassCount22 entries i : ℝ) * ((Nat.choose 16 j.val : ℝ)*15^j.val))) := by
  have hh := divInt_mem (profileFixClassSum22_mem entries i _ _
    (fun p => profileFixRawMoment22_mem p hz j.val))
    (mul_pos (show (0 : Int) < (profileListClassCount22 entries i : Int) by exact_mod_cast hc)
      (profileLayerDenominator22_pos j))
  simpa only [profileFixClassMomentMean22, Int.cast_mul, Int.cast_natCast,
    profileLayerDenominator22_real] using hh

theorem profileFixAllFourier22_mem (entries : List (Nat × List Nat))
    {box : Fix} {z : ℝ} (hz : Mem box z) (j : Fin 17) :
    Mem (profileFixAllFourier22 entries box j.val)
      (((entries.map (fun e => (e.1 : ℝ) *
        (profileBankPolynomial22 e.2 (profileRealFourierBank22 z)).coeff j.val)).sum) /
          ((4194304 : ℝ) * ((Nat.choose 16 j.val : ℝ)*15^j.val))) := by
  have hh := divInt_mem (profileFixAllSum22_mem entries _ _
    (fun p => profileFixRawFourier22_mem p hz j.val))
    (mul_pos (show (0 : Int) < 4194304 by decide) (profileLayerDenominator22_pos j))
  simpa only [profileFixAllFourier22, Int.cast_mul, Int.cast_ofNat,
    profileLayerDenominator22_real] using hh

#print axioms profileFixClassMean22_mem
#print axioms profileFixClassUpperCheck22_sound
end Spin.Turbo
