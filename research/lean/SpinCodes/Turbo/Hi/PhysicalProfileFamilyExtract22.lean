import SpinCodes.Turbo.Hi.PhysicalProfileFamilyDefs22

namespace Spin.Turbo.Hi
open Spin.NumericHi

theorem profileFamilyCheck22_row (entries : List (Nat × List Nat)) (box : Fix)
    (first second returns : List (List Int)) (active zeroHold : List Int) (den : Int)
    (h : profileFamilyCheck22 entries box first second returns active zeroHold den = true)
    (j : Fin 17) :
    profileFamilyRowCheck22 entries
      (profileFixCache22 entries (profileFixMomentBank22 box))
      (profileFixCache22 entries (profileFixMomentBank22 (Fix.mul box box)))
      (profileFixCache22 entries (profileFixFourierBank22 box))
      (profileFixCache22 entries (profileFixFourierBank22 (Fix.ofInt 1)))
      (profileFixMomentBank22 box) (profileFixMomentBank22 (Fix.mul box box))
      (profileFixFourierBank22 box) (profileFixFourierBank22 (Fix.ofInt 1))
      j.val (first.getD j.val []) (second.getD j.val []) (returns.getD j.val [])
      (active.getD j.val 0) (zeroHold.getD j.val 0) den = true := by
  simp only [profileFamilyCheck22, Bool.and_eq_true] at h
  exact List.all_eq_true.mp h.2 j.val (List.mem_range.mpr j.isLt)

theorem profileFamilyRowCheck22_parts (entries : List (Nat × List Nat))
    (fc sc zc bc : ProfileFixCache22) (fb sb zb bb : List (Fix × Fix))
    (j : Nat) (first second returns : List Int) (active zeroHold den : Int)
    (h : profileFamilyRowCheck22 entries fc sc zc bc fb sb zb bb j
      first second returns active zeroHold den = true) :
    (0 < den ∧ 0 ≤ active ∧ active ≤ den ∧ 0 ≤ zeroHold) ∧
    profileCachedCapsCheck22 entries fc fb j first den = true ∧
    profileCachedCapsCheck22 entries sc sb j second den = true ∧
    profileFixLowerCheck22 (profileFixCachedAllFourier22 entries bc bb j) (den-active) den = true ∧
    profileFixUpperCheck22 (profileFixCachedAllFourier22 entries zc zb j) zeroHold den = true ∧
    profileReturnCapsCheck22 first second returns active = true := by
  simpa only [profileFamilyRowCheck22, Bool.and_eq_true, decide_eq_true_eq, and_assoc] using h

#print axioms profileFamilyCheck22_row
#print axioms profileFamilyRowCheck22_parts
end Spin.Turbo.Hi

