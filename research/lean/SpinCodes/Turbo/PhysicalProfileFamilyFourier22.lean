import SpinCodes.Turbo.PhysicalProfileNumericCache22
import SpinCodes.Turbo.PhysicalProfileMoment22

noncomputable section
namespace Spin.Turbo
open Spin.Numeric Spin.Numeric.Fix

theorem profileFixCachedAllFourier22_eq (entries : List (Nat × List Nat))
    (box : Fix) (j : Nat) :
    profileFixCachedAllFourier22 entries
      (profileFixCache22 entries (profileFixFourierBank22 box))
      (profileFixFourierBank22 box) j = profileFixAllFourier22 entries box j := by
  unfold profileFixCachedAllFourier22 profileFixAllFourier22
  simp_rw [profileFixCachedRaw22_eq]
  rfl

theorem physicalZeroHold22_cached_mem {box : Fix} {z : ℝ} (hz : Mem box z)
    (j : Fin 17) :
    Mem (profileFixCachedAllFourier22 profileNumericEntries22
      (profileFixCache22 profileNumericEntries22 (profileFixFourierBank22 box))
      (profileFixFourierBank22 box) j.val)
      (physicalZeroHold22 (packetStepLaw22 j) z) := by
  rw [profileFixCachedAllFourier22_eq]
  exact physicalZeroHold22_fixed_mem hz j

theorem physicalZeroHold22_cached_upper {box : Fix} {z : ℝ} (hz : Mem box z)
    (j : Fin 17) (num den : Int) (hd : 0 < den)
    (h : profileFixUpperCheck22 (profileFixCachedAllFourier22 profileNumericEntries22
      (profileFixCache22 profileNumericEntries22 (profileFixFourierBank22 box))
      (profileFixFourierBank22 box) j.val) num den = true) :
    physicalZeroHold22 (packetStepLaw22 j) z ≤ (num : ℝ)/(den : ℝ) :=
  profileFixUpperCheck22_sound (physicalZeroHold22_cached_mem hz j) hd h

theorem physicalZeroHold22_cached_lower {box : Fix} {z : ℝ} (hz : Mem box z)
    (j : Fin 17) (num den : Int) (hd : 0 < den)
    (h : profileFixLowerCheck22 (profileFixCachedAllFourier22 profileNumericEntries22
      (profileFixCache22 profileNumericEntries22 (profileFixFourierBank22 box))
      (profileFixFourierBank22 box) j.val) num den = true) :
    (num : ℝ)/(den : ℝ) ≤ physicalZeroHold22 (packetStepLaw22 j) z :=
  profileFixLowerCheck22_sound (physicalZeroHold22_cached_mem hz j) hd h

#print axioms physicalZeroHold22_cached_upper
#print axioms physicalZeroHold22_cached_lower
end Spin.Turbo
