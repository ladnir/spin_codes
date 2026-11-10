import SpinCodes.Turbo.Hi.PhysicalProfileRows22
import SpinCodes.Turbo.Hi.PhysicalProfileFamily22
import SpinCodes.Turbo.Hi.PhysicalProfileFamilyExtract22
import SpinCodes.Turbo.Hi.PhysicalProfileFamilyFourier22

/-! The checked first-moment and Fourier bounds supply coarse birth cones.
Exact small-occupancy cones can then replace precisely rows zero, one and two. -/
noncomputable section
namespace Spin.Turbo.Hi
open Spin.NumericHi Spin.NumericHi.Fix

def profileFamilyBirth22 (small : Fin 3 → Fin 11 → ℝ)
    (first : List (List Int)) (zeroHold : List Int) (den : Int)
    (j : Fin 17) : Fin 11 → ℝ :=
  if h : j.val < 3 then small ⟨j.val, h⟩ else
    physicalCoarseBirthRow22 ((zeroHold.getD j.val 0 : ℝ)/(den : ℝ))
      (((first.getD j.val []).getD 0 0 : ℝ)/(den : ℝ))

theorem profileFamilyBirth22_coarse {box : Fix} {z : ℝ}
    (hzbox : Mem box z) (hz : 0 ≤ z)
    (first second returns : List (List Int)) (active zeroHold : List Int) (den : Int)
    (h : profileFamilyCheck22 profileNumericEntries22 box first second returns active zeroHold den = true)
    (j : Fin 17) :
    physicalCone22 (physicalCoarseBirthRow22
      ((zeroHold.getD j.val 0 : ℝ)/(den : ℝ))
      (((first.getD j.val []).getD 0 0 : ℝ)/(den : ℝ)))
      (turboExactBirthMass22 j z) := by
  have hj := profileFamilyCheck22_row _ _ _ _ _ _ _ _ h j
  obtain ⟨hs, hf, _, _, hzero, _⟩ :=
    profileFamilyRowCheck22_parts _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ hj
  apply turboExactBirthMass22_coarse_cone j z hz
  · exact physicalZeroHold22_cached_upper hzbox j _ _ hs.1 hzero
  · rw [physicalInputMoment22_eq_profile]
    exact profileCachedCapsCheck22_actual _ hs.1 hzbox j hf 0 (by decide)
      0 (by simp [physicalClassAllows22])

theorem profileFamilyBirth22_sound {box : Fix} {z : ℝ}
    (hzbox : Mem box z) (hz : 0 ≤ z)
    (first second returns : List (List Int)) (active zeroHold : List Int) (den : Int)
    (h : profileFamilyCheck22 profileNumericEntries22 box first second returns active zeroHold den = true)
    (small : Fin 3 → Fin 11 → ℝ)
    (hsmall : ∀ j : Fin 3, physicalCone22 (small j)
      (turboExactBirthMass22 ⟨j.val, by omega⟩ z)) (j : Fin 17) :
    physicalCone22 (profileFamilyBirth22 small first zeroHold den j)
      (turboExactBirthMass22 j z) := by
  unfold profileFamilyBirth22
  split_ifs with hj
  · exact hsmall ⟨j.val, hj⟩
  · exact profileFamilyBirth22_coarse hzbox hz first second returns active zeroHold den h j

#print axioms profileFamilyBirth22_coarse
#print axioms profileFamilyBirth22_sound
end Spin.Turbo.Hi
