import SpinCodes.Turbo.PhysicalProfileCompleteRows22
import SpinCodes.Turbo.PhysicalProfileBirthFamilyK20Tilt222

/-!
# Authenticated local transfer family at the retained rational parameter

Every occupation and every row is connected to the actual physical kernel,
including the exact small-input birth calculations. The selected-step
family includes the proven binomial thinning. This is one local family,
not an end-to-end distance certificate or coverage of all required tilts.
-/
noncomputable section
namespace Spin.Turbo.ProfileFamilyK20Tilt222

def physicalRows : Fin 17 → Matrix (Fin 11) (Fin 11) ℝ :=
  profileFamilyRows22 birthSmokeFamily22 firstCaps returnCaps endpointDen

def selectedRows (j : Fin 17) : Matrix (Fin 11) (Fin 11) ℝ :=
  turboThinnedRow22 j birthSmokeZ22 physicalRows

theorem physicalRows_sound (j : Fin 17) :
    physicalRowSound22 (Matrix.of (turboPhysicalKernel22 j birthSmokeZ22)) (physicalRows j) :=
  profileFamilyRows22_sound birthSmokeZ22_mem birthSmokeZ22_nonneg
    firstCaps secondCaps returnCaps feedbackActiveCaps zeroHoldCaps endpointDen familyChecked
    birthSmokeFamily22 birthSmokeFamily22_sound j

theorem selectedRows_sound (j : Fin 17) :
    physicalRowSound22 (Matrix.of (turboSelectedKernel22 j birthSmokeZ22)) (selectedRows j) :=
  profileFamilySelectedRows22_sound birthSmokeZ22_mem birthSmokeZ22_nonneg
    firstCaps secondCaps returnCaps feedbackActiveCaps zeroHoldCaps endpointDen familyChecked
    birthSmokeFamily22 birthSmokeFamily22_sound j

theorem physicalRows_nonneg (j : Fin 17) (i k : Fin 11) : 0 ≤ physicalRows j i k :=
  profileFamilyRows22_nonneg firstCaps secondCaps returnCaps feedbackActiveCaps zeroHoldCaps
    endpointDen familyChecked birthSmokeFamily22 birthSmokeFamily22_sound j i k

theorem selectedRows_nonneg (j : Fin 17) (i k : Fin 11) : 0 ≤ selectedRows j i k :=
  profileFamilySelectedRows22_nonneg firstCaps secondCaps returnCaps feedbackActiveCaps zeroHoldCaps
    endpointDen familyChecked birthSmokeFamily22 birthSmokeFamily22_sound j i k

#print axioms physicalRows_sound
#print axioms selectedRows_sound
#print axioms physicalRows_nonneg
#print axioms selectedRows_nonneg
end Spin.Turbo.ProfileFamilyK20Tilt222
