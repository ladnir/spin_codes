import SpinCodes.Turbo.PhysicalProfileBirthFamily22
import SpinCodes.Turbo.PhysicalSmallBirth22K20Tilt2BudgetSmoke
import SpinCodes.Turbo.PhysicalProfileFamilyK20Tilt222Check

noncomputable section
namespace Spin.Turbo.ProfileFamilyK20Tilt222

def smallBirthSmokeCaps22 (j : Fin 3) (i : Fin 11) : ℝ :=
  if j = 0 then (birthZeroCaps22.getD i.val 0 : ℝ)/(endpointDen : ℝ)
  else if j = 1 then (birthOneCaps.getD i.val 0 : ℝ)/(endpointDen : ℝ)
  else (birthTwoCaps.getD i.val 0 : ℝ)/(endpointDen : ℝ)

theorem smallBirthSmokeCaps22_sound (j : Fin 3) :
    physicalCone22 (smallBirthSmokeCaps22 j)
      (turboExactBirthMass22 ⟨j.val, by omega⟩ birthSmokeZ22) := by
  fin_cases j
  · change physicalCone22 (fun i => (birthZeroCaps22.getD i.val 0 : ℝ)/(endpointDen : ℝ))
      (turboExactBirthMass22 0 birthSmokeZ22)
    exact smallBirthSmokeZero22
  · change physicalCone22 (fun i => (birthOneCaps.getD i.val 0 : ℝ)/(endpointDen : ℝ))
      (turboExactBirthMass22 1 birthSmokeZ22)
    exact smallBirthSmokeOne22
  · change physicalCone22 (fun i => (birthTwoCaps.getD i.val 0 : ℝ)/(endpointDen : ℝ))
      (turboExactBirthMass22 2 birthSmokeZ22)
    exact smallBirthSmokeTwo22

def birthSmokeFamily22 : Fin 17 → Fin 11 → ℝ :=
  profileFamilyBirth22 smallBirthSmokeCaps22 firstCaps zeroHoldCaps endpointDen

theorem birthSmokeFamily22_sound (j : Fin 17) :
    physicalCone22 (birthSmokeFamily22 j) (turboExactBirthMass22 j birthSmokeZ22) :=
  profileFamilyBirth22_sound birthSmokeZ22_mem birthSmokeZ22_nonneg
    firstCaps secondCaps returnCaps feedbackActiveCaps zeroHoldCaps endpointDen familyChecked
    smallBirthSmokeCaps22 smallBirthSmokeCaps22_sound j

#print axioms birthSmokeFamily22_sound
end Spin.Turbo.ProfileFamilyK20Tilt222
