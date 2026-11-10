import SpinCodes.Turbo.PhysicalSmallBirth22Budget
import SpinCodes.Turbo.PhysicalSmallBirth22BudgetCheck

noncomputable section
namespace Spin.Turbo.ProfileFamilySmoke22
open Spin.Numeric Spin.Numeric.Fix

def birthSmokeZ22 : ℝ := (zNum : ℝ) / (zDen : ℝ)

theorem birthSmokeZ22_mem : Mem zbox birthSmokeZ22 := by
  unfold zbox birthSmokeZ22
  exact ofFrac_mem (show (0 : Int) < zDen by decide +kernel)

theorem birthSmokeZ22_nonneg : 0 ≤ birthSmokeZ22 := by
  unfold birthSmokeZ22
  exact div_nonneg (by exact_mod_cast (show (0 : Int) ≤ zNum by decide +kernel))
    (by exact_mod_cast (show (0 : Int) ≤ zDen by decide +kernel))

theorem birthZeroCaps22_exact :
    (fun i : Fin 11 => (birthZeroCaps22.getD i.val 0 : ℝ) / (endpointDen : ℝ)) =
      physicalInitialEnvelope22 := by
  funext i
  fin_cases i <;> norm_num [birthZeroCaps22, endpointDen, physicalInitialEnvelope22]

theorem smallBirthSmokeZero22 :
    physicalCone22 (fun i => (birthZeroCaps22.getD i.val 0 : ℝ) / (endpointDen : ℝ))
      (turboExactBirthMass22 0 birthSmokeZ22) := by
  rw [birthZeroCaps22_exact]
  exact physicalSmallBirth22_zero_exact_cone birthSmokeZ22

theorem smallBirthSmokeOne22 :
    physicalCone22 (fun i => (birthOneCaps.getD i.val 0 : ℝ) / (endpointDen : ℝ))
      (turboExactBirthMass22 1 birthSmokeZ22) :=
  smallBirthFixOneBudgetCheck22_cone birthSmokeZ22_mem smallBirthOneCaps_checked22

theorem smallBirthSmokeTwo22 :
    physicalCone22 (fun i => (birthTwoCaps.getD i.val 0 : ℝ) / (endpointDen : ℝ))
      (turboExactBirthMass22 2 birthSmokeZ22) :=
  smallBirthFixTwoBudgetCheck22_cone birthSmokeZ22_mem smallBirthTwoCaps_checked22

#print axioms smallBirthSmokeZero22
#print axioms smallBirthSmokeOne22
#print axioms smallBirthSmokeTwo22
end Spin.Turbo.ProfileFamilySmoke22
