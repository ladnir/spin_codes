import SpinCodes.Turbo.PhysicalSmallBirth22BudgetDefs
import SpinCodes.Turbo.PhysicalProfileFamilyK18Tilt322Data
import Lean.Elab.Tactic.Decide

namespace Spin.Turbo.ProfileFamilyK18Tilt322
set_option maxRecDepth 100000
set_option maxHeartbeats 2000000

def birthZeroCaps22 : List Int := [endpointDen,0,0,0,0,0,0,0,0,0,0]

theorem smallBirthZeroCaps_checked22 :
    smallBirthFixZeroBudgetCheck22 birthZeroCaps22 endpointDen = true := by
  decide +kernel

theorem smallBirthOneCaps_checked22 :
    smallBirthFixOneBudgetCheck22 zbox birthOneCaps endpointDen = true := by
  decide +kernel

theorem smallBirthTwoCaps_checked22 :
    smallBirthFixTwoBudgetCheck22 zbox birthTwoCaps endpointDen = true := by
  decide +kernel

end Spin.Turbo.ProfileFamilyK18Tilt322
