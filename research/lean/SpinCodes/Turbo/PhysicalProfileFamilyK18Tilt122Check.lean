import SpinCodes.Turbo.PhysicalProfileFamilyDefs22
import SpinCodes.Turbo.PhysicalProfileFamilyK18Tilt122Data

namespace Spin.Turbo.ProfileFamilyK18Tilt122

set_option maxHeartbeats 0 in
theorem familyChecked : profileFamilyCheck22 profileNumericEntries22 zbox
    firstCaps secondCaps returnCaps feedbackActiveCaps zeroHoldCaps endpointDen = true := by
  decide +kernel

#print axioms familyChecked
end Spin.Turbo.ProfileFamilyK18Tilt122
