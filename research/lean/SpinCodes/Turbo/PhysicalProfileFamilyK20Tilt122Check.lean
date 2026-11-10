import SpinCodes.Turbo.PhysicalProfileFamilyDefs22
import SpinCodes.Turbo.PhysicalProfileFamilyK20Tilt122Data

namespace Spin.Turbo.ProfileFamilyK20Tilt122

set_option maxHeartbeats 0 in
theorem familyChecked : profileFamilyCheck22 profileNumericEntries22 zbox
    firstCaps secondCaps returnCaps feedbackActiveCaps zeroHoldCaps endpointDen = true := by
  decide +kernel

#print axioms familyChecked
end Spin.Turbo.ProfileFamilyK20Tilt122
