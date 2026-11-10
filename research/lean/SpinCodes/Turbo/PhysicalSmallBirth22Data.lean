import SpinCodes.Turbo.PhysicalSmallBirth22Encoding

namespace Spin.Turbo
open Spin.UniformMixer Spin.Structured

set_option maxRecDepth 100000 in
set_option maxHeartbeats 1000000 in
theorem smallBirthSingleImages22_checked : ∀ a : SmallBirthOne22,
    PackedMap.eval productionPackedRows (PackedMap.eval smallBirthFeedbackRows22
      (smallBirthOnePacked22 a)) = smallBirthSingleImages22.getD (15*a.1.val+a.2.val) 0 := by
  decide +kernel

end Spin.Turbo
