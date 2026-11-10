import SpinCodes.Turbo.PhysicalSmallBirth22NumericCensusDefs
import Lean.Elab.Tactic.Decide

namespace Spin.Turbo

set_option maxRecDepth 100000
set_option maxHeartbeats 2000000

theorem smallBirthPair14_15Histogram22_checked :
    smallBirthCodeHistogram22 (smallBirthPairCodes22 14 15) =
      smallBirthPair14_15Histogram22 := by decide +kernel

end Spin.Turbo
