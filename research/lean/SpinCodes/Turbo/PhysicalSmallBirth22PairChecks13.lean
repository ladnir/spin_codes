import SpinCodes.Turbo.PhysicalSmallBirth22NumericCensusDefs
import Lean.Elab.Tactic.Decide

namespace Spin.Turbo

set_option maxRecDepth 100000
set_option maxHeartbeats 2000000

theorem smallBirthPair13_14Histogram22_checked :
    smallBirthCodeHistogram22 (smallBirthPairCodes22 13 14) =
      smallBirthPair13_14Histogram22 := by decide +kernel

theorem smallBirthPair13_15Histogram22_checked :
    smallBirthCodeHistogram22 (smallBirthPairCodes22 13 15) =
      smallBirthPair13_15Histogram22 := by decide +kernel

end Spin.Turbo
