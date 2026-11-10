import SpinCodes.Turbo.PhysicalSmallBirth22NumericCensusDefs
import Lean.Elab.Tactic.Decide

namespace Spin.Turbo

set_option maxRecDepth 100000
set_option maxHeartbeats 2000000

theorem smallBirthPair12_13Histogram22_checked :
    smallBirthCodeHistogram22 (smallBirthPairCodes22 12 13) =
      smallBirthPair12_13Histogram22 := by decide +kernel

theorem smallBirthPair12_14Histogram22_checked :
    smallBirthCodeHistogram22 (smallBirthPairCodes22 12 14) =
      smallBirthPair12_14Histogram22 := by decide +kernel

theorem smallBirthPair12_15Histogram22_checked :
    smallBirthCodeHistogram22 (smallBirthPairCodes22 12 15) =
      smallBirthPair12_15Histogram22 := by decide +kernel

end Spin.Turbo
