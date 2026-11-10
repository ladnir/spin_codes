import SpinCodes.Turbo.PhysicalSmallBirth22NumericCensusDefs
import Lean.Elab.Tactic.Decide

namespace Spin.Turbo

set_option maxRecDepth 100000
set_option maxHeartbeats 2000000

theorem smallBirthPair10_11Histogram22_checked :
    smallBirthCodeHistogram22 (smallBirthPairCodes22 10 11) =
      smallBirthPair10_11Histogram22 := by decide +kernel

theorem smallBirthPair10_12Histogram22_checked :
    smallBirthCodeHistogram22 (smallBirthPairCodes22 10 12) =
      smallBirthPair10_12Histogram22 := by decide +kernel

theorem smallBirthPair10_13Histogram22_checked :
    smallBirthCodeHistogram22 (smallBirthPairCodes22 10 13) =
      smallBirthPair10_13Histogram22 := by decide +kernel

theorem smallBirthPair10_14Histogram22_checked :
    smallBirthCodeHistogram22 (smallBirthPairCodes22 10 14) =
      smallBirthPair10_14Histogram22 := by decide +kernel

theorem smallBirthPair10_15Histogram22_checked :
    smallBirthCodeHistogram22 (smallBirthPairCodes22 10 15) =
      smallBirthPair10_15Histogram22 := by decide +kernel

end Spin.Turbo
