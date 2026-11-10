import SpinCodes.Turbo.PhysicalSmallBirth22NumericCensusDefs
import Lean.Elab.Tactic.Decide

namespace Spin.Turbo

set_option maxRecDepth 100000
set_option maxHeartbeats 2000000

theorem smallBirthPair9_10Histogram22_checked :
    smallBirthCodeHistogram22 (smallBirthPairCodes22 9 10) =
      smallBirthPair9_10Histogram22 := by decide +kernel

theorem smallBirthPair9_11Histogram22_checked :
    smallBirthCodeHistogram22 (smallBirthPairCodes22 9 11) =
      smallBirthPair9_11Histogram22 := by decide +kernel

theorem smallBirthPair9_12Histogram22_checked :
    smallBirthCodeHistogram22 (smallBirthPairCodes22 9 12) =
      smallBirthPair9_12Histogram22 := by decide +kernel

theorem smallBirthPair9_13Histogram22_checked :
    smallBirthCodeHistogram22 (smallBirthPairCodes22 9 13) =
      smallBirthPair9_13Histogram22 := by decide +kernel

theorem smallBirthPair9_14Histogram22_checked :
    smallBirthCodeHistogram22 (smallBirthPairCodes22 9 14) =
      smallBirthPair9_14Histogram22 := by decide +kernel

theorem smallBirthPair9_15Histogram22_checked :
    smallBirthCodeHistogram22 (smallBirthPairCodes22 9 15) =
      smallBirthPair9_15Histogram22 := by decide +kernel

end Spin.Turbo
