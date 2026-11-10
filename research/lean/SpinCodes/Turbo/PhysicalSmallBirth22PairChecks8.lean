import SpinCodes.Turbo.PhysicalSmallBirth22NumericCensusDefs
import Lean.Elab.Tactic.Decide

namespace Spin.Turbo

set_option maxRecDepth 100000
set_option maxHeartbeats 2000000

theorem smallBirthPair8_9Histogram22_checked :
    smallBirthCodeHistogram22 (smallBirthPairCodes22 8 9) =
      smallBirthPair8_9Histogram22 := by decide +kernel

theorem smallBirthPair8_10Histogram22_checked :
    smallBirthCodeHistogram22 (smallBirthPairCodes22 8 10) =
      smallBirthPair8_10Histogram22 := by decide +kernel

theorem smallBirthPair8_11Histogram22_checked :
    smallBirthCodeHistogram22 (smallBirthPairCodes22 8 11) =
      smallBirthPair8_11Histogram22 := by decide +kernel

theorem smallBirthPair8_12Histogram22_checked :
    smallBirthCodeHistogram22 (smallBirthPairCodes22 8 12) =
      smallBirthPair8_12Histogram22 := by decide +kernel

theorem smallBirthPair8_13Histogram22_checked :
    smallBirthCodeHistogram22 (smallBirthPairCodes22 8 13) =
      smallBirthPair8_13Histogram22 := by decide +kernel

theorem smallBirthPair8_14Histogram22_checked :
    smallBirthCodeHistogram22 (smallBirthPairCodes22 8 14) =
      smallBirthPair8_14Histogram22 := by decide +kernel

theorem smallBirthPair8_15Histogram22_checked :
    smallBirthCodeHistogram22 (smallBirthPairCodes22 8 15) =
      smallBirthPair8_15Histogram22 := by decide +kernel

end Spin.Turbo
