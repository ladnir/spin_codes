import SpinCodes.Turbo.PhysicalSmallBirth22NumericCensusDefs
import Lean.Elab.Tactic.Decide

namespace Spin.Turbo

set_option maxRecDepth 100000
set_option maxHeartbeats 2000000

theorem smallBirthPair11_12Histogram22_checked :
    smallBirthCodeHistogram22 (smallBirthPairCodes22 11 12) =
      smallBirthPair11_12Histogram22 := by decide +kernel

theorem smallBirthPair11_13Histogram22_checked :
    smallBirthCodeHistogram22 (smallBirthPairCodes22 11 13) =
      smallBirthPair11_13Histogram22 := by decide +kernel

theorem smallBirthPair11_14Histogram22_checked :
    smallBirthCodeHistogram22 (smallBirthPairCodes22 11 14) =
      smallBirthPair11_14Histogram22 := by decide +kernel

theorem smallBirthPair11_15Histogram22_checked :
    smallBirthCodeHistogram22 (smallBirthPairCodes22 11 15) =
      smallBirthPair11_15Histogram22 := by decide +kernel

end Spin.Turbo
