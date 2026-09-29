import SpinCodes.Structured.FiberNumericsData.Tables
namespace Spin.Structured.FiberNumerics.Data
set_option maxRecDepth 1000000
set_option maxHeartbeats 0
def values0 : List Int := [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
theorem values0_checked : (List.range 129).map (kraw 0) = values0 := by decide +kernel
theorem signed0_checked : signedSum spectrum values0 = 524288 * (kernels.getD 0 0 : Int) := by decide +kernel
theorem square0_checked : squareSum spectrum values0 = 524288 * (1 : Int) := by decide +kernel
theorem absolute0_checked : absSum spectrum values0 = 524288 := by decide +kernel
theorem cap0_checked : polyChoose 128 0 - kernels.getD 0 0 ≤ caps.getD 0 0 := by decide +kernel
end Spin.Structured.FiberNumerics.Data
