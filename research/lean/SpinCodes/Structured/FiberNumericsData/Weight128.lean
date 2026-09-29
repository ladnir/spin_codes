import SpinCodes.Structured.FiberNumericsData.Tables
namespace Spin.Structured.FiberNumerics.Data
set_option maxRecDepth 1000000
set_option maxHeartbeats 0
def values128 : List Int := [1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1]
theorem values128_checked : (List.range 129).map (kraw 128) = values128 := by decide +kernel
theorem signed128_checked : signedSum spectrum values128 = 524288 * (kernels.getD 128 0 : Int) := by decide +kernel
theorem square128_checked : squareSum spectrum values128 = 524288 * (1 : Int) := by decide +kernel
theorem absolute128_checked : absSum spectrum values128 = 524288 := by decide +kernel
theorem cap128_checked : polyChoose 128 128 - kernels.getD 128 0 ≤ caps.getD 128 0 := by decide +kernel
end Spin.Structured.FiberNumerics.Data
