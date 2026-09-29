import SpinCodes.Structured.FiberNumericsData.Tables
namespace Spin.Structured.FiberNumerics.Data
set_option maxRecDepth 1000000
set_option maxHeartbeats 0
def values3 : List Int := [341376, 325374, 309876, 294874, 280360, 266326, 252764, 239666, 227024, 214830, 203076, 191754, 180856, 170374, 160300, 150626, 141344, 132446, 123924, 115770, 107976, 100534, 93436, 86674, 80240, 74126, 68324, 62826, 57624, 52710, 48076, 43714, 39616, 35774, 32180, 28826, 25704, 22806, 20124, 17650, 15376, 13294, 11396, 9674, 8120, 6726, 5484, 4386, 3424, 2590, 1876, 1274, 776, 374, 60, -174, -336, -434, -476, -470, -424, -346, -244, -126, 0, 126, 244, 346, 424, 470, 476, 434, 336, 174, -60, -374, -776, -1274, -1876, -2590, -3424, -4386, -5484, -6726, -8120, -9674, -11396, -13294, -15376, -17650, -20124, -22806, -25704, -28826, -32180, -35774, -39616, -43714, -48076, -52710, -57624, -62826, -68324, -74126, -80240, -86674, -93436, -100534, -107976, -115770, -123924, -132446, -141344, -150626, -160300, -170374, -180856, -191754, -203076, -214830, -227024, -239666, -252764, -266326, -280360, -294874, -309876, -325374, -341376]
theorem values3_checked : (List.range 129).map (kraw 3) = values3 := by decide +kernel
theorem signed3_checked : signedSum spectrum values3 = 524288 * (kernels.getD 3 0 : Int) := by decide +kernel
theorem square3_checked : squareSum spectrum values3 = 524288 * (896452 : Int) := by decide +kernel
theorem absolute3_checked : absSum spectrum values3 = 190732352 := by decide +kernel
theorem cap3_checked : (190732352 : Int) < 524288 * ((caps.getD 3 0 : Int) + 1) := by decide +kernel
end Spin.Structured.FiberNumerics.Data
