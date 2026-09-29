import SpinCodes.Structured.FiberNumericsData.Tables
namespace Spin.Structured.FiberNumerics.Data
set_option maxRecDepth 1000000
set_option maxHeartbeats 0
def values126 : List Int := [8128, -7874, 7624, -7378, 7136, -6898, 6664, -6434, 6208, -5986, 5768, -5554, 5344, -5138, 4936, -4738, 4544, -4354, 4168, -3986, 3808, -3634, 3464, -3298, 3136, -2978, 2824, -2674, 2528, -2386, 2248, -2114, 1984, -1858, 1736, -1618, 1504, -1394, 1288, -1186, 1088, -994, 904, -818, 736, -658, 584, -514, 448, -386, 328, -274, 224, -178, 136, -98, 64, -34, 8, 14, -32, 46, -56, 62, -64, 62, -56, 46, -32, 14, 8, -34, 64, -98, 136, -178, 224, -274, 328, -386, 448, -514, 584, -658, 736, -818, 904, -994, 1088, -1186, 1288, -1394, 1504, -1618, 1736, -1858, 1984, -2114, 2248, -2386, 2528, -2674, 2824, -2978, 3136, -3298, 3464, -3634, 3808, -3986, 4168, -4354, 4544, -4738, 4936, -5138, 5344, -5554, 5768, -5986, 6208, -6434, 6664, -6898, 7136, -7378, 7624, -7874, 8128]
theorem values126_checked : (List.range 129).map (kraw 126) = values126 := by decide +kernel
theorem signed126_checked : signedSum spectrum values126 = 524288 * (kernels.getD 126 0 : Int) := by decide +kernel
theorem square126_checked : squareSum spectrum values126 = 524288 * (8572 : Int) := by decide +kernel
theorem absolute126_checked : absSum spectrum values126 = 32477552 := by decide +kernel
theorem cap126_checked : (32477552 : Int) < 524288 * ((caps.getD 126 0 : Int) + 1) := by decide +kernel
end Spin.Structured.FiberNumerics.Data
