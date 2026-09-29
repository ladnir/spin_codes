import SpinCodes.Structured.FiberNumericsData.Tables
namespace Spin.Structured.FiberNumerics.Data
set_option maxRecDepth 1000000
set_option maxHeartbeats 0
def values1 : List Int := [128, 126, 124, 122, 120, 118, 116, 114, 112, 110, 108, 106, 104, 102, 100, 98, 96, 94, 92, 90, 88, 86, 84, 82, 80, 78, 76, 74, 72, 70, 68, 66, 64, 62, 60, 58, 56, 54, 52, 50, 48, 46, 44, 42, 40, 38, 36, 34, 32, 30, 28, 26, 24, 22, 20, 18, 16, 14, 12, 10, 8, 6, 4, 2, 0, -2, -4, -6, -8, -10, -12, -14, -16, -18, -20, -22, -24, -26, -28, -30, -32, -34, -36, -38, -40, -42, -44, -46, -48, -50, -52, -54, -56, -58, -60, -62, -64, -66, -68, -70, -72, -74, -76, -78, -80, -82, -84, -86, -88, -90, -92, -94, -96, -98, -100, -102, -104, -106, -108, -110, -112, -114, -116, -118, -120, -122, -124, -126, -128]
theorem values1_checked : (List.range 129).map (kraw 1) = values1 := by decide +kernel
theorem signed1_checked : signedSum spectrum values1 = 524288 * (kernels.getD 1 0 : Int) := by decide +kernel
theorem square1_checked : squareSum spectrum values1 = 524288 * (128 : Int) := by decide +kernel
theorem absolute1_checked : absSum spectrum values1 = 4723464 := by decide +kernel
theorem cap1_checked : polyChoose 128 (1 - 1) < (caps.getD 1 0 + 1) * polyChoose 1 (1 - 1) := by decide +kernel
end Spin.Structured.FiberNumerics.Data
