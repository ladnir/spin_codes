import SpinCodes.Structured.DenseScalarExactDefs
import SpinCodes.Structured.DenseOccupationFixedVertexDefs

namespace Spin.Structured.DenseScalarExact.B017
open Spin.Structured.DenseOccupationFixed
set_option maxHeartbeats 0
set_option maxRecDepth 100000
def m : QInput := ⟨1, 20⟩
def c : QInput := ⟨109535552428011023053662565657799459079, 340282366920938463463374607431768211456⟩
def p : QInput := ⟨1125866225087155, 1125899906842624⟩
def y : QInput := ⟨2251867179211453, 4503599627370496⟩
def radius : QInput := ⟨22089186050063953819393210903837891424160990007679430767157352719427, 2500000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000⟩
def z : QInput := ⟨61797752808988751801433432259, 500000000000000000000000000000⟩
def a0 : QInput := ⟨655330003, 655360000⟩
def a1 : QInput := ⟨327670001, 327680000⟩
def x0 : QInput := ⟨487324427761897, 1000000000000000⟩
def x1 : QInput := ⟨41684710442490683359023934633467251675, 85070591730234615865843651857942052864⟩
def qn : Int := 2535301200456458549879099186215
def qd : Int := 5070602400912917605986812821504
theorem check0 : check qn qd z.num z.den radius.num radius.den 0 = true := by decide
theorem check48 : check qn qd z.num z.den radius.num radius.den 48 = true := by decide
theorem check56 : check qn qd z.num z.den radius.num radius.den 56 = true := by decide
theorem check64 : check qn qd z.num z.den radius.num radius.den 64 = true := by decide
theorem check72 : check qn qd z.num z.den radius.num radius.den 72 = true := by decide
theorem check80 : check qn qd z.num z.den radius.num radius.den 80 = true := by decide
def upper : Int := -400000000000000000000000
theorem check_00 : vertexCheck m c p y radius z a0 x0 50 upper = true := by decide
theorem check_01 : vertexCheck m c p y radius z a0 x1 50 upper = true := by decide
theorem check_10 : vertexCheck m c p y radius z a1 x0 50 upper = true := by decide
theorem check_11 : vertexCheck m c p y radius z a1 x1 50 upper = true := by decide
end Spin.Structured.DenseScalarExact.B017
