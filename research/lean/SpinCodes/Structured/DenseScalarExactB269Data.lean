import SpinCodes.Structured.DenseScalarExactDefs
import SpinCodes.Structured.DenseOccupationFixedVertexDefs

namespace Spin.Structured.DenseScalarExact.B269
open Spin.Structured.DenseOccupationFixed
set_option maxHeartbeats 0
set_option maxRecDepth 100000
def m : QInput := ⟨-2, 5⟩
def c : QInput := ⟨14161516955371, 25000000000000⟩
def p : QInput := ⟨3876927340355661, 4503599627370496⟩
def y : QInput := ⟨2615784076287461, 4503599627370496⟩
def radius : QInput := ⟨1380574128128984774587701579882385371979670928304752348529863831823, 156250000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000⟩
def z : QInput := ⟨61797752808988751801433432259, 500000000000000000000000000000⟩
def a0 : QInput := ⟨30001, 40000⟩
def a1 : QInput := ⟨70001, 80000⟩
def x0 : QInput := ⟨14359431847653, 25000000000000⟩
def x1 : QInput := ⟨1244734237944081, 2000000000000000⟩
def qn : Int := 10141204801825835630440214666721
def qd : Int := 20282409603651670423947251286016
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
end Spin.Structured.DenseScalarExact.B269
