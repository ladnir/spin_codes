import SpinCodes.Structured.DenseScalarExactDefs
import SpinCodes.Structured.DenseOccupationFixedVertexDefs

namespace Spin.Structured.DenseScalarExact.B278
open Spin.Structured.DenseOccupationFixed
set_option maxHeartbeats 0
set_option maxRecDepth 100000
def m : QInput := ⟨-8, 5⟩
def c : QInput := ⟨359352910655327, 250000000000000⟩
def p : QInput := ⟨369596523894813, 562949953421312⟩
def y : QInput := ⟨53590981916819, 70368744177664⟩
def radius : QInput := ⟨88356744200254992901881972608164275106581052014752335787063697777997, 10000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000⟩
def z : QInput := ⟨61797752808988751801433432259, 500000000000000000000000000000⟩
def a0 : QInput := ⟨10001, 20000⟩
def a1 : QInput := ⟨30001, 40000⟩
def x0 : QInput := ⟨20728422929781, 25000000000000⟩
def x1 : QInput := ⟨12530211486671, 15000000000000⟩
def qn : Int := 19807040628566084922271559847
def qd : Int := 39614081257132168796771975168
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
end Spin.Structured.DenseScalarExact.B278
