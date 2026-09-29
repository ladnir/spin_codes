import SpinCodes.Structured.DenseOccupationFixedVertexDefs

namespace Spin.Structured.DenseOccupationFixed.B080
set_option maxHeartbeats 0
set_option maxRecDepth 100000
def m : QInput := ⟨7, 5⟩
def c : QInput := ⟨-31533292681423, 250000000000000⟩
def p : QInput := ⟨847458784944225, 18014398509481984⟩
def y : QInput := ⟨4367999538656629, 9007199254740992⟩
def radius : QInput := ⟨55617715715190807, 1000000000000000000⟩
def z : QInput := ⟨945414332506776949, 1000000000000000000⟩
def a0 : QInput := ⟨10031, 320000⟩
def a1 : QInput := ⟨30061, 640000⟩
def x0 : QInput := ⟨13956016256281, 75000000000000⟩
def x1 : QInput := ⟨428613422794653, 2000000000000000⟩
def upper : Int := -400000000000000000000000
theorem check_00 : vertexCheck m c p y radius z a0 x0 50 upper = true := by decide
theorem check_01 : vertexCheck m c p y radius z a0 x1 50 upper = true := by decide
theorem check_10 : vertexCheck m c p y radius z a1 x0 50 upper = true := by decide
theorem check_11 : vertexCheck m c p y radius z a1 x1 50 upper = true := by decide
end Spin.Structured.DenseOccupationFixed.B080
