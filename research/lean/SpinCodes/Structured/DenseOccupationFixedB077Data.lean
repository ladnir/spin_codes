import SpinCodes.Structured.DenseOccupationFixedVertexDefs

namespace Spin.Structured.DenseOccupationFixed.B077
set_option maxHeartbeats 0
set_option maxRecDepth 100000
def m : QInput := ⟨7, 5⟩
def c : QInput := ⟨-31533292681423, 250000000000000⟩
def p : QInput := ⟨5024844544481359, 288230376151711744⟩
def y : QInput := ⟨8303147941204327, 18014398509481984⟩
def radius : QInput := ⟨18231841431496869, 50000000000000000⟩
def z : QInput := ⟨61328681698187213, 62500000000000000⟩
def a0 : QInput := ⟨10127, 1280000⟩
def a1 : QInput := ⟨10063, 640000⟩
def x0 : QInput := ⟨13956016256281, 75000000000000⟩
def x1 : QInput := ⟨428613422794653, 2000000000000000⟩
def upper : Int := -400000000000000000000000
theorem check_00 : vertexCheck m c p y radius z a0 x0 50 upper = true := by decide
theorem check_01 : vertexCheck m c p y radius z a0 x1 50 upper = true := by decide
theorem check_10 : vertexCheck m c p y radius z a1 x0 50 upper = true := by decide
theorem check_11 : vertexCheck m c p y radius z a1 x1 50 upper = true := by decide
end Spin.Structured.DenseOccupationFixed.B077
