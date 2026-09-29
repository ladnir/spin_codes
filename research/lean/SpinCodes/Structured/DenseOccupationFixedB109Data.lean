import SpinCodes.Structured.DenseOccupationFixedVertexDefs

namespace Spin.Structured.DenseOccupationFixed.B109
set_option maxHeartbeats 0
set_option maxRecDepth 100000
def m : QInput := ⟨4, 5⟩
def c : QInput := ⟨30683842683431, 1250000000000000⟩
def p : QInput := ⟨5748114866138067, 18446744073709551616⟩
def y : QInput := ⟨700962585195439, 1125899906842624⟩
def radius : QInput := ⟨975653087829892349, 1000000000000000000⟩
def z : QInput := ⟨499803674618173987, 500000000000000000⟩
def a0 : QInput := ⟨1, 10000⟩
def a1 : QInput := ⟨18191, 81920000⟩
def x0 : QInput := ⟨144603079418107, 500000000000000⟩
def x1 : QInput := ⟨663870278625033, 2000000000000000⟩
def upper : Int := -400000000000000000000000
theorem check_00 : vertexCheck m c p y radius z a0 x0 50 upper = true := by decide
theorem check_01 : vertexCheck m c p y radius z a0 x1 50 upper = true := by decide
theorem check_10 : vertexCheck m c p y radius z a1 x0 50 upper = true := by decide
theorem check_11 : vertexCheck m c p y radius z a1 x1 50 upper = true := by decide
end Spin.Structured.DenseOccupationFixed.B109
