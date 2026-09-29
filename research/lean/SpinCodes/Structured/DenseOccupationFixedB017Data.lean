import SpinCodes.Structured.DenseOccupationFixedVertexDefs

namespace Spin.Structured.DenseOccupationFixed.B017
set_option maxHeartbeats 0
set_option maxRecDepth 100000
def m : QInput := ⟨83, 50⟩
def c : QInput := ⟨-3446583973879, 20000000000000⟩
def p : QInput := ⟨91863826404957, 18014398509481984⟩
def y : QInput := ⟨5106789601943699, 9007199254740992⟩
def radius : QInput := ⟨695129352285639459, 1000000000000000000⟩
def z : QInput := ⟨993593053369497721, 1000000000000000000⟩
def a0 : QInput := ⟨10511, 5120000⟩
def a1 : QInput := ⟨2051, 512000⟩
def x0 : QInput := ⟨18296026121, 120000000000⟩
def x1 : QInput := ⟨791599208623, 5000000000000⟩
def upper : Int := -400000000000000000000000
theorem check_00 : vertexCheck m c p y radius z a0 x0 50 upper = true := by decide
theorem check_01 : vertexCheck m c p y radius z a0 x1 50 upper = true := by decide
theorem check_10 : vertexCheck m c p y radius z a1 x0 50 upper = true := by decide
theorem check_11 : vertexCheck m c p y radius z a1 x1 50 upper = true := by decide
end Spin.Structured.DenseOccupationFixed.B017
