import SpinCodes.Structured.DenseOccupationFixedVertexDefs

namespace Spin.Structured.DenseOccupationFixed.B029
set_option maxHeartbeats 0
set_option maxRecDepth 100000
def m : QInput := ⟨33, 20⟩
def c : QInput := ⟨-5335812508647, 31250000000000⟩
def p : QInput := ⟨4555512470969881, 576460752303423488⟩
def y : QInput := ⟨6195950020969639, 18014398509481984⟩
def radius : QInput := ⟨28414558065447697, 40000000000000000⟩
def z : QInput := ⟨496998138849668249, 500000000000000000⟩
def a0 : QInput := ⟨2051, 512000⟩
def a1 : QInput := ⟨10127, 1280000⟩
def x0 : QInput := ⟨791599208623, 5000000000000⟩
def x1 : QInput := ⟨1609032935677, 10000000000000⟩
def upper : Int := -400000000000000000000000
theorem check_00 : vertexCheck m c p y radius z a0 x0 50 upper = true := by decide
theorem check_01 : vertexCheck m c p y radius z a0 x1 50 upper = true := by decide
theorem check_10 : vertexCheck m c p y radius z a1 x0 50 upper = true := by decide
theorem check_11 : vertexCheck m c p y radius z a1 x1 50 upper = true := by decide
end Spin.Structured.DenseOccupationFixed.B029
