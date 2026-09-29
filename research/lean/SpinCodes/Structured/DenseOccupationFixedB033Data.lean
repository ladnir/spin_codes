import SpinCodes.Structured.DenseOccupationFixedVertexDefs

namespace Spin.Structured.DenseOccupationFixed.B033
set_option maxHeartbeats 0
set_option maxRecDepth 100000
def m : QInput := ⟨33, 20⟩
def c : QInput := ⟨-5335812508647, 31250000000000⟩
def p : QInput := ⟨2278504734861077, 18014398509481984⟩
def y : QInput := ⟨455741800043343, 1125899906842624⟩
def radius : QInput := ⟨1379591990555021, 1000000000000000000⟩
def z : QInput := ⟨219252283879008421, 250000000000000000⟩
def a0 : QInput := ⟨2003, 32000⟩
def a1 : QInput := ⟨10007, 80000⟩
def x0 : QInput := ⟨791599208623, 5000000000000⟩
def x1 : QInput := ⟨1609032935677, 10000000000000⟩
def upper : Int := -400000000000000000000000
theorem check_00 : vertexCheck m c p y radius z a0 x0 50 upper = true := by decide
theorem check_01 : vertexCheck m c p y radius z a0 x1 50 upper = true := by decide
theorem check_10 : vertexCheck m c p y radius z a1 x0 50 upper = true := by decide
theorem check_11 : vertexCheck m c p y radius z a1 x1 50 upper = true := by decide
end Spin.Structured.DenseOccupationFixed.B033
