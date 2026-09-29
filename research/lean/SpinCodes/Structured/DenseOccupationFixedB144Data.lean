import SpinCodes.Structured.DenseOccupationFixedVertexDefs

namespace Spin.Structured.DenseOccupationFixed.B144
set_option maxHeartbeats 0
set_option maxRecDepth 100000
def m : QInput := ⟨1, 5⟩
def c : QInput := ⟨15724076464601, 62500000000000⟩
def p : QInput := ⟨8016278118172317, 4611686018427387904⟩
def y : QInput := ⟨1509092462470363, 2251799813685248⟩
def radius : QInput := ⟨866571231344504417, 1000000000000000000⟩
def z : QInput := ⟨498829233529035243, 500000000000000000⟩
def a0 : QInput := ⟨12047, 20480000⟩
def a1 : QInput := ⟨11023, 10240000⟩
def x0 : QInput := ⟨10640568152347, 25000000000000⟩
def x1 : QInput := ⟨46255924070167, 100000000000000⟩
def upper : Int := -400000000000000000000000
theorem check_00 : vertexCheck m c p y radius z a0 x0 50 upper = true := by decide
theorem check_01 : vertexCheck m c p y radius z a0 x1 50 upper = true := by decide
theorem check_10 : vertexCheck m c p y radius z a1 x0 50 upper = true := by decide
theorem check_11 : vertexCheck m c p y radius z a1 x1 50 upper = true := by decide
end Spin.Structured.DenseOccupationFixed.B144
