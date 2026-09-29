import SpinCodes.Structured.DenseOccupationFixedVertexDefs

namespace Spin.Structured.DenseOccupationFixed.B105
set_option maxHeartbeats 0
set_option maxRecDepth 100000
def m : QInput := ⟨1, 1⟩
def c : QInput := ⟨-16647078810249, 500000000000000⟩
def p : QInput := ⟨656333912706797, 18014398509481984⟩
def y : QInput := ⟨4598273613442935, 9007199254740992⟩
def radius : QInput := ⟨19084683862117739, 200000000000000000⟩
def z : QInput := ⟨5972549877036271, 6250000000000000⟩
def a0 : QInput := ⟨10063, 640000⟩
def a1 : QInput := ⟨10031, 320000⟩
def x0 : QInput := ⟨499776708257287, 2000000000000000⟩
def x1 : QInput := ⟨144603079418107, 500000000000000⟩
def upper : Int := -400000000000000000000000
theorem check_00 : vertexCheck m c p y radius z a0 x0 50 upper = true := by decide
theorem check_01 : vertexCheck m c p y radius z a0 x1 50 upper = true := by decide
theorem check_10 : vertexCheck m c p y radius z a1 x0 50 upper = true := by decide
theorem check_11 : vertexCheck m c p y radius z a1 x1 50 upper = true := by decide
end Spin.Structured.DenseOccupationFixed.B105
