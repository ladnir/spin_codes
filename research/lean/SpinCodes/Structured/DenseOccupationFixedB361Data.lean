import SpinCodes.Structured.DenseOccupationFixedVertexDefs

namespace Spin.Structured.DenseOccupationFixed.B361
set_option maxHeartbeats 0
set_option maxRecDepth 100000
def m : QInput := ⟨-1, 1⟩
def c : QInput := ⟨483352921189751, 500000000000000⟩
def p : QInput := ⟨5582818842050355, 1152921504606846976⟩
def y : QInput := ⟨7842102082714381, 9007199254740992⟩
def radius : QInput := ⟨122276118556298217, 200000000000000000⟩
def z : QInput := ⟨991543744621664727, 1000000000000000000⟩
def a0 : QInput := ⟨10511, 5120000⟩
def a1 : QInput := ⟨2051, 512000⟩
def x0 : QInput := ⟨355396920581893, 500000000000000⟩
def x1 : QInput := ⟨1500223291742713, 2000000000000000⟩
def upper : Int := -400000000000000000000000
theorem check_00 : vertexCheck m c p y radius z a0 x0 50 upper = true := by decide
theorem check_01 : vertexCheck m c p y radius z a0 x1 50 upper = true := by decide
theorem check_10 : vertexCheck m c p y radius z a1 x0 50 upper = true := by decide
theorem check_11 : vertexCheck m c p y radius z a1 x1 50 upper = true := by decide
end Spin.Structured.DenseOccupationFixed.B361
