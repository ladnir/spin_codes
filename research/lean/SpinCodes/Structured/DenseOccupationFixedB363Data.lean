import SpinCodes.Structured.DenseOccupationFixedVertexDefs

namespace Spin.Structured.DenseOccupationFixed.B363
set_option maxHeartbeats 0
set_option maxRecDepth 100000
def m : QInput := ⟨-1, 1⟩
def c : QInput := ⟨483352921189751, 500000000000000⟩
def p : QInput := ⟨8028335195517133, 288230376151711744⟩
def y : QInput := ⟨8213345682318813, 9007199254740992⟩
def radius : QInput := ⟨797725876831829, 20000000000000000⟩
def z : QInput := ⟨46958368624699149, 50000000000000000⟩
def a0 : QInput := ⟨10127, 1280000⟩
def a1 : QInput := ⟨10063, 640000⟩
def x0 : QInput := ⟨355396920581893, 500000000000000⟩
def x1 : QInput := ⟨1500223291742713, 2000000000000000⟩
def upper : Int := -400000000000000000000000
theorem check_00 : vertexCheck m c p y radius z a0 x0 50 upper = true := by decide
theorem check_01 : vertexCheck m c p y radius z a0 x1 50 upper = true := by decide
theorem check_10 : vertexCheck m c p y radius z a1 x0 50 upper = true := by decide
theorem check_11 : vertexCheck m c p y radius z a1 x1 50 upper = true := by decide
end Spin.Structured.DenseOccupationFixed.B363
