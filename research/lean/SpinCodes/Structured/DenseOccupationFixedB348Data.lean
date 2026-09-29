import SpinCodes.Structured.DenseOccupationFixedVertexDefs

namespace Spin.Structured.DenseOccupationFixed.B348
set_option maxHeartbeats 0
set_option maxRecDepth 100000
def m : QInput := ⟨-4, 5⟩
def c : QInput := ⟨1030683842683431, 1250000000000000⟩
def p : QInput := ⟨8443108539554299, 18446744073709551616⟩
def y : QInput := ⟨1878038348029621, 2251799813685248⟩
def radius : QInput := ⟨29782059060146459, 31250000000000000⟩
def z : QInput := ⟨999229990462954797, 1000000000000000000⟩
def a0 : QInput := ⟨1, 10000⟩
def a1 : QInput := ⟨2819, 8192000⟩
def x0 : QInput := ⟨1336129721374967, 2000000000000000⟩
def x1 : QInput := ⟨355396920581893, 500000000000000⟩
def upper : Int := -400000000000000000000000
theorem check_00 : vertexCheck m c p y radius z a0 x0 50 upper = true := by decide
theorem check_01 : vertexCheck m c p y radius z a0 x1 50 upper = true := by decide
theorem check_10 : vertexCheck m c p y radius z a1 x0 50 upper = true := by decide
theorem check_11 : vertexCheck m c p y radius z a1 x1 50 upper = true := by decide
end Spin.Structured.DenseOccupationFixed.B348
