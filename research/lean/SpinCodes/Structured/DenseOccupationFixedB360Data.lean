import SpinCodes.Structured.DenseOccupationFixedVertexDefs

namespace Spin.Structured.DenseOccupationFixed.B360
set_option maxHeartbeats 0
set_option maxRecDepth 100000
def m : QInput := ⟨-1, 1⟩
def c : QInput := ⟨483352921189751, 500000000000000⟩
def p : QInput := ⟨6735083947156095, 2305843009213693952⟩
def y : QInput := ⟨1029718593646931, 1125899906842624⟩
def radius : QInput := ⟨142905463235814773, 200000000000000000⟩
def z : QInput := ⟨248526174362229249, 250000000000000000⟩
def a0 : QInput := ⟨11023, 10240000⟩
def a1 : QInput := ⟨10511, 5120000⟩
def x0 : QInput := ⟨355396920581893, 500000000000000⟩
def x1 : QInput := ⟨1500223291742713, 2000000000000000⟩
def upper : Int := -400000000000000000000000
theorem check_00 : vertexCheck m c p y radius z a0 x0 50 upper = true := by decide
theorem check_01 : vertexCheck m c p y radius z a0 x1 50 upper = true := by decide
theorem check_10 : vertexCheck m c p y radius z a1 x0 50 upper = true := by decide
theorem check_11 : vertexCheck m c p y radius z a1 x1 50 upper = true := by decide
end Spin.Structured.DenseOccupationFixed.B360
