import SpinCodes.Structured.DenseOccupationFixedVertexDefs

namespace Spin.Structured.DenseOccupationFixed.B069
set_option maxHeartbeats 0
set_option maxRecDepth 100000
def m : QInput := ⟨31, 20⟩
def c : QInput := ⟨-77022601619127, 500000000000000⟩
def p : QInput := ⟨2393778308360057, 18014398509481984⟩
def y : QInput := ⟨1760693787441903, 4503599627370496⟩
def radius : QInput := ⟨1249065029511317, 1000000000000000000⟩
def z : QInput := ⟨437600365093662993, 500000000000000000⟩
def a0 : QInput := ⟨30029, 320000⟩
def a1 : QInput := ⟨10007, 80000⟩
def x0 : QInput := ⟨4271577070219, 25000000000000⟩
def x1 : QInput := ⟨13956016256281, 75000000000000⟩
def upper : Int := -400000000000000000000000
theorem check_00 : vertexCheck m c p y radius z a0 x0 50 upper = true := by decide
theorem check_01 : vertexCheck m c p y radius z a0 x1 50 upper = true := by decide
theorem check_10 : vertexCheck m c p y radius z a1 x0 50 upper = true := by decide
theorem check_11 : vertexCheck m c p y radius z a1 x1 50 upper = true := by decide
end Spin.Structured.DenseOccupationFixed.B069
