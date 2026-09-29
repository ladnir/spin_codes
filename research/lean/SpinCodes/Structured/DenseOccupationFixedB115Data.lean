import SpinCodes.Structured.DenseOccupationFixedVertexDefs

namespace Spin.Structured.DenseOccupationFixed.B115
set_option maxHeartbeats 0
set_option maxRecDepth 100000
def m : QInput := ⟨4, 5⟩
def c : QInput := ⟨30683842683431, 1250000000000000⟩
def p : QInput := ⟨1416392720674615, 144115188075855872⟩
def y : QInput := ⟨2516854545074777, 4503599627370496⟩
def radius : QInput := ⟨501872902558936359, 1000000000000000000⟩
def z : QInput := ⟨987411147655164997, 1000000000000000000⟩
def a0 : QInput := ⟨2051, 512000⟩
def a1 : QInput := ⟨10127, 1280000⟩
def x0 : QInput := ⟨144603079418107, 500000000000000⟩
def x1 : QInput := ⟨663870278625033, 2000000000000000⟩
def upper : Int := -400000000000000000000000
theorem check_00 : vertexCheck m c p y radius z a0 x0 50 upper = true := by decide
theorem check_01 : vertexCheck m c p y radius z a0 x1 50 upper = true := by decide
theorem check_10 : vertexCheck m c p y radius z a1 x0 50 upper = true := by decide
theorem check_11 : vertexCheck m c p y radius z a1 x1 50 upper = true := by decide
end Spin.Structured.DenseOccupationFixed.B115
