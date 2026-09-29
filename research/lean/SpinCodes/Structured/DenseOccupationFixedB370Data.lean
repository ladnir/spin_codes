import SpinCodes.Structured.DenseOccupationFixedVertexDefs

namespace Spin.Structured.DenseOccupationFixed.B370
set_option maxHeartbeats 0
set_option maxRecDepth 100000
def m : QInput := ⟨-6, 5⟩
def c : QInput := ⟨11167281715537733, 10000000000000000⟩
def p : QInput := ⟨5877428163472163, 1152921504606846976⟩
def y : QInput := ⟨4014916974657849, 4503599627370496⟩
def radius : QInput := ⟨147356476440371319, 250000000000000000⟩
def z : QInput := ⟨990885090766869147, 1000000000000000000⟩
def a0 : QInput := ⟨10511, 5120000⟩
def a1 : QInput := ⟨2051, 512000⟩
def x0 : QInput := ⟨1500223291742713, 2000000000000000⟩
def x1 : QInput := ⟨1571386577205347, 2000000000000000⟩
def upper : Int := -400000000000000000000000
theorem check_00 : vertexCheck m c p y radius z a0 x0 50 upper = true := by decide
theorem check_01 : vertexCheck m c p y radius z a0 x1 50 upper = true := by decide
theorem check_10 : vertexCheck m c p y radius z a1 x0 50 upper = true := by decide
theorem check_11 : vertexCheck m c p y radius z a1 x1 50 upper = true := by decide
end Spin.Structured.DenseOccupationFixed.B370
