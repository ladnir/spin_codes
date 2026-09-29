import SpinCodes.Structured.DenseOccupationFixedVertexDefs

namespace Spin.Structured.DenseOccupationFixed.B101
set_option maxHeartbeats 0
set_option maxRecDepth 100000
def m : QInput := ⟨1, 1⟩
def c : QInput := ⟨-16647078810249, 500000000000000⟩
def p : QInput := ⟨2931632439007657, 1152921504606846976⟩
def y : QInput := ⟨8916579636254071, 18014398509481984⟩
def radius : QInput := ⟨214243213565521611, 250000000000000000⟩
def z : QInput := ⟨249367656860749659, 250000000000000000⟩
def a0 : QInput := ⟨11023, 10240000⟩
def a1 : QInput := ⟨10511, 5120000⟩
def x0 : QInput := ⟨499776708257287, 2000000000000000⟩
def x1 : QInput := ⟨144603079418107, 500000000000000⟩
def upper : Int := -400000000000000000000000
theorem check_00 : vertexCheck m c p y radius z a0 x0 50 upper = true := by decide
theorem check_01 : vertexCheck m c p y radius z a0 x1 50 upper = true := by decide
theorem check_10 : vertexCheck m c p y radius z a1 x0 50 upper = true := by decide
theorem check_11 : vertexCheck m c p y radius z a1 x1 50 upper = true := by decide
end Spin.Structured.DenseOccupationFixed.B101
