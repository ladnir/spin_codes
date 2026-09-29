import SpinCodes.Structured.DenseOccupationFixedVertexDefs

namespace Spin.Structured.DenseOccupationFixed.B087
set_option maxHeartbeats 0
set_option maxRecDepth 100000
def m : QInput := ⟨6, 5⟩
def c : QInput := ⟨-832718284462267, 10000000000000000⟩
def p : QInput := ⟨6234582211406925, 4611686018427387904⟩
def y : QInput := ⟨8907424680921167, 18014398509481984⟩
def radius : QInput := ⟨114999093717063149, 125000000000000000⟩
def z : QInput := ⟨99865433968626971, 100000000000000000⟩
def a0 : QInput := ⟨12047, 20480000⟩
def a1 : QInput := ⟨11023, 10240000⟩
def x0 : QInput := ⟨428613422794653, 2000000000000000⟩
def x1 : QInput := ⟨499776708257287, 2000000000000000⟩
def upper : Int := -400000000000000000000000
theorem check_00 : vertexCheck m c p y radius z a0 x0 50 upper = true := by decide
theorem check_01 : vertexCheck m c p y radius z a0 x1 50 upper = true := by decide
theorem check_10 : vertexCheck m c p y radius z a1 x0 50 upper = true := by decide
theorem check_11 : vertexCheck m c p y radius z a1 x1 50 upper = true := by decide
end Spin.Structured.DenseOccupationFixed.B087
