import SpinCodes.Structured.DenseOccupationFixedVertexDefs

namespace Spin.Structured.DenseOccupationFixed.B128
set_option maxHeartbeats 0
set_option maxRecDepth 100000
def m : QInput := ⟨3, 5⟩
def c : QInput := ⟨909341020092481, 10000000000000000⟩
def p : QInput := ⟨5867204301919767, 144115188075855872⟩
def y : QInput := ⟨2740253334654435, 4503599627370496⟩
def radius : QInput := ⟨8648306053745819, 200000000000000000⟩
def z : QInput := ⟨188136601309538867, 200000000000000000⟩
def a0 : QInput := ⟨10063, 640000⟩
def a1 : QInput := ⟨10031, 320000⟩
def x0 : QInput := ⟨663870278625033, 2000000000000000⟩
def x1 : QInput := ⟨755265762055919, 2000000000000000⟩
def upper : Int := -400000000000000000000000
theorem check_00 : vertexCheck m c p y radius z a0 x0 50 upper = true := by decide
theorem check_01 : vertexCheck m c p y radius z a0 x1 50 upper = true := by decide
theorem check_10 : vertexCheck m c p y radius z a1 x0 50 upper = true := by decide
theorem check_11 : vertexCheck m c p y radius z a1 x1 50 upper = true := by decide
end Spin.Structured.DenseOccupationFixed.B128
