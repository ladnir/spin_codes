import SpinCodes.Structured.DenseOccupationFixedVertexDefs

namespace Spin.Structured.DenseOccupationFixed.B402
set_option maxHeartbeats 0
set_option maxRecDepth 100000
def m : QInput := ⟨-163, 100⟩
def c : QInput := ⟨29249441311893, 20000000000000⟩
def p : QInput := ⟨3970861265985637, 4611686018427387904⟩
def y : QInput := ⟨8592891644846165, 9007199254740992⟩
def radius : QInput := ⟨903051364244867331, 1000000000000000000⟩
def z : QInput := ⟨62396713257614457, 62500000000000000⟩
def a0 : QInput := ⟨2819, 8192000⟩
def a1 : QInput := ⟨12047, 20480000⟩
def x0 : QInput := ⟨12530211486671, 15000000000000⟩
def x1 : QInput := ⟨8390967064323, 10000000000000⟩
def upper : Int := -400000000000000000000000
theorem check_00 : vertexCheck m c p y radius z a0 x0 50 upper = true := by decide
theorem check_01 : vertexCheck m c p y radius z a0 x1 50 upper = true := by decide
theorem check_10 : vertexCheck m c p y radius z a1 x0 50 upper = true := by decide
theorem check_11 : vertexCheck m c p y radius z a1 x1 50 upper = true := by decide
end Spin.Structured.DenseOccupationFixed.B402
