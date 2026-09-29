import SpinCodes.Structured.DenseOccupationFixedVertexDefs

namespace Spin.Structured.DenseOccupationFixed.B047
set_option maxHeartbeats 0
set_option maxRecDepth 100000
def m : QInput := ⟨8, 5⟩
def c : QInput := ⟨-40647089344673, 250000000000000⟩
def p : QInput := ⟨6627332330147145, 9223372036854775808⟩
def y : QInput := ⟨2036402027419661, 4503599627370496⟩
def radius : QInput := ⟨959786189066151951, 1000000000000000000⟩
def z : QInput := ⟨249836040548743199, 250000000000000000⟩
def a0 : QInput := ⟨2819, 8192000⟩
def a1 : QInput := ⟨12047, 20480000⟩
def x0 : QInput := ⟨2469788513329, 15000000000000⟩
def x1 : QInput := ⟨4271577070219, 25000000000000⟩
def upper : Int := -400000000000000000000000
theorem check_00 : vertexCheck m c p y radius z a0 x0 50 upper = true := by decide
theorem check_01 : vertexCheck m c p y radius z a0 x1 50 upper = true := by decide
theorem check_10 : vertexCheck m c p y radius z a1 x0 50 upper = true := by decide
theorem check_11 : vertexCheck m c p y radius z a1 x1 50 upper = true := by decide
end Spin.Structured.DenseOccupationFixed.B047
