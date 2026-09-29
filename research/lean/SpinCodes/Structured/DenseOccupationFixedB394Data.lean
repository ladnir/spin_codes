import SpinCodes.Structured.DenseOccupationFixedVertexDefs

namespace Spin.Structured.DenseOccupationFixed.B394
set_option maxHeartbeats 0
set_option maxRecDepth 100000
def m : QInput := ⟨-8, 5⟩
def c : QInput := ⟨359352910655327, 250000000000000⟩
def p : QInput := ⟨7785299271488899, 9223372036854775808⟩
def y : QInput := ⟨8533120385928193, 9007199254740992⟩
def radius : QInput := ⟨45271811707643113, 50000000000000000⟩
def z : QInput := ⟨499195554386816231, 500000000000000000⟩
def a0 : QInput := ⟨2819, 8192000⟩
def a1 : QInput := ⟨12047, 20480000⟩
def x0 : QInput := ⟨20728422929781, 25000000000000⟩
def x1 : QInput := ⟨12530211486671, 15000000000000⟩
def upper : Int := -400000000000000000000000
theorem check_00 : vertexCheck m c p y radius z a0 x0 50 upper = true := by decide
theorem check_01 : vertexCheck m c p y radius z a0 x1 50 upper = true := by decide
theorem check_10 : vertexCheck m c p y radius z a1 x0 50 upper = true := by decide
theorem check_11 : vertexCheck m c p y radius z a1 x1 50 upper = true := by decide
end Spin.Structured.DenseOccupationFixed.B394
