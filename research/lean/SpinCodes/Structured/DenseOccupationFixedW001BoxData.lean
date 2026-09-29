import SpinCodes.Structured.DenseOccupationFixedVertexDefs

namespace Spin.Structured.DenseOccupationFixed.W001Box
set_option maxHeartbeats 0
set_option maxRecDepth 100000
def m : QInput := ⟨833, 500⟩
def c : QInput := ⟨-43311, 250000⟩
def p : QInput := ⟨7540393141205639, 18446744073709551616⟩
def y : QInput := ⟨3569924084993373, 9007199254740992⟩
def radius : QInput := ⟨61224800996618679, 62500000000000000⟩
def z : QInput := ⟨499835864583549449, 500000000000000000⟩
def a0 : QInput := ⟨18191, 81920000⟩
def a1 : QInput := ⟨2819, 8192000⟩
def x0 : QInput := ⟨13, 125⟩
def x1 : QInput := ⟨18296026121, 120000000000⟩
def upper : Int := -400000000000000000000000
theorem check_00 : vertexCheck m c p y radius z a0 x0 50 upper = true := by decide
theorem check_01 : vertexCheck m c p y radius z a0 x1 50 upper = true := by decide
theorem check_10 : vertexCheck m c p y radius z a1 x0 50 upper = true := by decide
theorem check_11 : vertexCheck m c p y radius z a1 x1 50 upper = true := by decide
end Spin.Structured.DenseOccupationFixed.W001Box
