import SpinCodes.Structured.DenseOccupationFixedVertexDefs

namespace Spin.Structured.DenseOccupationFixed.B431
set_option maxHeartbeats 0
set_option maxRecDepth 100000
def m : QInput := ⟨-833, 500⟩
def c : QInput := ⟨373189, 250000⟩
def p : QInput := ⟨538722190003689, 18014398509481984⟩
def y : QInput := ⟨8534980030140455, 9007199254740992⟩
def radius : QInput := ⟨13653203490707487, 500000000000000000⟩
def z : QInput := ⟨233018207922273079, 250000000000000000⟩
def a0 : QInput := ⟨10127, 1280000⟩
def a1 : QInput := ⟨10063, 640000⟩
def x0 : QInput := ⟨101703973879, 120000000000⟩
def x1 : QInput := ⟨112, 125⟩
def upper : Int := -400000000000000000000000
theorem check_00 : vertexCheck m c p y radius z a0 x0 50 upper = true := by decide
theorem check_01 : vertexCheck m c p y radius z a0 x1 50 upper = true := by decide
theorem check_10 : vertexCheck m c p y radius z a1 x0 50 upper = true := by decide
theorem check_11 : vertexCheck m c p y radius z a1 x1 50 upper = true := by decide
end Spin.Structured.DenseOccupationFixed.B431
