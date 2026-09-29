import SpinCodes.Structured.DenseOccupationFixedVertexDefs

namespace Spin.Structured.DenseFourierExact.B008
open Spin.Structured.DenseOccupationFixed
set_option maxHeartbeats 0
set_option maxRecDepth 100000
def m : QInput := ⟨833, 500⟩
def c : QInput := ⟨-43311, 250000⟩
def p : QInput := ⟨4311748468409033, 4503599627370496⟩
def y : QInput := ⟨7531809796518669, 18014398509481984⟩
def radius : QInput := ⟨8730523732031219078396629560637140297878066034797938498377780447952241647, 10000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000⟩
def z : QInput := ⟨24463354009787175409700147631, 125000000000000000000000000000⟩
def a0 : QInput := ⟨70001, 80000⟩
def a1 : QInput := ⟨1, 1⟩
def x0 : QInput := ⟨13, 125⟩
def x1 : QInput := ⟨18296026121, 120000000000⟩
def upper : Int := -400000000000000000000000
theorem check_00 : vertexCheck m c p y radius z a0 x0 50 upper=true := by decide
theorem check_01 : vertexCheck m c p y radius z a0 x1 50 upper=true := by decide
theorem check_10 : vertexCheck m c p y radius z a1 x0 50 upper=true := by decide
theorem check_11 : vertexCheck m c p y radius z a1 x1 50 upper=true := by decide
end Spin.Structured.DenseFourierExact.B008
