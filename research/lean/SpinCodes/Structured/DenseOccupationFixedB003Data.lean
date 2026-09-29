import SpinCodes.Structured.DenseOccupationFixedVertexDefs

namespace Spin.Structured.DenseOccupationFixed.B003
set_option maxHeartbeats 0
set_option maxRecDepth 100000
def m : QInput := ⟨833, 500⟩
def c : QInput := ⟨-43311, 250000⟩
def p : QInput := ⟨1389286033866889, 1152921504606846976⟩
def y : QInput := ⟨7170850984949569, 18014398509481984⟩
def radius : QInput := ⟨941546053044503591, 1000000000000000000⟩
def z : QInput := ⟨249758325674315079, 250000000000000000⟩
def a0 : QInput := ⟨12047, 20480000⟩
def a1 : QInput := ⟨11023, 10240000⟩
def x0 : QInput := ⟨13, 125⟩
def x1 : QInput := ⟨18296026121, 120000000000⟩
def upper : Int := -400000000000000000000000
theorem check_00 : vertexCheck m c p y radius z a0 x0 50 upper = true := by decide
theorem check_01 : vertexCheck m c p y radius z a0 x1 50 upper = true := by decide
theorem check_10 : vertexCheck m c p y radius z a1 x0 50 upper = true := by decide
theorem check_11 : vertexCheck m c p y radius z a1 x1 50 upper = true := by decide
end Spin.Structured.DenseOccupationFixed.B003
