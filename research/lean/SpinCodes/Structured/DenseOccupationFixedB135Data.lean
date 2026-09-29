import SpinCodes.Structured.DenseOccupationFixedVertexDefs

namespace Spin.Structured.DenseOccupationFixed.B135
set_option maxHeartbeats 0
set_option maxRecDepth 100000
def m : QInput := ⟨2, 5⟩
def c : QInput := ⟨4161516955371, 25000000000000⟩
def p : QInput := ⟨3508088791386291, 1152921504606846976⟩
def y : QInput := ⟨5628901625793451, 9007199254740992⟩
def radius : QInput := ⟨794714617649193959, 1000000000000000000⟩
def z : QInput := ⟨996181689204542863, 1000000000000000000⟩
def a0 : QInput := ⟨11023, 10240000⟩
def a1 : QInput := ⟨10511, 5120000⟩
def x0 : QInput := ⟨755265762055919, 2000000000000000⟩
def x1 : QInput := ⟨10640568152347, 25000000000000⟩
def upper : Int := -400000000000000000000000
theorem check_00 : vertexCheck m c p y radius z a0 x0 50 upper = true := by decide
theorem check_01 : vertexCheck m c p y radius z a0 x1 50 upper = true := by decide
theorem check_10 : vertexCheck m c p y radius z a1 x0 50 upper = true := by decide
theorem check_11 : vertexCheck m c p y radius z a1 x1 50 upper = true := by decide
end Spin.Structured.DenseOccupationFixed.B135
