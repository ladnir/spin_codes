import SpinCodes.Structured.DenseOccupationFixedVertexDefs

namespace Spin.Structured.DenseOccupationFixed.B379
set_option maxHeartbeats 0
set_option maxRecDepth 100000
def m : QInput := ⟨-7, 5⟩
def c : QInput := ⟨318466707318577, 250000000000000⟩
def p : QInput := ⟨3077632761682905, 576460752303423488⟩
def y : QInput := ⟨8190449866801671, 9007199254740992⟩
def radius : QInput := ⟨569494515032899949, 1000000000000000000⟩
def z : QInput := ⟨990263899332569973, 1000000000000000000⟩
def a0 : QInput := ⟨10511, 5120000⟩
def a1 : QInput := ⟨2051, 512000⟩
def x0 : QInput := ⟨1571386577205347, 2000000000000000⟩
def x1 : QInput := ⟨61043983743719, 75000000000000⟩
def upper : Int := -400000000000000000000000
theorem check_00 : vertexCheck m c p y radius z a0 x0 50 upper = true := by decide
theorem check_01 : vertexCheck m c p y radius z a0 x1 50 upper = true := by decide
theorem check_10 : vertexCheck m c p y radius z a1 x0 50 upper = true := by decide
theorem check_11 : vertexCheck m c p y radius z a1 x1 50 upper = true := by decide
end Spin.Structured.DenseOccupationFixed.B379
