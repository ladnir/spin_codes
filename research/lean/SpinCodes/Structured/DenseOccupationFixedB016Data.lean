import SpinCodes.Structured.DenseOccupationFixedVertexDefs

namespace Spin.Structured.DenseOccupationFixed.B016
set_option maxHeartbeats 0
set_option maxRecDepth 100000
def m : QInput := ⟨83, 50⟩
def c : QInput := ⟨-3446583973879, 20000000000000⟩
def p : QInput := ⟨653360498700495, 288230376151711744⟩
def y : QInput := ⟨3590912238627717, 9007199254740992⟩
def radius : QInput := ⟨89412267416196007, 100000000000000000⟩
def z : QInput := ⟨998182373560825953, 1000000000000000000⟩
def a0 : QInput := ⟨11023, 10240000⟩
def a1 : QInput := ⟨10511, 5120000⟩
def x0 : QInput := ⟨18296026121, 120000000000⟩
def x1 : QInput := ⟨791599208623, 5000000000000⟩
def upper : Int := -400000000000000000000000
theorem check_00 : vertexCheck m c p y radius z a0 x0 50 upper = true := by decide
theorem check_01 : vertexCheck m c p y radius z a0 x1 50 upper = true := by decide
theorem check_10 : vertexCheck m c p y radius z a1 x0 50 upper = true := by decide
theorem check_11 : vertexCheck m c p y radius z a1 x1 50 upper = true := by decide
end Spin.Structured.DenseOccupationFixed.B016
