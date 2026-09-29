import SpinCodes.Structured.DenseOccupationFixedVertexDefs

namespace Spin.Structured.DenseFourierExact.B107
open Spin.Structured.DenseOccupationFixed
set_option maxHeartbeats 0
set_option maxRecDepth 100000
def m : QInput := ⟨3, 5⟩
def c : QInput := ⟨909341020092481, 10000000000000000⟩
def p : QInput := ⟨8647679663855109, 9007199254740992⟩
def y : QInput := ⟨4093913618003789, 9007199254740992⟩
def radius : QInput := ⟨1649667967354308707547242832331479311630012919078782378469466335943607, 10000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000⟩
def z : QInput := ⟨139163099236396518512507940269, 1000000000000000000000000000000⟩
def a0 : QInput := ⟨150001, 160000⟩
def a1 : QInput := ⟨310001, 320000⟩
def x0 : QInput := ⟨663870278625033, 2000000000000000⟩
def x1 : QInput := ⟨755265762055919, 2000000000000000⟩
def upper : Int := -400000000000000000000000
theorem check_00 : vertexCheck m c p y radius z a0 x0 50 upper=true := by decide
theorem check_01 : vertexCheck m c p y radius z a0 x1 50 upper=true := by decide
theorem check_10 : vertexCheck m c p y radius z a1 x0 50 upper=true := by decide
theorem check_11 : vertexCheck m c p y radius z a1 x1 50 upper=true := by decide
end Spin.Structured.DenseFourierExact.B107
