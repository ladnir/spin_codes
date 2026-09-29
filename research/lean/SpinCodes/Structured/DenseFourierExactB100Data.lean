import SpinCodes.Structured.DenseOccupationFixedVertexDefs

namespace Spin.Structured.DenseFourierExact.B100
open Spin.Structured.DenseOccupationFixed
set_option maxHeartbeats 0
set_option maxRecDepth 100000
def m : QInput := ⟨3, 5⟩
def c : QInput := ⟨909341020092481, 10000000000000000⟩
def p : QInput := ⟨5879720969148561, 18014398509481984⟩
def y : QInput := ⟨67378130503675, 140737488355328⟩
def radius : QInput := ⟨9243429491788077166606364581428006950353831045027948693415211311459581965711100918084072441, 2000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000⟩
def z : QInput := ⟨168700512164821796735716864601, 250000000000000000000000000000⟩
def a0 : QInput := ⟨10003, 40000⟩
def a1 : QInput := ⟨50011, 160000⟩
def x0 : QInput := ⟨663870278625033, 2000000000000000⟩
def x1 : QInput := ⟨755265762055919, 2000000000000000⟩
def upper : Int := -400000000000000000000000
theorem check_00 : vertexCheck m c p y radius z a0 x0 50 upper=true := by decide
theorem check_01 : vertexCheck m c p y radius z a0 x1 50 upper=true := by decide
theorem check_10 : vertexCheck m c p y radius z a1 x0 50 upper=true := by decide
theorem check_11 : vertexCheck m c p y radius z a1 x1 50 upper=true := by decide
end Spin.Structured.DenseFourierExact.B100
