import SpinCodes.Structured.DenseOccupationFixedVertexDefs

namespace Spin.Structured.DenseFourierExact.B098
open Spin.Structured.DenseOccupationFixed
set_option maxHeartbeats 0
set_option maxRecDepth 100000
def m : QInput := ⟨3, 5⟩
def c : QInput := ⟨909341020092481, 10000000000000000⟩
def p : QInput := ⟨6107708308652481, 36028797018963968⟩
def y : QInput := ⟨758188060117363, 1125899906842624⟩
def radius : QInput := ⟨797665094007456611009247078822630631653328810946373594041110550717171592422312271415912034773, 2500000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000⟩
def z : QInput := ⟨734028908017088197304625964779, 1000000000000000000000000000000⟩
def a0 : QInput := ⟨30029, 320000⟩
def a1 : QInput := ⟨10007, 80000⟩
def x0 : QInput := ⟨663870278625033, 2000000000000000⟩
def x1 : QInput := ⟨755265762055919, 2000000000000000⟩
def upper : Int := -400000000000000000000000
theorem check_00 : vertexCheck m c p y radius z a0 x0 50 upper=true := by decide
theorem check_01 : vertexCheck m c p y radius z a0 x1 50 upper=true := by decide
theorem check_10 : vertexCheck m c p y radius z a1 x0 50 upper=true := by decide
theorem check_11 : vertexCheck m c p y radius z a1 x1 50 upper=true := by decide
end Spin.Structured.DenseFourierExact.B098
