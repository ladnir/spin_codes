import SpinCodes.Structured.DenseOccupationFixedVertexDefs

namespace Spin.Structured.DenseFourierExact.B061
open Spin.Structured.DenseOccupationFixed
set_option maxHeartbeats 0
set_option maxRecDepth 100000
def m : QInput := ⟨6, 5⟩
def c : QInput := ⟨-832718284462267, 10000000000000000⟩
def p : QInput := ⟨8814024575881773, 18014398509481984⟩
def y : QInput := ⟨755933253628365, 2251799813685248⟩
def radius : QInput := ⟨1043667339675262170145543964747671300732936648755369707246574359543048417041170243952528613, 625000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000⟩
def z : QInput := ⟨654077080551653369762611097397, 1000000000000000000000000000000⟩
def a0 : QInput := ⟨70009, 160000⟩
def a1 : QInput := ⟨150017, 320000⟩
def x0 : QInput := ⟨428613422794653, 2000000000000000⟩
def x1 : QInput := ⟨499776708257287, 2000000000000000⟩
def upper : Int := -400000000000000000000000
theorem check_00 : vertexCheck m c p y radius z a0 x0 50 upper=true := by decide
theorem check_01 : vertexCheck m c p y radius z a0 x1 50 upper=true := by decide
theorem check_10 : vertexCheck m c p y radius z a1 x0 50 upper=true := by decide
theorem check_11 : vertexCheck m c p y radius z a1 x1 50 upper=true := by decide
end Spin.Structured.DenseFourierExact.B061
