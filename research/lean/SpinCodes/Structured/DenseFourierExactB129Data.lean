import SpinCodes.Structured.DenseOccupationFixedVertexDefs

namespace Spin.Structured.DenseFourierExact.B129
open Spin.Structured.DenseOccupationFixed
set_option maxHeartbeats 0
set_option maxRecDepth 100000
def m : QInput := ⟨1, 5⟩
def c : QInput := ⟨15724076464601, 62500000000000⟩
def p : QInput := ⟨8585446259202839, 9007199254740992⟩
def y : QInput := ⟨8788432720174529, 18014398509481984⟩
def radius : QInput := ⟨138308204133989368685187375815212277675408492435252023415088047842633, 10000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000⟩
def z : QInput := ⟨61809742059976736250288177081, 500000000000000000000000000000⟩
def a0 : QInput := ⟨150001, 160000⟩
def a1 : QInput := ⟨310001, 320000⟩
def x0 : QInput := ⟨10640568152347, 25000000000000⟩
def x1 : QInput := ⟨46255924070167, 100000000000000⟩
def upper : Int := -400000000000000000000000
theorem check_00 : vertexCheck m c p y radius z a0 x0 50 upper=true := by decide
theorem check_01 : vertexCheck m c p y radius z a0 x1 50 upper=true := by decide
theorem check_10 : vertexCheck m c p y radius z a1 x0 50 upper=true := by decide
theorem check_11 : vertexCheck m c p y radius z a1 x1 50 upper=true := by decide
end Spin.Structured.DenseFourierExact.B129
