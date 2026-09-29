import SpinCodes.Structured.DenseOccupationFixedVertexDefs

namespace Spin.Structured.DenseFourierExact.B069
open Spin.Structured.DenseOccupationFixed
set_option maxHeartbeats 0
set_option maxRecDepth 100000
def m : QInput := ⟨1, 1⟩
def c : QInput := ⟨-16647078810249, 500000000000000⟩
def p : QInput := ⟨4508382918673613, 18014398509481984⟩
def y : QInput := ⟨8445944435854149, 18014398509481984⟩
def radius : QInput := ⟨2097128403179466924668705376047236427358820910170529033272571862996891784572294634949427689789, 10000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000⟩
def z : QInput := ⟨727509670879533099545798794729, 1000000000000000000000000000000⟩
def a0 : QInput := ⟨30013, 160000⟩
def a1 : QInput := ⟨10003, 40000⟩
def x0 : QInput := ⟨499776708257287, 2000000000000000⟩
def x1 : QInput := ⟨144603079418107, 500000000000000⟩
def upper : Int := -400000000000000000000000
theorem check_00 : vertexCheck m c p y radius z a0 x0 50 upper=true := by decide
theorem check_01 : vertexCheck m c p y radius z a0 x1 50 upper=true := by decide
theorem check_10 : vertexCheck m c p y radius z a1 x0 50 upper=true := by decide
theorem check_11 : vertexCheck m c p y radius z a1 x1 50 upper=true := by decide
end Spin.Structured.DenseFourierExact.B069
