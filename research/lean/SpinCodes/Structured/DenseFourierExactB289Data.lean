import SpinCodes.Structured.DenseOccupationFixedVertexDefs

namespace Spin.Structured.DenseFourierExact.B289
open Spin.Structured.DenseOccupationFixed
set_option maxHeartbeats 0
set_option maxRecDepth 100000
def m : QInput := ⟨-33, 20⟩
def c : QInput := ⟨46226687491353, 31250000000000⟩
def p : QInput := ⟨2041405660094129, 18014398509481984⟩
def y : QInput := ⟨4280891984665623, 4503599627370496⟩
def radius : QInput := ⟨4964811301507357600640598277402803503069849600324427111811818118774527523160562963028293977909, 2500000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000⟩
def z : QInput := ⟨779127604144557393674790673843, 1000000000000000000000000000000⟩
def a0 : QInput := ⟨10031, 320000⟩
def a1 : QInput := ⟨2003, 32000⟩
def x0 : QInput := ⟨8390967064323, 10000000000000⟩
def x1 : QInput := ⟨4208400791377, 5000000000000⟩
def upper : Int := -400000000000000000000000
theorem check_00 : vertexCheck m c p y radius z a0 x0 50 upper=true := by decide
theorem check_01 : vertexCheck m c p y radius z a0 x1 50 upper=true := by decide
theorem check_10 : vertexCheck m c p y radius z a1 x0 50 upper=true := by decide
theorem check_11 : vertexCheck m c p y radius z a1 x1 50 upper=true := by decide
end Spin.Structured.DenseFourierExact.B289
