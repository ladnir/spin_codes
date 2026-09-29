import SpinCodes.Structured.DenseOccupationFixedVertexDefs

namespace Spin.Structured.DenseFourierExact.B019
open Spin.Structured.DenseOccupationFixed
set_option maxHeartbeats 0
set_option maxRecDepth 100000
def m : QInput := ⟨33, 20⟩
def c : QInput := ⟨-5335812508647, 31250000000000⟩
def p : QInput := ⟨4264963108698761, 9007199254740992⟩
def y : QInput := ⟨2606388136873927, 9007199254740992⟩
def radius : QInput := ⟨17094777305611659726923639659099494811522178265763228104564402482039399636818462233894509211, 400000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000⟩
def z : QInput := ⟨716822628492396305761323849149, 1000000000000000000000000000000⟩
def a0 : QInput := ⟨6001, 16000⟩
def a1 : QInput := ⟨10001, 20000⟩
def x0 : QInput := ⟨791599208623, 5000000000000⟩
def x1 : QInput := ⟨1609032935677, 10000000000000⟩
def upper : Int := -400000000000000000000000
theorem check_00 : vertexCheck m c p y radius z a0 x0 50 upper=true := by decide
theorem check_01 : vertexCheck m c p y radius z a0 x1 50 upper=true := by decide
theorem check_10 : vertexCheck m c p y radius z a1 x0 50 upper=true := by decide
theorem check_11 : vertexCheck m c p y radius z a1 x1 50 upper=true := by decide
end Spin.Structured.DenseFourierExact.B019
