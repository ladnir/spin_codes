import SpinCodes.Structured.DenseOccupationFixedVertexDefs

namespace Spin.Structured.DenseFourierExact.B121
open Spin.Structured.DenseOccupationFixed
set_option maxHeartbeats 0
set_option maxRecDepth 100000
def m : QInput := ⟨2, 5⟩
def c : QInput := ⟨4161516955371, 25000000000000⟩
def p : QInput := ⟨8990173165116617, 9007199254740992⟩
def y : QInput := ⟨8321922620456365, 18014398509481984⟩
def radius : QInput := ⟨161692158644144152003292020206135035172141909141910187917564932116019, 2500000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000⟩
def z : QInput := ⟨139947924148734641672489112037, 1000000000000000000000000000000⟩
def a0 : QInput := ⟨1270001, 1280000⟩
def a1 : QInput := ⟨1, 1⟩
def x0 : QInput := ⟨755265762055919, 2000000000000000⟩
def x1 : QInput := ⟨10640568152347, 25000000000000⟩
def upper : Int := -400000000000000000000000
theorem check_00 : vertexCheck m c p y radius z a0 x0 50 upper=true := by decide
theorem check_01 : vertexCheck m c p y radius z a0 x1 50 upper=true := by decide
theorem check_10 : vertexCheck m c p y radius z a1 x0 50 upper=true := by decide
theorem check_11 : vertexCheck m c p y radius z a1 x1 50 upper=true := by decide
end Spin.Structured.DenseFourierExact.B121
