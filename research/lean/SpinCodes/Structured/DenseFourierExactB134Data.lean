import SpinCodes.Structured.DenseOccupationFixedVertexDefs

namespace Spin.Structured.DenseFourierExact.B134
open Spin.Structured.DenseOccupationFixed
set_option maxHeartbeats 0
set_option maxRecDepth 100000
def m : QInput := ⟨1, 5⟩
def c : QInput := ⟨15724076464601, 62500000000000⟩
def p : QInput := ⟨2250718002527641, 2251799813685248⟩
def y : QInput := ⟨4326852849591391, 9007199254740992⟩
def radius : QInput := ⟨40532412814737869077750591183413660747546602719725177978649935623573, 1250000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000⟩
def z : QInput := ⟨134344346016041137816724063577, 1000000000000000000000000000000⟩
def a0 : QInput := ⟨5110001, 5120000⟩
def a1 : QInput := ⟨1, 1⟩
def x0 : QInput := ⟨10640568152347, 25000000000000⟩
def x1 : QInput := ⟨46255924070167, 100000000000000⟩
def upper : Int := -400000000000000000000000
theorem check_00 : vertexCheck m c p y radius z a0 x0 50 upper=true := by decide
theorem check_01 : vertexCheck m c p y radius z a0 x1 50 upper=true := by decide
theorem check_10 : vertexCheck m c p y radius z a1 x0 50 upper=true := by decide
theorem check_11 : vertexCheck m c p y radius z a1 x1 50 upper=true := by decide
end Spin.Structured.DenseFourierExact.B134
