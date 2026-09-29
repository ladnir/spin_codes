import SpinCodes.Structured.DenseOccupationFixedVertexDefs

namespace Spin.Structured.DenseFourierExact.B117
open Spin.Structured.DenseOccupationFixed
set_option maxHeartbeats 0
set_option maxRecDepth 100000
def m : QInput := ⟨2, 5⟩
def c : QInput := ⟨4161516955371, 25000000000000⟩
def p : QInput := ⟨4151870571704245, 4503599627370496⟩
def y : QInput := ⟨66366420288691, 140737488355328⟩
def radius : QInput := ⟨1218619883833225501920407611056466526626755833729032041663784628743049, 5000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000⟩
def z : QInput := ⟨71062275730418897364172289921, 500000000000000000000000000000⟩
def a0 : QInput := ⟨70001, 80000⟩
def a1 : QInput := ⟨150001, 160000⟩
def x0 : QInput := ⟨755265762055919, 2000000000000000⟩
def x1 : QInput := ⟨10640568152347, 25000000000000⟩
def upper : Int := -400000000000000000000000
theorem check_00 : vertexCheck m c p y radius z a0 x0 50 upper=true := by decide
theorem check_01 : vertexCheck m c p y radius z a0 x1 50 upper=true := by decide
theorem check_10 : vertexCheck m c p y radius z a1 x0 50 upper=true := by decide
theorem check_11 : vertexCheck m c p y radius z a1 x1 50 upper=true := by decide
end Spin.Structured.DenseFourierExact.B117
