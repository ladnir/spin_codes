import SpinCodes.Structured.DenseOccupationFixedVertexDefs

namespace Spin.Structured.DenseFourierExact.B112
open Spin.Structured.DenseOccupationFixed
set_option maxHeartbeats 0
set_option maxRecDepth 100000
def m : QInput := ⟨2, 5⟩
def c : QInput := ⟨4161516955371, 25000000000000⟩
def p : QInput := ⟨2344639505756275, 4503599627370496⟩
def y : QInput := ⟨3455537353331205, 4503599627370496⟩
def radius : QInput := ⟨10329083692095926720796998368762519907957069698736683489245804063334080601, 10000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000⟩
def z : QInput := ⟨9830168630173641774147740207, 50000000000000000000000000000⟩
def a0 : QInput := ⟨90023, 320000⟩
def a1 : QInput := ⟨50011, 160000⟩
def x0 : QInput := ⟨755265762055919, 2000000000000000⟩
def x1 : QInput := ⟨10640568152347, 25000000000000⟩
def upper : Int := -400000000000000000000000
theorem check_00 : vertexCheck m c p y radius z a0 x0 50 upper=true := by decide
theorem check_01 : vertexCheck m c p y radius z a0 x1 50 upper=true := by decide
theorem check_10 : vertexCheck m c p y radius z a1 x0 50 upper=true := by decide
theorem check_11 : vertexCheck m c p y radius z a1 x1 50 upper=true := by decide
end Spin.Structured.DenseFourierExact.B112
