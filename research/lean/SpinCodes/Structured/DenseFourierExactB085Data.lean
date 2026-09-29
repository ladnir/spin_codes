import SpinCodes.Structured.DenseOccupationFixedVertexDefs

namespace Spin.Structured.DenseFourierExact.B085
open Spin.Structured.DenseOccupationFixed
set_option maxHeartbeats 0
set_option maxRecDepth 100000
def m : QInput := ⟨4, 5⟩
def c : QInput := ⟨30683842683431, 1250000000000000⟩
def p : QInput := ⟨7228206018573285, 18014398509481984⟩
def y : QInput := ⟨7963494514188715, 18014398509481984⟩
def radius : QInput := ⟨1209634286068883774938325866772688773949518740802873216568366265271096360447382042232979033, 5000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000⟩
def z : QInput := ⟨152703764196646332145502129861, 250000000000000000000000000000⟩
def a0 : QInput := ⟨110021, 320000⟩
def a1 : QInput := ⟨230041, 640000⟩
def x0 : QInput := ⟨144603079418107, 500000000000000⟩
def x1 : QInput := ⟨663870278625033, 2000000000000000⟩
def upper : Int := -400000000000000000000000
theorem check_00 : vertexCheck m c p y radius z a0 x0 50 upper=true := by decide
theorem check_01 : vertexCheck m c p y radius z a0 x1 50 upper=true := by decide
theorem check_10 : vertexCheck m c p y radius z a1 x0 50 upper=true := by decide
theorem check_11 : vertexCheck m c p y radius z a1 x1 50 upper=true := by decide
end Spin.Structured.DenseFourierExact.B085
