import SpinCodes.Structured.DenseOccupationFixedVertexDefs

namespace Spin.Structured.DenseFourierExact.B097
open Spin.Structured.DenseOccupationFixed
set_option maxHeartbeats 0
set_option maxRecDepth 100000
def m : QInput := ⟨4, 5⟩
def c : QInput := ⟨30683842683431, 1250000000000000⟩
def p : QInput := ⟨8942891430850125, 9007199254740992⟩
def y : QInput := ⟨8169198738063629, 18014398509481984⟩
def radius : QInput := ⟨287244090634974209000572122268512839050413752083904940194603327845101, 2500000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000⟩
def z : QInput := ⟨144005096727870717787454821967, 1000000000000000000000000000000⟩
def a0 : QInput := ⟨310001, 320000⟩
def a1 : QInput := ⟨1, 1⟩
def x0 : QInput := ⟨144603079418107, 500000000000000⟩
def x1 : QInput := ⟨663870278625033, 2000000000000000⟩
def upper : Int := -400000000000000000000000
theorem check_00 : vertexCheck m c p y radius z a0 x0 50 upper=true := by decide
theorem check_01 : vertexCheck m c p y radius z a0 x1 50 upper=true := by decide
theorem check_10 : vertexCheck m c p y radius z a1 x0 50 upper=true := by decide
theorem check_11 : vertexCheck m c p y radius z a1 x1 50 upper=true := by decide
end Spin.Structured.DenseFourierExact.B097
