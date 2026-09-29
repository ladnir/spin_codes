import SpinCodes.Structured.DenseOccupationFixedVertexDefs

namespace Spin.Structured.DenseFourierExact.B250
open Spin.Structured.DenseOccupationFixed
set_option maxHeartbeats 0
set_option maxRecDepth 100000
def m : QInput := ⟨-4, 5⟩
def c : QInput := ⟨1030683842683431, 1250000000000000⟩
def p : QInput := ⟨8889321548273121, 18014398509481984⟩
def y : QInput := ⟨3948150018169933, 4503599627370496⟩
def radius : QInput := ⟨4137523218079900317807461647298854423036418421429725753374593223394697, 10000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000⟩
def z : QInput := ⟨36517657921075685689843698701, 250000000000000000000000000000⟩
def a0 : QInput := ⟨10003, 40000⟩
def a1 : QInput := ⟨10001, 20000⟩
def x0 : QInput := ⟨1336129721374967, 2000000000000000⟩
def x1 : QInput := ⟨355396920581893, 500000000000000⟩
def upper : Int := -400000000000000000000000
theorem check_00 : vertexCheck m c p y radius z a0 x0 50 upper=true := by decide
theorem check_01 : vertexCheck m c p y radius z a0 x1 50 upper=true := by decide
theorem check_10 : vertexCheck m c p y radius z a1 x0 50 upper=true := by decide
theorem check_11 : vertexCheck m c p y radius z a1 x1 50 upper=true := by decide
end Spin.Structured.DenseFourierExact.B250
