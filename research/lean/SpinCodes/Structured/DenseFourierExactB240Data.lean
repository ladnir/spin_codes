import SpinCodes.Structured.DenseOccupationFixedVertexDefs

namespace Spin.Structured.DenseFourierExact.B240
open Spin.Structured.DenseOccupationFixed
set_option maxHeartbeats 0
set_option maxRecDepth 100000
def m : QInput := ⟨-3, 5⟩
def c : QInput := ⟨6909341020092481, 10000000000000000⟩
def p : QInput := ⟨2770717136971445, 18014398509481984⟩
def y : QInput := ⟨7623923779171179, 9007199254740992⟩
def radius : QInput := ⟨29764630106772185940615773792302755677813347449435925180673960280476255944492752246457327561, 312500000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000⟩
def z : QInput := ⟨365576028424162091352714126547, 500000000000000000000000000000⟩
def a0 : QInput := ⟨2003, 32000⟩
def a1 : QInput := ⟨10007, 80000⟩
def x0 : QInput := ⟨1244734237944081, 2000000000000000⟩
def x1 : QInput := ⟨1336129721374967, 2000000000000000⟩
def upper : Int := -400000000000000000000000
theorem check_00 : vertexCheck m c p y radius z a0 x0 50 upper=true := by decide
theorem check_01 : vertexCheck m c p y radius z a0 x1 50 upper=true := by decide
theorem check_10 : vertexCheck m c p y radius z a1 x0 50 upper=true := by decide
theorem check_11 : vertexCheck m c p y radius z a1 x1 50 upper=true := by decide
end Spin.Structured.DenseFourierExact.B240
