import SpinCodes.Structured.DenseOccupationFixedVertexDefs

namespace Spin.Structured.DenseFourierExact.B027
open Spin.Structured.DenseOccupationFixed
set_option maxHeartbeats 0
set_option maxRecDepth 100000
def m : QInput := ⟨163, 100⟩
def c : QInput := ⟨-3350558688107, 20000000000000⟩
def p : QInput := ⟨5314823625655089, 9007199254740992⟩
def y : QInput := ⟨2396444708353143, 9007199254740992⟩
def radius : QInput := ⟨21138063203464848155164257505128425222797082132350763035630434223125604319371973372978516033, 5000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000⟩
def z : QInput := ⟨673036769763149052622117888201, 1000000000000000000000000000000⟩
def a0 : QInput := ⟨10001, 20000⟩
def a1 : QInput := ⟨50003, 80000⟩
def x0 : QInput := ⟨1609032935677, 10000000000000⟩
def x1 : QInput := ⟨2469788513329, 15000000000000⟩
def upper : Int := -400000000000000000000000
theorem check_00 : vertexCheck m c p y radius z a0 x0 50 upper=true := by decide
theorem check_01 : vertexCheck m c p y radius z a0 x1 50 upper=true := by decide
theorem check_10 : vertexCheck m c p y radius z a1 x0 50 upper=true := by decide
theorem check_11 : vertexCheck m c p y radius z a1 x1 50 upper=true := by decide
end Spin.Structured.DenseFourierExact.B027
