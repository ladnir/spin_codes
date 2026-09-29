import SpinCodes.Structured.DenseOccupationFixedVertexDefs

namespace Spin.Structured.DenseFourierExact.B032
open Spin.Structured.DenseOccupationFixed
set_option maxHeartbeats 0
set_option maxRecDepth 100000
def m : QInput := ⟨8, 5⟩
def c : QInput := ⟨-40647089344673, 250000000000000⟩
def p : QInput := ⟨4503553296132189, 18014398509481984⟩
def y : QInput := ⟨3783632554279913, 9007199254740992⟩
def radius : QInput := ⟨10640915657438657150776743915865135360954599468047588247315000507055253096432073245330919779171, 10000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000⟩
def z : QInput := ⟨750451989816520167880247345847, 1000000000000000000000000000000⟩
def a0 : QInput := ⟨30013, 160000⟩
def a1 : QInput := ⟨10003, 40000⟩
def x0 : QInput := ⟨2469788513329, 15000000000000⟩
def x1 : QInput := ⟨4271577070219, 25000000000000⟩
def upper : Int := -400000000000000000000000
theorem check_00 : vertexCheck m c p y radius z a0 x0 50 upper=true := by decide
theorem check_01 : vertexCheck m c p y radius z a0 x1 50 upper=true := by decide
theorem check_10 : vertexCheck m c p y radius z a1 x0 50 upper=true := by decide
theorem check_11 : vertexCheck m c p y radius z a1 x1 50 upper=true := by decide
end Spin.Structured.DenseFourierExact.B032
