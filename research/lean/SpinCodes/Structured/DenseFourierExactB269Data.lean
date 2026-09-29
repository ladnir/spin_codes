import SpinCodes.Structured.DenseOccupationFixedVertexDefs

namespace Spin.Structured.DenseFourierExact.B269
open Spin.Structured.DenseOccupationFixed
set_option maxHeartbeats 0
set_option maxRecDepth 100000
def m : QInput := ⟨-7, 5⟩
def c : QInput := ⟨318466707318577, 250000000000000⟩
def p : QInput := ⟨8026065456636501, 18014398509481984⟩
def y : QInput := ⟨512712222824705, 562949953421312⟩
def radius : QInput := ⟨699546618592847885535376737068669977591586810056757388830256149704998231, 2500000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000⟩
def z : QInput := ⟨189398664564211694713966292461, 1000000000000000000000000000000⟩
def a0 : QInput := ⟨10003, 40000⟩
def a1 : QInput := ⟨10001, 20000⟩
def x0 : QInput := ⟨1571386577205347, 2000000000000000⟩
def x1 : QInput := ⟨61043983743719, 75000000000000⟩
def upper : Int := -400000000000000000000000
theorem check_00 : vertexCheck m c p y radius z a0 x0 50 upper=true := by decide
theorem check_01 : vertexCheck m c p y radius z a0 x1 50 upper=true := by decide
theorem check_10 : vertexCheck m c p y radius z a1 x0 50 upper=true := by decide
theorem check_11 : vertexCheck m c p y radius z a1 x1 50 upper=true := by decide
end Spin.Structured.DenseFourierExact.B269
