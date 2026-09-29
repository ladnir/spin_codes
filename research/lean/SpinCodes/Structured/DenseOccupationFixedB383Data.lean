import SpinCodes.Structured.DenseOccupationFixedVertexDefs

namespace Spin.Structured.DenseOccupationFixed.B383
set_option maxHeartbeats 0
set_option maxRecDepth 100000
def m : QInput := ⟨-7, 5⟩
def c : QInput := ⟨318466707318577, 250000000000000⟩
def p : QInput := ⟨2618750247837145, 36028797018963968⟩
def y : QInput := ⟨1032579159424779, 1125899906842624⟩
def radius : QInput := ⟨122197359495383, 250000000000000000⟩
def z : QInput := ⟨1081953645390783, 1250000000000000⟩
def a0 : QInput := ⟨10031, 320000⟩
def a1 : QInput := ⟨2003, 32000⟩
def x0 : QInput := ⟨1571386577205347, 2000000000000000⟩
def x1 : QInput := ⟨61043983743719, 75000000000000⟩
def upper : Int := -400000000000000000000000
theorem check_00 : vertexCheck m c p y radius z a0 x0 50 upper = true := by decide
theorem check_01 : vertexCheck m c p y radius z a0 x1 50 upper = true := by decide
theorem check_10 : vertexCheck m c p y radius z a1 x0 50 upper = true := by decide
theorem check_11 : vertexCheck m c p y radius z a1 x1 50 upper = true := by decide
end Spin.Structured.DenseOccupationFixed.B383
