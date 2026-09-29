import SpinCodes.Structured.DenseOccupationFixedVertexDefs

namespace Spin.Structured.DenseFourierExact.B234
open Spin.Structured.DenseOccupationFixed
set_option maxHeartbeats 0
set_option maxRecDepth 100000
def m : QInput := ⟨-2, 5⟩
def c : QInput := ⟨14161516955371, 25000000000000⟩
def p : QInput := ⟨4201953848697771, 9007199254740992⟩
def y : QInput := ⟨3924165104999497, 4503599627370496⟩
def radius : QInput := ⟨2403793593969702679630908632069952025522171253417999442584873230686735839, 10000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000⟩
def z : QInput := ⟨23565361038062345912536605661, 125000000000000000000000000000⟩
def a0 : QInput := ⟨30013, 160000⟩
def a1 : QInput := ⟨10003, 40000⟩
def x0 : QInput := ⟨14359431847653, 25000000000000⟩
def x1 : QInput := ⟨1244734237944081, 2000000000000000⟩
def upper : Int := -400000000000000000000000
theorem check_00 : vertexCheck m c p y radius z a0 x0 50 upper=true := by decide
theorem check_01 : vertexCheck m c p y radius z a0 x1 50 upper=true := by decide
theorem check_10 : vertexCheck m c p y radius z a1 x0 50 upper=true := by decide
theorem check_11 : vertexCheck m c p y radius z a1 x1 50 upper=true := by decide
end Spin.Structured.DenseFourierExact.B234
