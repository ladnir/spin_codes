import SpinCodes.Structured.DenseOccupationFixedVertexDefs

namespace Spin.Structured.DenseOccupationFixed.B392
set_option maxHeartbeats 0
set_option maxRecDepth 100000
def m : QInput := ⟨-31, 20⟩
def c : QInput := ⟨697977398380873, 500000000000000⟩
def p : QInput := ⟨1339177652678037, 18014398509481984⟩
def y : QInput := ⟨8379785276675909, 9007199254740992⟩
def radius : QInput := ⟨9039846643563, 25000000000000000⟩
def z : QInput := ⟨86042465242671661, 100000000000000000⟩
def a0 : QInput := ⟨10031, 320000⟩
def a1 : QInput := ⟨2003, 32000⟩
def x0 : QInput := ⟨61043983743719, 75000000000000⟩
def x1 : QInput := ⟨20728422929781, 25000000000000⟩
def upper : Int := -400000000000000000000000
theorem check_00 : vertexCheck m c p y radius z a0 x0 50 upper = true := by decide
theorem check_01 : vertexCheck m c p y radius z a0 x1 50 upper = true := by decide
theorem check_10 : vertexCheck m c p y radius z a1 x0 50 upper = true := by decide
theorem check_11 : vertexCheck m c p y radius z a1 x1 50 upper = true := by decide
end Spin.Structured.DenseOccupationFixed.B392
