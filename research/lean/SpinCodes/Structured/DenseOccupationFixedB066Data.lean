import SpinCodes.Structured.DenseOccupationFixedVertexDefs

namespace Spin.Structured.DenseOccupationFixed.B066
set_option maxHeartbeats 0
set_option maxRecDepth 100000
def m : QInput := ⟨31, 20⟩
def c : QInput := ⟨-77022601619127, 500000000000000⟩
def p : QInput := ⟨6486000215290855, 144115188075855872⟩
def y : QInput := ⟨3972993769894303, 9007199254740992⟩
def radius : QInput := ⟨20326765402718051, 250000000000000000⟩
def z : QInput := ⟨476288882392710319, 500000000000000000⟩
def a0 : QInput := ⟨10031, 320000⟩
def a1 : QInput := ⟨30061, 640000⟩
def x0 : QInput := ⟨4271577070219, 25000000000000⟩
def x1 : QInput := ⟨13956016256281, 75000000000000⟩
def upper : Int := -400000000000000000000000
theorem check_00 : vertexCheck m c p y radius z a0 x0 50 upper = true := by decide
theorem check_01 : vertexCheck m c p y radius z a0 x1 50 upper = true := by decide
theorem check_10 : vertexCheck m c p y radius z a1 x0 50 upper = true := by decide
theorem check_11 : vertexCheck m c p y radius z a1 x1 50 upper = true := by decide
end Spin.Structured.DenseOccupationFixed.B066
