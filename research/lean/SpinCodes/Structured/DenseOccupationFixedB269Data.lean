import SpinCodes.Structured.DenseOccupationFixedVertexDefs

namespace Spin.Structured.DenseOccupationFixed.B269
set_option maxHeartbeats 0
set_option maxRecDepth 100000
def m : QInput := ⟨-3, 100⟩
def c : QInput := ⟨3076970006775960790657142966285747138917, 8507059173023461586584365185794205286400⟩
def p : QInput := ⟨7062211145660689, 9223372036854775808⟩
def y : QInput := ⟨912594572419959, 1125899906842624⟩
def radius : QInput := ⟨231346326113961457, 250000000000000000⟩
def z : QInput := ⟨998750352847232903, 1000000000000000000⟩
def a0 : QInput := ⟨2819, 8192000⟩
def a1 : QInput := ⟨12047, 20480000⟩
def x0 : QInput := ⟨253319723456283, 500000000000000⟩
def x1 : QInput := ⟨43385881287743932506819717224474801189, 85070591730234615865843651857942052864⟩
def upper : Int := -400000000000000000000000
theorem check_00 : vertexCheck m c p y radius z a0 x0 50 upper = true := by decide
theorem check_01 : vertexCheck m c p y radius z a0 x1 50 upper = true := by decide
theorem check_10 : vertexCheck m c p y radius z a1 x0 50 upper = true := by decide
theorem check_11 : vertexCheck m c p y radius z a1 x1 50 upper = true := by decide
end Spin.Structured.DenseOccupationFixed.B269
