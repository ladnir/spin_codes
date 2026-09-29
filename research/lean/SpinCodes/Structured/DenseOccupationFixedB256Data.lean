import SpinCodes.Structured.DenseOccupationFixedVertexDefs

namespace Spin.Structured.DenseOccupationFixed.B256
set_option maxHeartbeats 0
set_option maxRecDepth 100000
def m : QInput := ⟨-3, 100⟩
def c : QInput := ⟨3076970006775960790657142966285747138917, 8507059173023461586584365185794205286400⟩
def p : QInput := ⟨1655419989492613, 18014398509481984⟩
def y : QInput := ⟨3468465479991229, 4503599627370496⟩
def radius : QInput := ⟨103366818434587, 1000000000000000000⟩
def z : QInput := ⟨829675516563886211, 1000000000000000000⟩
def a0 : QInput := ⟨10031, 320000⟩
def a1 : QInput := ⟨2003, 32000⟩
def x0 : QInput := ⟨85921262203607322891619886036512203103, 170141183460469231731687303715884105728⟩
def x1 : QInput := ⟨253319723456283, 500000000000000⟩
def upper : Int := -400000000000000000000000
theorem check_00 : vertexCheck m c p y radius z a0 x0 50 upper = true := by decide
theorem check_01 : vertexCheck m c p y radius z a0 x1 50 upper = true := by decide
theorem check_10 : vertexCheck m c p y radius z a1 x0 50 upper = true := by decide
theorem check_11 : vertexCheck m c p y radius z a1 x1 50 upper = true := by decide
end Spin.Structured.DenseOccupationFixed.B256
