import SpinCodes.Structured.DenseOccupationFixedVertexDefs

namespace Spin.Structured.DenseOccupationFixed.B218
set_option maxHeartbeats 0
set_option maxRecDepth 100000
def m : QInput := ⟨0, 1⟩
def c : QInput := ⟨3465735902799727, 10000000000000000⟩
def p : QInput := ⟨6029098534259095, 9223372036854775808⟩
def y : QInput := ⟨478059630163223, 562949953421312⟩
def radius : QInput := ⟨932845770814462529, 1000000000000000000⟩
def z : QInput := ⟨998881861735586129, 1000000000000000000⟩
def a0 : QInput := ⟨18191, 81920000⟩
def a1 : QInput := ⟨2819, 8192000⟩
def x0 : QInput := ⟨129223289418938324344194200630890939964085310021, 259614842926741381426524816461004800000000000000⟩
def x1 : QInput := ⟨130391553507803057082330615830113860035914689979, 259614842926741381426524816461004800000000000000⟩
def upper : Int := -400000000000000000000000
theorem check_00 : vertexCheck m c p y radius z a0 x0 50 upper = true := by decide
theorem check_01 : vertexCheck m c p y radius z a0 x1 50 upper = true := by decide
theorem check_10 : vertexCheck m c p y radius z a1 x0 50 upper = true := by decide
theorem check_11 : vertexCheck m c p y radius z a1 x1 50 upper = true := by decide
end Spin.Structured.DenseOccupationFixed.B218
