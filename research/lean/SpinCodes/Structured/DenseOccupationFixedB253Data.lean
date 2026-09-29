import SpinCodes.Structured.DenseOccupationFixedVertexDefs

namespace Spin.Structured.DenseOccupationFixed.B253
set_option maxHeartbeats 0
set_option maxRecDepth 100000
def m : QInput := ⟨0, 1⟩
def c : QInput := ⟨3465735902799727, 10000000000000000⟩
def p : QInput := ⟨959837913988893, 18014398509481984⟩
def y : QInput := ⟨7294169837842663, 9007199254740992⟩
def radius : QInput := ⟨3976574297585327, 1000000000000000000⟩
def z : QInput := ⟨11204571114096203, 12500000000000000⟩
def a0 : QInput := ⟨10063, 640000⟩
def a1 : QInput := ⟨10031, 320000⟩
def x0 : QInput := ⟨129223289418938324344194200630890939964085310021, 259614842926741381426524816461004800000000000000⟩
def x1 : QInput := ⟨130391553507803057082330615830113860035914689979, 259614842926741381426524816461004800000000000000⟩
def upper : Int := -400000000000000000000000
theorem check_00 : vertexCheck m c p y radius z a0 x0 50 upper = true := by decide
theorem check_01 : vertexCheck m c p y radius z a0 x1 50 upper = true := by decide
theorem check_10 : vertexCheck m c p y radius z a1 x0 50 upper = true := by decide
theorem check_11 : vertexCheck m c p y radius z a1 x1 50 upper = true := by decide
end Spin.Structured.DenseOccupationFixed.B253
