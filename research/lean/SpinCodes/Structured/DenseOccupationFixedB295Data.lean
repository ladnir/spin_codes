import SpinCodes.Structured.DenseOccupationFixedVertexDefs

namespace Spin.Structured.DenseOccupationFixed.B295
set_option maxHeartbeats 0
set_option maxRecDepth 100000
def m : QInput := ⟨-9, 100⟩
def c : QInput := ⟨3339835265656588966103687050071259581601, 8507059173023461586584365185794205286400⟩
def p : QInput := ⟨6029098534259095, 9223372036854775808⟩
def y : QInput := ⟨478059630163223, 562949953421312⟩
def radius : QInput := ⟨932845770814462529, 1000000000000000000⟩
def z : QInput := ⟨998881861735586129, 1000000000000000000⟩
def a0 : QInput := ⟨18191, 81920000⟩
def a1 : QInput := ⟨2819, 8192000⟩
def x0 : QInput := ⟨88471573424730015812668800004660364803, 170141183460469231731687303715884105728⟩
def x1 : QInput := ⟨54471334987363575166937446516365236680666043211, 103845937170696552570609926584401920000000000000⟩
def upper : Int := -400000000000000000000000
theorem check_00 : vertexCheck m c p y radius z a0 x0 50 upper = true := by decide
theorem check_01 : vertexCheck m c p y radius z a0 x1 50 upper = true := by decide
theorem check_10 : vertexCheck m c p y radius z a1 x0 50 upper = true := by decide
theorem check_11 : vertexCheck m c p y radius z a1 x1 50 upper = true := by decide
end Spin.Structured.DenseOccupationFixed.B295
