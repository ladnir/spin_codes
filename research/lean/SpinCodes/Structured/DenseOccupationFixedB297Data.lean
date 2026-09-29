import SpinCodes.Structured.DenseOccupationFixedVertexDefs

namespace Spin.Structured.DenseOccupationFixed.B297
set_option maxHeartbeats 0
set_option maxRecDepth 100000
def m : QInput := ⟨-1, 10⟩
def c : QInput := ⟨397841147503783, 1000000000000000⟩
def p : QInput := ⟨897319249100983, 1152921504606846976⟩
def y : QInput := ⟨921242437605643, 1125899906842624⟩
def radius : QInput := ⟨923553931666586849, 1000000000000000000⟩
def z : QInput := ⟨998717836027014471, 1000000000000000000⟩
def a0 : QInput := ⟨2819, 8192000⟩
def a1 : QInput := ⟨12047, 20480000⟩
def x0 : QInput := ⟨54471334987363575166937446516365236680666043211, 103845937170696552570609926584401920000000000000⟩
def x1 : QInput := ⟨53744075929833, 100000000000000⟩
def upper : Int := -400000000000000000000000
theorem check_00 : vertexCheck m c p y radius z a0 x0 50 upper = true := by decide
theorem check_01 : vertexCheck m c p y radius z a0 x1 50 upper = true := by decide
theorem check_10 : vertexCheck m c p y radius z a1 x0 50 upper = true := by decide
theorem check_11 : vertexCheck m c p y radius z a1 x1 50 upper = true := by decide
end Spin.Structured.DenseOccupationFixed.B297
