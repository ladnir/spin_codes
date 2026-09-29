import SpinCodes.Structured.DenseOccupationFixedVertexDefs

namespace Spin.Structured.DenseFourierExact.B060
open Spin.Structured.DenseOccupationFixed
set_option maxHeartbeats 0
set_option maxRecDepth 100000
def m : QInput := ⟨6, 5⟩
def c : QInput := ⟨-832718284462267, 10000000000000000⟩
def p : QInput := ⟨3988264048149027, 9007199254740992⟩
def y : QInput := ⟨7396803716810457, 18014398509481984⟩
def radius : QInput := ⟨79108463780419843474990620677374414593516202329253383098334771254172228099256928468194007, 1000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000⟩
def z : QInput := ⟨577880326639773235549538724383, 1000000000000000000000000000000⟩
def a0 : QInput := ⟨6001, 16000⟩
def a1 : QInput := ⟨70009, 160000⟩
def x0 : QInput := ⟨428613422794653, 2000000000000000⟩
def x1 : QInput := ⟨499776708257287, 2000000000000000⟩
def upper : Int := -400000000000000000000000
theorem check_00 : vertexCheck m c p y radius z a0 x0 50 upper=true := by decide
theorem check_01 : vertexCheck m c p y radius z a0 x1 50 upper=true := by decide
theorem check_10 : vertexCheck m c p y radius z a1 x0 50 upper=true := by decide
theorem check_11 : vertexCheck m c p y radius z a1 x1 50 upper=true := by decide
end Spin.Structured.DenseFourierExact.B060
