import SpinCodes.Structured.DenseOccupationFixedVertexDefs

namespace Spin.Structured.DenseFourierExact.B210
open Spin.Structured.DenseOccupationFixed
set_option maxHeartbeats 0
set_option maxRecDepth 100000
def m : QInput := ⟨-1, 10⟩
def c : QInput := ⟨397841147503783, 1000000000000000⟩
def p : QInput := ⟨4557672428423, 8796093022208⟩
def y : QInput := ⟨3722964314139027, 4503599627370496⟩
def radius : QInput := ⟨12358993241762213709760320435033430098140351095429819626195655081395551, 10000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000⟩
def z : QInput := ⟨2406450320492416543689332129, 15625000000000000000000000000⟩
def a0 : QInput := ⟨10003, 40000⟩
def a1 : QInput := ⟨10001, 20000⟩
def x0 : QInput := ⟨54471334987363575166937446516365236680666043211, 103845937170696552570609926584401920000000000000⟩
def x1 : QInput := ⟨53744075929833, 100000000000000⟩
def upper : Int := -400000000000000000000000
theorem check_00 : vertexCheck m c p y radius z a0 x0 50 upper=true := by decide
theorem check_01 : vertexCheck m c p y radius z a0 x1 50 upper=true := by decide
theorem check_10 : vertexCheck m c p y radius z a1 x0 50 upper=true := by decide
theorem check_11 : vertexCheck m c p y radius z a1 x1 50 upper=true := by decide
end Spin.Structured.DenseFourierExact.B210
