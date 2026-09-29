import SpinCodes.Structured.DenseOccupationFixedVertexDefs

namespace Spin.Structured.DenseFourierExact.B304
open Spin.Structured.DenseOccupationFixed
set_option maxHeartbeats 0
set_option maxRecDepth 100000
def m : QInput := ⟨-833, 500⟩
def c : QInput := ⟨373189, 250000⟩
def p : QInput := ⟨7898872143836449, 18014398509481984⟩
def y : QInput := ⟨8474020851551243, 9007199254740992⟩
def radius : QInput := ⟨316110245518880297727394223217063501059138860419352387311649251587907761, 5000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000⟩
def z : QInput := ⟨180483675024133742363734526109, 1000000000000000000000000000000⟩
def a0 : QInput := ⟨10003, 40000⟩
def a1 : QInput := ⟨10001, 20000⟩
def x0 : QInput := ⟨101703973879, 120000000000⟩
def x1 : QInput := ⟨112, 125⟩
def upper : Int := -400000000000000000000000
theorem check_00 : vertexCheck m c p y radius z a0 x0 50 upper=true := by decide
theorem check_01 : vertexCheck m c p y radius z a0 x1 50 upper=true := by decide
theorem check_10 : vertexCheck m c p y radius z a1 x0 50 upper=true := by decide
theorem check_11 : vertexCheck m c p y radius z a1 x1 50 upper=true := by decide
end Spin.Structured.DenseFourierExact.B304
