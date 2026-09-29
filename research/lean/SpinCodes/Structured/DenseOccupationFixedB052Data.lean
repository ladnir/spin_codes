import SpinCodes.Structured.DenseOccupationFixedVertexDefs

namespace Spin.Structured.DenseOccupationFixed.B052
set_option maxHeartbeats 0
set_option maxRecDepth 100000
def m : QInput := ⟨8, 5⟩
def c : QInput := ⟨-40647089344673, 250000000000000⟩
def p : QInput := ⟨7977919565650517, 576460752303423488⟩
def y : QInput := ⟨7122801545191761, 18014398509481984⟩
def radius : QInput := ⟨251581202835488383, 500000000000000000⟩
def z : QInput := ⟨123432549285059611, 125000000000000000⟩
def a0 : QInput := ⟨10127, 1280000⟩
def a1 : QInput := ⟨10063, 640000⟩
def x0 : QInput := ⟨2469788513329, 15000000000000⟩
def x1 : QInput := ⟨4271577070219, 25000000000000⟩
def upper : Int := -400000000000000000000000
theorem check_00 : vertexCheck m c p y radius z a0 x0 50 upper = true := by decide
theorem check_01 : vertexCheck m c p y radius z a0 x1 50 upper = true := by decide
theorem check_10 : vertexCheck m c p y radius z a1 x0 50 upper = true := by decide
theorem check_11 : vertexCheck m c p y radius z a1 x1 50 upper = true := by decide
end Spin.Structured.DenseOccupationFixed.B052
