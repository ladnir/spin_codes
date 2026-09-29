import SpinCodes.Structured.DenseOccupationFixedVertexDefs

namespace Spin.Structured.DenseOccupationFixed.B034
set_option maxHeartbeats 0
set_option maxRecDepth 100000
def m : QInput := ⟨163, 100⟩
def c : QInput := ⟨-3350558688107, 20000000000000⟩
def p : QInput := ⟨4605167708181261, 18446744073709551616⟩
def y : QInput := ⟨8250128820775807, 18014398509481984⟩
def radius : QInput := ⟨985526468309376367, 1000000000000000000⟩
def z : QInput := ⟨999767766661628143, 1000000000000000000⟩
def a0 : QInput := ⟨1, 10000⟩
def a1 : QInput := ⟨18191, 81920000⟩
def x0 : QInput := ⟨1609032935677, 10000000000000⟩
def x1 : QInput := ⟨2469788513329, 15000000000000⟩
def upper : Int := -400000000000000000000000
theorem check_00 : vertexCheck m c p y radius z a0 x0 50 upper = true := by decide
theorem check_01 : vertexCheck m c p y radius z a0 x1 50 upper = true := by decide
theorem check_10 : vertexCheck m c p y radius z a1 x0 50 upper = true := by decide
theorem check_11 : vertexCheck m c p y radius z a1 x1 50 upper = true := by decide
end Spin.Structured.DenseOccupationFixed.B034
