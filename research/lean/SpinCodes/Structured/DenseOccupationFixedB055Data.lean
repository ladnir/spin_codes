import SpinCodes.Structured.DenseOccupationFixedVertexDefs

namespace Spin.Structured.DenseOccupationFixed.B055
set_option maxHeartbeats 0
set_option maxRecDepth 100000
def m : QInput := ⟨8, 5⟩
def c : QInput := ⟨-40647089344673, 250000000000000⟩
def p : QInput := ⟨6216294641663195, 72057594037927936⟩
def y : QInput := ⟨7450700983996881, 18014398509481984⟩
def radius : QInput := ⟨2635486807327633, 250000000000000000⟩
def z : QInput := ⟨91435826982901871, 100000000000000000⟩
def a0 : QInput := ⟨2003, 32000⟩
def a1 : QInput := ⟨30029, 320000⟩
def x0 : QInput := ⟨2469788513329, 15000000000000⟩
def x1 : QInput := ⟨4271577070219, 25000000000000⟩
def upper : Int := -400000000000000000000000
theorem check_00 : vertexCheck m c p y radius z a0 x0 50 upper = true := by decide
theorem check_01 : vertexCheck m c p y radius z a0 x1 50 upper = true := by decide
theorem check_10 : vertexCheck m c p y radius z a1 x0 50 upper = true := by decide
theorem check_11 : vertexCheck m c p y radius z a1 x1 50 upper = true := by decide
end Spin.Structured.DenseOccupationFixed.B055
