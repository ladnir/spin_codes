import SpinCodes.Structured.DenseOccupationFixedVertexDefs

namespace Spin.Structured.DenseOccupationFixed.B322
set_option maxHeartbeats 0
set_option maxRecDepth 100000
def m : QInput := ⟨-1, 5⟩
def c : QInput := ⟨28224076464601, 62500000000000⟩
def p : QInput := ⟨1823726364909603, 2305843009213693952⟩
def y : QInput := ⟨929614684570483, 1125899906842624⟩
def radius : QInput := ⟨184345546323889217, 200000000000000000⟩
def z : QInput := ⟨199737063717954677, 200000000000000000⟩
def a0 : QInput := ⟨2819, 8192000⟩
def a1 : QInput := ⟨12047, 20480000⟩
def x0 : QInput := ⟨53744075929833, 100000000000000⟩
def x1 : QInput := ⟨14359431847653, 25000000000000⟩
def upper : Int := -400000000000000000000000
theorem check_00 : vertexCheck m c p y radius z a0 x0 50 upper = true := by decide
theorem check_01 : vertexCheck m c p y radius z a0 x1 50 upper = true := by decide
theorem check_10 : vertexCheck m c p y radius z a1 x0 50 upper = true := by decide
theorem check_11 : vertexCheck m c p y radius z a1 x1 50 upper = true := by decide
end Spin.Structured.DenseOccupationFixed.B322
