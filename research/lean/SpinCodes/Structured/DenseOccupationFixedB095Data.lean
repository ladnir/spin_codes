import SpinCodes.Structured.DenseOccupationFixedVertexDefs

namespace Spin.Structured.DenseOccupationFixed.B095
set_option maxHeartbeats 0
set_option maxRecDepth 100000
def m : QInput := ⟨6, 5⟩
def c : QInput := ⟨-832718284462267, 10000000000000000⟩
def p : QInput := ⟨5150272555664811, 72057594037927936⟩
def y : QInput := ⟨2289351432152455, 4503599627370496⟩
def radius : QInput := ⟨4841740775093529, 500000000000000000⟩
def z : QInput := ⟨91278442712743337, 100000000000000000⟩
def a0 : QInput := ⟨30061, 640000⟩
def a1 : QInput := ⟨2003, 32000⟩
def x0 : QInput := ⟨428613422794653, 2000000000000000⟩
def x1 : QInput := ⟨499776708257287, 2000000000000000⟩
def upper : Int := -400000000000000000000000
theorem check_00 : vertexCheck m c p y radius z a0 x0 50 upper = true := by decide
theorem check_01 : vertexCheck m c p y radius z a0 x1 50 upper = true := by decide
theorem check_10 : vertexCheck m c p y radius z a1 x0 50 upper = true := by decide
theorem check_11 : vertexCheck m c p y radius z a1 x1 50 upper = true := by decide
end Spin.Structured.DenseOccupationFixed.B095
