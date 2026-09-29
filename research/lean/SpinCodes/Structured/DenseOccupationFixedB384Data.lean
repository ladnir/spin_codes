import SpinCodes.Structured.DenseOccupationFixedVertexDefs

namespace Spin.Structured.DenseOccupationFixed.B384
set_option maxHeartbeats 0
set_option maxRecDepth 100000
def m : QInput := ⟨-31, 20⟩
def c : QInput := ⟨697977398380873, 500000000000000⟩
def p : QInput := ⟨6718804664920955, 18446744073709551616⟩
def y : QInput := ⟨4305814936677867, 4503599627370496⟩
def radius : QInput := ⟨956707446014747587, 1000000000000000000⟩
def z : QInput := ⟨249820978036117291, 250000000000000000⟩
def a0 : QInput := ⟨1, 10000⟩
def a1 : QInput := ⟨2819, 8192000⟩
def x0 : QInput := ⟨61043983743719, 75000000000000⟩
def x1 : QInput := ⟨20728422929781, 25000000000000⟩
def upper : Int := -400000000000000000000000
theorem check_00 : vertexCheck m c p y radius z a0 x0 50 upper = true := by decide
theorem check_01 : vertexCheck m c p y radius z a0 x1 50 upper = true := by decide
theorem check_10 : vertexCheck m c p y radius z a1 x0 50 upper = true := by decide
theorem check_11 : vertexCheck m c p y radius z a1 x1 50 upper = true := by decide
end Spin.Structured.DenseOccupationFixed.B384
