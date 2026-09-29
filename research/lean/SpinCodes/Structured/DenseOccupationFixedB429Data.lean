import SpinCodes.Structured.DenseOccupationFixedVertexDefs

namespace Spin.Structured.DenseOccupationFixed.B429
set_option maxHeartbeats 0
set_option maxRecDepth 100000
def m : QInput := ⟨-833, 500⟩
def c : QInput := ⟨373189, 250000⟩
def p : QInput := ⟨830979024820423, 144115188075855872⟩
def y : QInput := ⟨8442289394766573, 9007199254740992⟩
def radius : QInput := ⟨535893138787010667, 1000000000000000000⟩
def z : QInput := ⟨247290628373232957, 250000000000000000⟩
def a0 : QInput := ⟨10511, 5120000⟩
def a1 : QInput := ⟨2051, 512000⟩
def x0 : QInput := ⟨101703973879, 120000000000⟩
def x1 : QInput := ⟨112, 125⟩
def upper : Int := -400000000000000000000000
theorem check_00 : vertexCheck m c p y radius z a0 x0 50 upper = true := by decide
theorem check_01 : vertexCheck m c p y radius z a0 x1 50 upper = true := by decide
theorem check_10 : vertexCheck m c p y radius z a1 x0 50 upper = true := by decide
theorem check_11 : vertexCheck m c p y radius z a1 x1 50 upper = true := by decide
end Spin.Structured.DenseOccupationFixed.B429
