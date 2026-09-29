import SpinCodes.Structured.DenseOccupationFixedVertexDefs

namespace Spin.Structured.DenseOccupationFixed.B162
set_option maxHeartbeats 0
set_option maxRecDepth 100000
def m : QInput := ⟨1, 10⟩
def c : QInput := ⟨297841147503783, 1000000000000000⟩
def p : QInput := ⟨8612477418834685, 4611686018427387904⟩
def y : QInput := ⟨1604923011409571, 2251799813685248⟩
def radius : QInput := ⟨849647159067508959, 1000000000000000000⟩
def z : QInput := ⟨199465073891749811, 200000000000000000⟩
def a0 : QInput := ⟨12047, 20480000⟩
def a1 : QInput := ⟨11023, 10240000⟩
def x0 : QInput := ⟨46255924070167, 100000000000000⟩
def x1 : QInput := ⟨49374602183332977403672480068036683319333956789, 103845937170696552570609926584401920000000000000⟩
def upper : Int := -400000000000000000000000
theorem check_00 : vertexCheck m c p y radius z a0 x0 50 upper = true := by decide
theorem check_01 : vertexCheck m c p y radius z a0 x1 50 upper = true := by decide
theorem check_10 : vertexCheck m c p y radius z a1 x0 50 upper = true := by decide
theorem check_11 : vertexCheck m c p y radius z a1 x1 50 upper = true := by decide
end Spin.Structured.DenseOccupationFixed.B162
