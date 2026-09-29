import SpinCodes.Structured.DenseOccupationFixedVertexDefs

namespace Spin.Structured.DenseOccupationFixed.B420
set_option maxHeartbeats 0
set_option maxRecDepth 100000
def m : QInput := ⟨-83, 50⟩
def c : QInput := ⟨29753416026121, 20000000000000⟩
def p : QInput := ⟨3606530255009465, 1152921504606846976⟩
def y : QInput := ⟨4273739439765363, 4503599627370496⟩
def radius : QInput := ⟨688438337421590637, 1000000000000000000⟩
def z : QInput := ⟨496706335540625369, 500000000000000000⟩
def a0 : QInput := ⟨11023, 10240000⟩
def a1 : QInput := ⟨10511, 5120000⟩
def x0 : QInput := ⟨4208400791377, 5000000000000⟩
def x1 : QInput := ⟨101703973879, 120000000000⟩
def upper : Int := -400000000000000000000000
theorem check_00 : vertexCheck m c p y radius z a0 x0 50 upper = true := by decide
theorem check_01 : vertexCheck m c p y radius z a0 x1 50 upper = true := by decide
theorem check_10 : vertexCheck m c p y radius z a1 x0 50 upper = true := by decide
theorem check_11 : vertexCheck m c p y radius z a1 x1 50 upper = true := by decide
end Spin.Structured.DenseOccupationFixed.B420
