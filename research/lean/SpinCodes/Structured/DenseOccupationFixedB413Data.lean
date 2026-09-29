import SpinCodes.Structured.DenseOccupationFixedVertexDefs

namespace Spin.Structured.DenseOccupationFixed.B413
set_option maxHeartbeats 0
set_option maxRecDepth 100000
def m : QInput := ⟨-33, 20⟩
def c : QInput := ⟨46226687491353, 31250000000000⟩
def p : QInput := ⟨6412853813896097, 1152921504606846976⟩
def y : QInput := ⟨520435057616331, 562949953421312⟩
def radius : QInput := ⟨275827164287343331, 500000000000000000⟩
def z : QInput := ⟨989687944403566779, 1000000000000000000⟩
def a0 : QInput := ⟨10511, 5120000⟩
def a1 : QInput := ⟨2051, 512000⟩
def x0 : QInput := ⟨8390967064323, 10000000000000⟩
def x1 : QInput := ⟨4208400791377, 5000000000000⟩
def upper : Int := -400000000000000000000000
theorem check_00 : vertexCheck m c p y radius z a0 x0 50 upper = true := by decide
theorem check_01 : vertexCheck m c p y radius z a0 x1 50 upper = true := by decide
theorem check_10 : vertexCheck m c p y radius z a1 x0 50 upper = true := by decide
theorem check_11 : vertexCheck m c p y radius z a1 x1 50 upper = true := by decide
end Spin.Structured.DenseOccupationFixed.B413
