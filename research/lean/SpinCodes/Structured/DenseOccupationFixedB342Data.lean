import SpinCodes.Structured.DenseOccupationFixedVertexDefs

namespace Spin.Structured.DenseOccupationFixed.B342
set_option maxHeartbeats 0
set_option maxRecDepth 100000
def m : QInput := ⟨-3, 5⟩
def c : QInput := ⟨6909341020092481, 10000000000000000⟩
def p : QInput := ⟨5930898936995693, 2305843009213693952⟩
def y : QInput := ⟨3804329281385513, 4503599627370496⟩
def radius : QInput := ⟨380281973519517571, 500000000000000000⟩
def z : QInput := ⟨995255006651210737, 1000000000000000000⟩
def a0 : QInput := ⟨11023, 10240000⟩
def a1 : QInput := ⟨10511, 5120000⟩
def x0 : QInput := ⟨1244734237944081, 2000000000000000⟩
def x1 : QInput := ⟨1336129721374967, 2000000000000000⟩
def upper : Int := -400000000000000000000000
theorem check_00 : vertexCheck m c p y radius z a0 x0 50 upper = true := by decide
theorem check_01 : vertexCheck m c p y radius z a0 x1 50 upper = true := by decide
theorem check_10 : vertexCheck m c p y radius z a1 x0 50 upper = true := by decide
theorem check_11 : vertexCheck m c p y radius z a1 x1 50 upper = true := by decide
end Spin.Structured.DenseOccupationFixed.B342
