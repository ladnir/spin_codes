import SpinCodes.Structured.DenseOccupationFixedVertexDefs

namespace Spin.Structured.DenseFourierExact.B293
open Spin.Structured.DenseOccupationFixed
set_option maxHeartbeats 0
set_option maxRecDepth 100000
def m : QInput := ⟨-33, 20⟩
def c : QInput := ⟨46226687491353, 31250000000000⟩
def p : QInput := ⟨3453082862849675, 4503599627370496⟩
def y : QInput := ⟨6629750246321639, 9007199254740992⟩
def radius : QInput := ⟨6946212725587510547174942130333153149963649028366927936839931643271, 156250000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000⟩
def z : QInput := ⟨65594999908555291432260744129, 500000000000000000000000000000⟩
def a0 : QInput := ⟨30001, 40000⟩
def a1 : QInput := ⟨70001, 80000⟩
def x0 : QInput := ⟨8390967064323, 10000000000000⟩
def x1 : QInput := ⟨4208400791377, 5000000000000⟩
def upper : Int := -400000000000000000000000
theorem check_00 : vertexCheck m c p y radius z a0 x0 50 upper=true := by decide
theorem check_01 : vertexCheck m c p y radius z a0 x1 50 upper=true := by decide
theorem check_10 : vertexCheck m c p y radius z a1 x0 50 upper=true := by decide
theorem check_11 : vertexCheck m c p y radius z a1 x1 50 upper=true := by decide
end Spin.Structured.DenseFourierExact.B293
