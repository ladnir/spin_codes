import SpinCodes.Structured.DenseOccupationFixedVertexDefs

namespace Spin.Structured.DenseOccupationFixed.B308
set_option maxHeartbeats 0
set_option maxRecDepth 100000
def m : QInput := ⟨-7, 100⟩
def c : QInput := ⟨1625681846115929475145509125033299608399, 4253529586511730793292182592897102643200⟩
def p : QInput := ⟨6934123379006515, 1152921504606846976⟩
def y : QInput := ⟨3150303217345729, 4503599627370496⟩
def radius : QInput := ⟨589568824387140823, 1000000000000000000⟩
def z : QInput := ⟨990489247114675271, 1000000000000000000⟩
def a0 : QInput := ⟨10511, 5120000⟩
def a1 : QInput := ⟨2051, 512000⟩
def x0 : QInput := ⟨12967792439091, 25000000000000⟩
def x1 : QInput := ⟨88471573424730015812668800004660364803, 170141183460469231731687303715884105728⟩
def upper : Int := -400000000000000000000000
theorem check_00 : vertexCheck m c p y radius z a0 x0 50 upper = true := by decide
theorem check_01 : vertexCheck m c p y radius z a0 x1 50 upper = true := by decide
theorem check_10 : vertexCheck m c p y radius z a1 x0 50 upper = true := by decide
theorem check_11 : vertexCheck m c p y radius z a1 x1 50 upper = true := by decide
end Spin.Structured.DenseOccupationFixed.B308
