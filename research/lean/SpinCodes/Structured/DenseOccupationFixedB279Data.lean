import SpinCodes.Structured.DenseOccupationFixedVertexDefs

namespace Spin.Structured.DenseOccupationFixed.B279
set_option maxHeartbeats 0
set_option maxRecDepth 100000
def m : QInput := ⟨-7, 100⟩
def c : QInput := ⟨1625681846115929475145509125033299608399, 4253529586511730793292182592897102643200⟩
def p : QInput := ⟨7539775244674395, 576460752303423488⟩
def y : QInput := ⟨6970124630314481, 9007199254740992⟩
def radius : QInput := ⟨280298832600148967, 1000000000000000000⟩
def z : QInput := ⟨976189201681259753, 1000000000000000000⟩
def a0 : QInput := ⟨2051, 512000⟩
def a1 : QInput := ⟨10127, 1280000⟩
def x0 : QInput := ⟨87621922880410294620235849331902475503, 170141183460469231731687303715884105728⟩
def x1 : QInput := ⟨12967792439091, 25000000000000⟩
def upper : Int := -400000000000000000000000
theorem check_00 : vertexCheck m c p y radius z a0 x0 50 upper = true := by decide
theorem check_01 : vertexCheck m c p y radius z a0 x1 50 upper = true := by decide
theorem check_10 : vertexCheck m c p y radius z a1 x0 50 upper = true := by decide
theorem check_11 : vertexCheck m c p y radius z a1 x1 50 upper = true := by decide
end Spin.Structured.DenseOccupationFixed.B279
