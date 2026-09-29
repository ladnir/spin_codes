import SpinCodes.Structured.DenseOccupationFixedVertexDefs

namespace Spin.Structured.DenseFourierExact.B272
open Spin.Structured.DenseOccupationFixed
set_option maxHeartbeats 0
set_option maxRecDepth 100000
def m : QInput := ⟨-31, 20⟩
def c : QInput := ⟨697977398380873, 500000000000000⟩
def p : QInput := ⟨85033740808161, 562949953421312⟩
def y : QInput := ⟨8334938389558539, 9007199254740992⟩
def radius : QInput := ⟨311343612414962374183656586894033144731010530126261514857487300691032579237389966807009418139, 10000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000⟩
def z : QInput := ⟨711037741313782822851429396129, 1000000000000000000000000000000⟩
def a0 : QInput := ⟨2003, 32000⟩
def a1 : QInput := ⟨10007, 80000⟩
def x0 : QInput := ⟨61043983743719, 75000000000000⟩
def x1 : QInput := ⟨20728422929781, 25000000000000⟩
def upper : Int := -400000000000000000000000
theorem check_00 : vertexCheck m c p y radius z a0 x0 50 upper=true := by decide
theorem check_01 : vertexCheck m c p y radius z a0 x1 50 upper=true := by decide
theorem check_10 : vertexCheck m c p y radius z a1 x0 50 upper=true := by decide
theorem check_11 : vertexCheck m c p y radius z a1 x1 50 upper=true := by decide
end Spin.Structured.DenseFourierExact.B272
