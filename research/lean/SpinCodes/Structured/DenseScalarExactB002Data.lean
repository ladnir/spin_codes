import SpinCodes.Structured.DenseScalarExactDefs
import SpinCodes.Structured.DenseOccupationFixedVertexDefs

namespace Spin.Structured.DenseScalarExact.B002
open Spin.Structured.DenseOccupationFixed
set_option maxHeartbeats 0
set_option maxRecDepth 100000
def m : QInput := ⟨7, 100⟩
def c : QInput := ⟨53117391002404332784602253741220096935, 170141183460469231731687303715884105728⟩
def p : QInput := ⟨9006913727979809, 9007199254740992⟩
def y : QInput := ⟨4503742395276807, 9007199254740992⟩
def radius : QInput := ⟨17671348840051012240027438102640428012619154974615448948420594120313, 2000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000⟩
def z : QInput := ⟨61797752808988751801433432259, 500000000000000000000000000000⟩
def a0 : QInput := ⟨81910001, 81920000⟩
def a1 : QInput := ⟨1, 1⟩
def x0 : QInput := ⟨81669610035739215919018503711223740925, 170141183460469231731687303715884105728⟩
def x1 : QInput := ⟨12032207560909, 25000000000000⟩
def qn : Int := 40564819207303340265603461989863
def qd : Int := 81129638414606681695789005144064
theorem check0 : check qn qd z.num z.den radius.num radius.den 0 = true := by decide
theorem check48 : check qn qd z.num z.den radius.num radius.den 48 = true := by decide
theorem check56 : check qn qd z.num z.den radius.num radius.den 56 = true := by decide
theorem check64 : check qn qd z.num z.den radius.num radius.den 64 = true := by decide
theorem check72 : check qn qd z.num z.den radius.num radius.den 72 = true := by decide
theorem check80 : check qn qd z.num z.den radius.num radius.den 80 = true := by decide
def upper : Int := -400000000000000000000000
theorem check_00 : vertexCheck m c p y radius z a0 x0 50 upper = true := by decide
theorem check_01 : vertexCheck m c p y radius z a0 x1 50 upper = true := by decide
theorem check_10 : vertexCheck m c p y radius z a1 x0 50 upper = true := by decide
theorem check_11 : vertexCheck m c p y radius z a1 x1 50 upper = true := by decide
end Spin.Structured.DenseScalarExact.B002
