import SpinCodes.Structured.DenseScalarExactDefs
import SpinCodes.Structured.DenseOccupationFixedVertexDefs

namespace Spin.Structured.DenseScalarExact.B156
open Spin.Structured.DenseOccupationFixed
set_option maxHeartbeats 0
set_option maxRecDepth 100000
def m : QInput := ⟨-1, 100⟩
def c : QInput := ⟨1495524372286176733882761540124617467907, 4253529586511730793292182592897102643200⟩
def p : QInput := ⟨4419997189385557, 4503599627370496⟩
def y : QInput := ⟨4588783371256219, 9007199254740992⟩
def radius : QInput := ⟨88356744200255584758768073507337993909500774591238135222099097685133, 10000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000⟩
def z : QInput := ⟨61797752808988751801433432259, 500000000000000000000000000000⟩
def a0 : QInput := ⟨310001, 320000⟩
def a1 : QInput := ⟨630001, 640000⟩
def x0 : QInput := ⟨130391553507803057082330615830113860035914689979, 259614842926741381426524816461004800000000000000⟩
def x1 : QInput := ⟨85921262203607322891619886036512203103, 170141183460469231731687303715884105728⟩
def qn : Int := 20282409603651668929039825028983
def qd : Int := 40564819207303340847894502572032
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
end Spin.Structured.DenseScalarExact.B156
