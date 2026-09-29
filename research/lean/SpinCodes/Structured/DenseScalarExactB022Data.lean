import SpinCodes.Structured.DenseScalarExactDefs
import SpinCodes.Structured.DenseOccupationFixedVertexDefs

namespace Spin.Structured.DenseScalarExact.B022
open Spin.Structured.DenseOccupationFixed
set_option maxHeartbeats 0
set_option maxRecDepth 100000
def m : QInput := ⟨1, 20⟩
def c : QInput := ⟨109535552428011023053662565657799459079, 340282366920938463463374607431768211456⟩
def p : QInput := ⟨2251783158420411, 2251799813685248⟩
def y : QInput := ⟨4503632938146551, 9007199254740992⟩
def radius : QInput := ⟨88356744200255128720052476886508320971955112528247310119573280488567, 10000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000⟩
def z : QInput := ⟨61797752808988751801433432259, 500000000000000000000000000000⟩
def a0 : QInput := ⟨655350001, 655360000⟩
def a1 : QInput := ⟨1, 1⟩
def x0 : QInput := ⟨82519260580058937111451454383981630225, 170141183460469231731687303715884105728⟩
def x1 : QInput := ⟨487324427761897, 1000000000000000⟩
def qn : Int := 10141204801825836104747187652461
def qd : Int := 20282409603651670423947251286016
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
end Spin.Structured.DenseScalarExact.B022
