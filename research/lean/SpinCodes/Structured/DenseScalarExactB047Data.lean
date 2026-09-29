import SpinCodes.Structured.DenseScalarExactDefs
import SpinCodes.Structured.DenseOccupationFixedVertexDefs

namespace Spin.Structured.DenseScalarExact.B047
open Spin.Structured.DenseOccupationFixed
set_option maxHeartbeats 0
set_option maxRecDepth 100000
def m : QInput := ⟨-3, 100⟩
def c : QInput := ⟨3076970006775960790657142966285747138917, 8507059173023461586584365185794205286400⟩
def p : QInput := ⟨4503460818407165, 4503599627370496⟩
def y : QInput := ⟨2251869220306149, 4503599627370496⟩
def radius : QInput := ⟨88356744200255713033454392730599798818834657907069507907659323886801, 10000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000⟩
def z : QInput := ⟨61797752808988751801433432259, 500000000000000000000000000000⟩
def a0 : QInput := ⟨655330003, 655360000⟩
def a1 : QInput := ⟨327670001, 327680000⟩
def x0 : QInput := ⟨85921262203607322891619886036512203103, 170141183460469231731687303715884105728⟩
def x1 : QInput := ⟨253319723456283, 500000000000000⟩
def qn : Int := 10141204801825834317055835157585
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
end Spin.Structured.DenseScalarExact.B047
