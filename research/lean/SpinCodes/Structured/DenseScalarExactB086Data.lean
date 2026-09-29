import SpinCodes.Structured.DenseScalarExactDefs
import SpinCodes.Structured.DenseOccupationFixedVertexDefs

namespace Spin.Structured.DenseScalarExact.B086
open Spin.Structured.DenseOccupationFixed
set_option maxHeartbeats 0
set_option maxRecDepth 100000
def m : QInput := ⟨9, 100⟩
def c : QInput := ⟨102967997603379096932443767333991244233, 340282366920938463463374607431768211456⟩
def p : QInput := ⟨2100848041461973, 2251799813685248⟩
def y : QInput := ⟨301699736299365, 562949953421312⟩
def radius : QInput := ⟨22089186050063798628827466131520919485681446126324844318346981602307, 2500000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000⟩
def z : QInput := ⟨61797752808988751801433432259, 500000000000000000000000000000⟩
def a0 : QInput := ⟨70001, 80000⟩
def a1 : QInput := ⟨150001, 160000⟩
def x0 : QInput := ⟨49374602183332977403672480068036683319333956789, 103845937170696552570609926584401920000000000000⟩
def x1 : QInput := ⟨81669610035739215919018503711223740925, 170141183460469231731687303715884105728⟩
def qn : Int := 633825300114114682071391547145
def qd : Int := 1267650600228229401496703205376
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
end Spin.Structured.DenseScalarExact.B086
