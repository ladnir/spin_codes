import SpinCodes.Structured.DenseScalarExactDefs
import SpinCodes.Structured.DenseOccupationFixedVertexDefs

namespace Spin.Structured.DenseScalarExact.B174
open Spin.Structured.DenseOccupationFixed
set_option maxHeartbeats 0
set_option maxRecDepth 100000
def m : QInput := ⟨3, 100⟩
def c : QInput := ⟨112870329263410277722384480428476839213, 340282366920938463463374607431768211456⟩
def p : QInput := ⟨8989407783734301, 9007199254740992⟩
def y : QInput := ⟨4512512968952473, 9007199254740992⟩
def radius : QInput := ⟨17671348840050987795066258658698062349895792153419758432881973828773, 2000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000⟩
def z : QInput := ⟨61797752808988751801433432259, 500000000000000000000000000000⟩
def a0 : QInput := ⟨2550001, 2560000⟩
def a1 : QInput := ⟨5110001, 5120000⟩
def x0 : QInput := ⟨246680276543717, 500000000000000⟩
def x1 : QInput := ⟨84219921256861908840067417679371902625, 170141183460469231731687303715884105728⟩
def qn : Int := 40564819207303340928912128876373
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
end Spin.Structured.DenseScalarExact.B174
