import SpinCodes.Structured.DenseScalarExactDefs
import SpinCodes.Structured.DenseOccupationFixedVertexDefs

namespace Spin.Structured.DenseScalarExact.B062
open Spin.Structured.DenseOccupationFixed
set_option maxHeartbeats 0
set_option maxRecDepth 100000
def m : QInput := ⟨-7, 100⟩
def c : QInput := ⟨1625681846115929475145509125033299608399, 4253529586511730793292182592897102643200⟩
def p : QInput := ⟨4503314109998881, 4503599627370496⟩
def y : QInput := ⟨1125971290711093, 2251799813685248⟩
def radius : QInput := ⟨44178372100127926088983268575232415626391275084199492055213573981181, 5000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000⟩
def z : QInput := ⟨61797752808988751801433432259, 500000000000000000000000000000⟩
def a0 : QInput := ⟨81910001, 81920000⟩
def a1 : QInput := ⟨163830001, 163840000⟩
def x0 : QInput := ⟨87621922880410294620235849331902475503, 170141183460469231731687303715884105728⟩
def x1 : QInput := ⟨12967792439091, 25000000000000⟩
def qn : Int := 5070602400912917078547924286933
def qd : Int := 10141204801825835211973625643008
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
end Spin.Structured.DenseScalarExact.B062
