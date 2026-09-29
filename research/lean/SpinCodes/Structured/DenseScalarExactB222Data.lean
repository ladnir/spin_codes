import SpinCodes.Structured.DenseScalarExactDefs
import SpinCodes.Structured.DenseOccupationFixedVertexDefs

namespace Spin.Structured.DenseScalarExact.B222
open Spin.Structured.DenseOccupationFixed
set_option maxHeartbeats 0
set_option maxRecDepth 100000
def m : QInput := ⟨-1, 20⟩
def c : QInput := ⟨632748353870289731134156480146939348259, 1701411834604692317316873037158841057280⟩
def p : QInput := ⟨562379481478419, 562949953421312⟩
def y : QInput := ⟨2254084016180229, 4503599627370496⟩
def radius : QInput := ⟨11044593025031868326506898504241446697597102514313250864185032873189, 1250000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000⟩
def z : QInput := ⟨61797752808988751801433432259, 500000000000000000000000000000⟩
def a0 : QInput := ⟨5110001, 5120000⟩
def a1 : QInput := ⟨10230001, 10240000⟩
def x0 : QInput := ⟨43385881287743932506819717224474801189, 85070591730234615865843651857942052864⟩
def x1 : QInput := ⟨87621922880410294620235849331902475503, 170141183460469231731687303715884105728⟩
def qn : Int := 1267650600228229408418077977951
def qd : Int := 2535301200456458802993406410752
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
end Spin.Structured.DenseScalarExact.B222
