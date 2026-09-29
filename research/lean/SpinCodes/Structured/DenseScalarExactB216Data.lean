import SpinCodes.Structured.DenseScalarExactDefs
import SpinCodes.Structured.DenseOccupationFixedVertexDefs

namespace Spin.Structured.DenseScalarExact.B216
open Spin.Structured.DenseOccupationFixed
set_option maxHeartbeats 0
set_option maxRecDepth 100000
def m : QInput := ⟨-1, 20⟩
def c : QInput := ⟨632748353870289731134156480146939348259, 1701411834604692317316873037158841057280⟩
def p : QInput := ⟨8970803836287323, 9007199254740992⟩
def y : QInput := ⟨4521871166463003, 9007199254740992⟩
def radius : QInput := ⟨88356744200255237775700184758228698745376641207331249730482427059851, 10000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000⟩
def z : QInput := ⟨61797752808988751801433432259, 500000000000000000000000000000⟩
def a0 : QInput := ⟨1270001, 1280000⟩
def a1 : QInput := ⟨2550001, 2560000⟩
def x0 : QInput := ⟨43385881287743932506819717224474801189, 85070591730234615865843651857942052864⟩
def x1 : QInput := ⟨87621922880410294620235849331902475503, 170141183460469231731687303715884105728⟩
def qn : Int := 40564819207303339453641157410969
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
end Spin.Structured.DenseScalarExact.B216
