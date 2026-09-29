import SpinCodes.Structured.DenseScalarExactDefs
import SpinCodes.Structured.DenseOccupationFixedVertexDefs

namespace Spin.Structured.DenseScalarExact.B005
open Spin.Structured.DenseOccupationFixed
set_option maxHeartbeats 0
set_option maxRecDepth 100000
def m : QInput := ⟨3, 100⟩
def c : QInput := ⟨112870329263410277722384480428476839213, 340282366920938463463374607431768211456⟩
def p : QInput := ⟨4503060766676629, 4503599627370496⟩
def y : QInput := ⟨4504138552547359, 9007199254740992⟩
def radius : QInput := ⟨44178372100127576592267433846857777984674494245456627234899847767567, 5000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000⟩
def z : QInput := ⟨61797752808988751801433432259, 500000000000000000000000000000⟩
def a0 : QInput := ⟨163810003, 163840000⟩
def a1 : QInput := ⟨81910001, 81920000⟩
def x0 : QInput := ⟨41684710442490683359023934633467251675, 85070591730234615865843651857942052864⟩
def x1 : QInput := ⟨246680276543717, 500000000000000⟩
def qn : Int := 20282409603651672434488560972811
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
end Spin.Structured.DenseScalarExact.B005
