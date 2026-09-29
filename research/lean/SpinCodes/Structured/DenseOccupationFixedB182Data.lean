import SpinCodes.Structured.DenseOccupationFixedVertexDefs

namespace Spin.Structured.DenseOccupationFixed.B182
set_option maxHeartbeats 0
set_option maxRecDepth 100000
def m : QInput := ⟨1, 20⟩
def c : QInput := ⟨109535552428011023053662565657799459079, 340282366920938463463374607431768211456⟩
def p : QInput := ⟨4831935820581855, 18446744073709551616⟩
def y : QInput := ⟨7288151109622443, 9007199254740992⟩
def radius : QInput := ⟨38935291265312167, 40000000000000000⟩
def z : QInput := ⟨999567439056017201, 1000000000000000000⟩
def a0 : QInput := ⟨1, 10000⟩
def a1 : QInput := ⟨18191, 81920000⟩
def x0 : QInput := ⟨82519260580058937111451454383981630225, 170141183460469231731687303715884105728⟩
def x1 : QInput := ⟨41684710442490683359023934633467251675, 85070591730234615865843651857942052864⟩
def upper : Int := -400000000000000000000000
theorem check_00 : vertexCheck m c p y radius z a0 x0 50 upper = true := by decide
theorem check_01 : vertexCheck m c p y radius z a0 x1 50 upper = true := by decide
theorem check_10 : vertexCheck m c p y radius z a1 x0 50 upper = true := by decide
theorem check_11 : vertexCheck m c p y radius z a1 x1 50 upper = true := by decide
end Spin.Structured.DenseOccupationFixed.B182
