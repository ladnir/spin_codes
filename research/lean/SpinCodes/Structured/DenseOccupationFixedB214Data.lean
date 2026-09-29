import SpinCodes.Structured.DenseOccupationFixedVertexDefs

namespace Spin.Structured.DenseOccupationFixed.B214
set_option maxHeartbeats 0
set_option maxRecDepth 100000
def m : QInput := ⟨1, 100⟩
def c : QInput := ⟨58119563056842377037993588567825857659, 170141183460469231731687303715884105728⟩
def p : QInput := ⟨4831935820581855, 18446744073709551616⟩
def y : QInput := ⟨7288151109622443, 9007199254740992⟩
def radius : QInput := ⟨38935291265312167, 40000000000000000⟩
def z : QInput := ⟨999567439056017201, 1000000000000000000⟩
def a0 : QInput := ⟨1, 10000⟩
def a1 : QInput := ⟨18191, 81920000⟩
def x0 : QInput := ⟨84219921256861908840067417679371902625, 170141183460469231731687303715884105728⟩
def x1 : QInput := ⟨129223289418938324344194200630890939964085310021, 259614842926741381426524816461004800000000000000⟩
def upper : Int := -400000000000000000000000
theorem check_00 : vertexCheck m c p y radius z a0 x0 50 upper = true := by decide
theorem check_01 : vertexCheck m c p y radius z a0 x1 50 upper = true := by decide
theorem check_10 : vertexCheck m c p y radius z a1 x0 50 upper = true := by decide
theorem check_11 : vertexCheck m c p y radius z a1 x1 50 upper = true := by decide
end Spin.Structured.DenseOccupationFixed.B214
