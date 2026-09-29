import SpinCodes.Structured.DenseOccupationFixedVertexDefs

namespace Spin.Structured.DenseOccupationFixed.B368
set_option maxHeartbeats 0
set_option maxRecDepth 100000
def m : QInput := ⟨-6, 5⟩
def c : QInput := ⟨11167281715537733, 10000000000000000⟩
def p : QInput := ⟨7348299824085129, 4611686018427387904⟩
def y : QInput := ⟨8189711292995049, 9007199254740992⟩
def radius : QInput := ⟨167583227283313407, 200000000000000000⟩
def z : QInput := ⟨498544644658212743, 500000000000000000⟩
def a0 : QInput := ⟨12047, 20480000⟩
def a1 : QInput := ⟨11023, 10240000⟩
def x0 : QInput := ⟨1500223291742713, 2000000000000000⟩
def x1 : QInput := ⟨1571386577205347, 2000000000000000⟩
def upper : Int := -400000000000000000000000
theorem check_00 : vertexCheck m c p y radius z a0 x0 50 upper = true := by decide
theorem check_01 : vertexCheck m c p y radius z a0 x1 50 upper = true := by decide
theorem check_10 : vertexCheck m c p y radius z a1 x0 50 upper = true := by decide
theorem check_11 : vertexCheck m c p y radius z a1 x1 50 upper = true := by decide
end Spin.Structured.DenseOccupationFixed.B368
