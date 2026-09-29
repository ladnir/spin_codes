import SpinCodes.Structured.DenseOccupationFixedVertexDefs

namespace Spin.Structured.DenseOccupationFixed.B373
set_option maxHeartbeats 0
set_option maxRecDepth 100000
def m : QInput := ⟨-6, 5⟩
def c : QInput := ⟨11167281715537733, 10000000000000000⟩
def p : QInput := ⟨5334824824586743, 144115188075855872⟩
def y : QInput := ⟨4006523395568699, 4503599627370496⟩
def radius : QInput := ⟨6277302835779813, 250000000000000000⟩
def z : QInput := ⟨186787592241933929, 200000000000000000⟩
def a0 : QInput := ⟨10063, 640000⟩
def a1 : QInput := ⟨10031, 320000⟩
def x0 : QInput := ⟨1500223291742713, 2000000000000000⟩
def x1 : QInput := ⟨1571386577205347, 2000000000000000⟩
def upper : Int := -400000000000000000000000
theorem check_00 : vertexCheck m c p y radius z a0 x0 50 upper = true := by decide
theorem check_01 : vertexCheck m c p y radius z a0 x1 50 upper = true := by decide
theorem check_10 : vertexCheck m c p y radius z a1 x0 50 upper = true := by decide
theorem check_11 : vertexCheck m c p y radius z a1 x1 50 upper = true := by decide
end Spin.Structured.DenseOccupationFixed.B373
