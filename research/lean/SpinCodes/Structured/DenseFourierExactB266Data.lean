import SpinCodes.Structured.DenseOccupationFixedVertexDefs

namespace Spin.Structured.DenseFourierExact.B266
open Spin.Structured.DenseOccupationFixed
set_option maxHeartbeats 0
set_option maxRecDepth 100000
def m : QInput := ⟨-6, 5⟩
def c : QInput := ⟨11167281715537733, 10000000000000000⟩
def p : QInput := ⟨2204388020114543, 2251799813685248⟩
def y : QInput := ⟨2608189365246017, 4503599627370496⟩
def radius : QInput := ⟨58058119024108778004704772379665260268068448606547372694823392420157, 1000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000⟩
def z : QInput := ⟨134205669425665989180916433867, 1000000000000000000000000000000⟩
def a0 : QInput := ⟨150001, 160000⟩
def a1 : QInput := ⟨1, 1⟩
def x0 : QInput := ⟨1500223291742713, 2000000000000000⟩
def x1 : QInput := ⟨1571386577205347, 2000000000000000⟩
def upper : Int := -400000000000000000000000
theorem check_00 : vertexCheck m c p y radius z a0 x0 50 upper=true := by decide
theorem check_01 : vertexCheck m c p y radius z a0 x1 50 upper=true := by decide
theorem check_10 : vertexCheck m c p y radius z a1 x0 50 upper=true := by decide
theorem check_11 : vertexCheck m c p y radius z a1 x1 50 upper=true := by decide
end Spin.Structured.DenseFourierExact.B266
