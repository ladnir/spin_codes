import SpinCodes.Structured.DenseOccupationFixedVertexDefs

namespace Spin.Structured.DenseFourierExact.B004
open Spin.Structured.DenseOccupationFixed
set_option maxHeartbeats 0
set_option maxRecDepth 100000
def m : QInput := ⟨833, 500⟩
def c : QInput := ⟨-43311, 250000⟩
def p : QInput := ⟨6438296908235089, 18014398509481984⟩
def y : QInput := ⟨6705157923346021, 18014398509481984⟩
def radius : QInput := ⟨679898891845300386972719415329788187085138769658264011560915459376851474965136363565170602553, 10000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000⟩
def z : QInput := ⟨181295472596993133903124642981, 250000000000000000000000000000⟩
def a0 : QInput := ⟨10003, 40000⟩
def a1 : QInput := ⟨6001, 16000⟩
def x0 : QInput := ⟨13, 125⟩
def x1 : QInput := ⟨18296026121, 120000000000⟩
def upper : Int := -400000000000000000000000
theorem check_00 : vertexCheck m c p y radius z a0 x0 50 upper=true := by decide
theorem check_01 : vertexCheck m c p y radius z a0 x1 50 upper=true := by decide
theorem check_10 : vertexCheck m c p y radius z a1 x0 50 upper=true := by decide
theorem check_11 : vertexCheck m c p y radius z a1 x1 50 upper=true := by decide
end Spin.Structured.DenseFourierExact.B004
