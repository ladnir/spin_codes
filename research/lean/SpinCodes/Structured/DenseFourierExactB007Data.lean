import SpinCodes.Structured.DenseOccupationFixedVertexDefs

namespace Spin.Structured.DenseFourierExact.B007
open Spin.Structured.DenseOccupationFixed
set_option maxHeartbeats 0
set_option maxRecDepth 100000
def m : QInput := ⟨833, 500⟩
def c : QInput := ⟨-43311, 250000⟩
def p : QInput := ⟨3927908298396545, 4503599627370496⟩
def y : QInput := ⟨8349925854794397, 18014398509481984⟩
def radius : QInput := ⟨3849529999616570012482186105593632168746379158673550315478232814068257137, 10000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000⟩
def z : QInput := ⟨95605631708104962096277475509, 500000000000000000000000000000⟩
def a0 : QInput := ⟨30001, 40000⟩
def a1 : QInput := ⟨70001, 80000⟩
def x0 : QInput := ⟨13, 125⟩
def x1 : QInput := ⟨18296026121, 120000000000⟩
def upper : Int := -400000000000000000000000
theorem check_00 : vertexCheck m c p y radius z a0 x0 50 upper=true := by decide
theorem check_01 : vertexCheck m c p y radius z a0 x1 50 upper=true := by decide
theorem check_10 : vertexCheck m c p y radius z a1 x0 50 upper=true := by decide
theorem check_11 : vertexCheck m c p y radius z a1 x1 50 upper=true := by decide
end Spin.Structured.DenseFourierExact.B007
