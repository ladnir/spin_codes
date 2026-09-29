import SpinCodes.Structured.DenseOccupationFixedVertexDefs

namespace Spin.Structured.DenseFourierExact.B045
open Spin.Structured.DenseOccupationFixed
set_option maxHeartbeats 0
set_option maxRecDepth 100000
def m : QInput := ⟨31, 20⟩
def c : QInput := ⟨-77022601619127, 500000000000000⟩
def p : QInput := ⟨7761352683183377, 9007199254740992⟩
def y : QInput := ⟨8860320953993497, 18014398509481984⟩
def radius : QInput := ⟨39837107984645418213980371460804149715207961951433463607902953727300717, 10000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000⟩
def z : QInput := ⟨20279238649322490279482906123, 125000000000000000000000000000⟩
def a0 : QInput := ⟨30001, 40000⟩
def a1 : QInput := ⟨70001, 80000⟩
def x0 : QInput := ⟨4271577070219, 25000000000000⟩
def x1 : QInput := ⟨13956016256281, 75000000000000⟩
def upper : Int := -400000000000000000000000
theorem check_00 : vertexCheck m c p y radius z a0 x0 50 upper=true := by decide
theorem check_01 : vertexCheck m c p y radius z a0 x1 50 upper=true := by decide
theorem check_10 : vertexCheck m c p y radius z a1 x0 50 upper=true := by decide
theorem check_11 : vertexCheck m c p y radius z a1 x1 50 upper=true := by decide
end Spin.Structured.DenseFourierExact.B045
