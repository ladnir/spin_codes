import SpinCodes.Structured.DenseOccupationFixedVertexDefs

namespace Spin.Structured.DenseFourierExact.B155
open Spin.Structured.DenseOccupationFixed
set_option maxHeartbeats 0
set_option maxRecDepth 100000
def m : QInput := ⟨1, 10⟩
def c : QInput := ⟨297841147503783, 1000000000000000⟩
def p : QInput := ⟨2250718002527641, 2251799813685248⟩
def y : QInput := ⟨4326852849591391, 9007199254740992⟩
def radius : QInput := ⟨40532412814737869077750591183413660747546602719725177978649935623573, 1250000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000⟩
def z : QInput := ⟨134344346016041137816724063577, 1000000000000000000000000000000⟩
def a0 : QInput := ⟨10230001, 10240000⟩
def a1 : QInput := ⟨20470001, 20480000⟩
def x0 : QInput := ⟨94384754313803, 200000000000000⟩
def x1 : QInput := ⟨49374602183332977403672480068036683319333956789, 103845937170696552570609926584401920000000000000⟩
def upper : Int := -400000000000000000000000
theorem check_00 : vertexCheck m c p y radius z a0 x0 50 upper=true := by decide
theorem check_01 : vertexCheck m c p y radius z a0 x1 50 upper=true := by decide
theorem check_10 : vertexCheck m c p y radius z a1 x0 50 upper=true := by decide
theorem check_11 : vertexCheck m c p y radius z a1 x1 50 upper=true := by decide
end Spin.Structured.DenseFourierExact.B155
