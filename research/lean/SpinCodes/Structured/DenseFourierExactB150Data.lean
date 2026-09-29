import SpinCodes.Structured.DenseOccupationFixedVertexDefs

namespace Spin.Structured.DenseFourierExact.B150
open Spin.Structured.DenseOccupationFixed
set_option maxHeartbeats 0
set_option maxRecDepth 100000
def m : QInput := ⟨1, 10⟩
def c : QInput := ⟨297841147503783, 1000000000000000⟩
def p : QInput := ⟨560795178958293, 562949953421312⟩
def y : QInput := ⟨8527337085287279, 18014398509481984⟩
def radius : QInput := ⟨27057525647649150063702156134911806760153586215022924037816735921857, 625000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000⟩
def z : QInput := ⟨136721146514200557757045665589, 1000000000000000000000000000000⟩
def a0 : QInput := ⟨1270001, 1280000⟩
def a1 : QInput := ⟨2550001, 2560000⟩
def x0 : QInput := ⟨46255924070167, 100000000000000⟩
def x1 : QInput := ⟨49374602183332977403672480068036683319333956789, 103845937170696552570609926584401920000000000000⟩
def upper : Int := -400000000000000000000000
theorem check_00 : vertexCheck m c p y radius z a0 x0 50 upper=true := by decide
theorem check_01 : vertexCheck m c p y radius z a0 x1 50 upper=true := by decide
theorem check_10 : vertexCheck m c p y radius z a1 x0 50 upper=true := by decide
theorem check_11 : vertexCheck m c p y radius z a1 x1 50 upper=true := by decide
end Spin.Structured.DenseFourierExact.B150
