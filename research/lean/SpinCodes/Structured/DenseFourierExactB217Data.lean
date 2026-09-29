import SpinCodes.Structured.DenseOccupationFixedVertexDefs

namespace Spin.Structured.DenseFourierExact.B217
open Spin.Structured.DenseOccupationFixed
set_option maxHeartbeats 0
set_option maxRecDepth 100000
def m : QInput := ⟨-1, 10⟩
def c : QInput := ⟨397841147503783, 1000000000000000⟩
def p : QInput := ⟨8993777074240001, 9007199254740992⟩
def y : QInput := ⟨1170619456446307, 2251799813685248⟩
def radius : QInput := ⟨158203392309202796562687919547201970848109694137228991731959619586933, 5000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000⟩
def z : QInput := ⟨67071607883642631428690298457, 500000000000000000000000000000⟩
def a0 : QInput := ⟨5110001, 5120000⟩
def a1 : QInput := ⟨10230001, 10240000⟩
def x0 : QInput := ⟨54471334987363575166937446516365236680666043211, 103845937170696552570609926584401920000000000000⟩
def x1 : QInput := ⟨105615245686197, 200000000000000⟩
def upper : Int := -400000000000000000000000
theorem check_00 : vertexCheck m c p y radius z a0 x0 50 upper=true := by decide
theorem check_01 : vertexCheck m c p y radius z a0 x1 50 upper=true := by decide
theorem check_10 : vertexCheck m c p y radius z a1 x0 50 upper=true := by decide
theorem check_11 : vertexCheck m c p y radius z a1 x1 50 upper=true := by decide
end Spin.Structured.DenseFourierExact.B217
