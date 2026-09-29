import SpinCodes.Structured.DenseOccupationFixedVertexDefs

namespace Spin.Structured.DenseFourierExact.B213
open Spin.Structured.DenseOccupationFixed
set_option maxHeartbeats 0
set_option maxRecDepth 100000
def m : QInput := ⟨-1, 10⟩
def c : QInput := ⟨397841147503783, 1000000000000000⟩
def p : QInput := ⟨6343963028079651, 9007199254740992⟩
def y : QInput := ⟨687903258322587, 1125899906842624⟩
def radius : QInput := ⟨3690033842062272904760115527347695672953396537379006059153081413997691, 5000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000⟩
def z : QInput := ⟨37576012902851399154830185989, 250000000000000000000000000000⟩
def a0 : QInput := ⟨10001, 20000⟩
def a1 : QInput := ⟨30001, 40000⟩
def x0 : QInput := ⟨54471334987363575166937446516365236680666043211, 103845937170696552570609926584401920000000000000⟩
def x1 : QInput := ⟨53744075929833, 100000000000000⟩
def upper : Int := -400000000000000000000000
theorem check_00 : vertexCheck m c p y radius z a0 x0 50 upper=true := by decide
theorem check_01 : vertexCheck m c p y radius z a0 x1 50 upper=true := by decide
theorem check_10 : vertexCheck m c p y radius z a1 x0 50 upper=true := by decide
theorem check_11 : vertexCheck m c p y radius z a1 x1 50 upper=true := by decide
end Spin.Structured.DenseFourierExact.B213
