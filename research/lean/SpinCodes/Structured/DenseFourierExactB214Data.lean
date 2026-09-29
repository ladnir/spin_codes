import SpinCodes.Structured.DenseOccupationFixedVertexDefs

namespace Spin.Structured.DenseFourierExact.B214
open Spin.Structured.DenseOccupationFixed
set_option maxHeartbeats 0
set_option maxRecDepth 100000
def m : QInput := ⟨-9, 100⟩
def c : QInput := ⟨3339835265656588966103687050071259581601, 8507059173023461586584365185794205286400⟩
def p : QInput := ⟨6343963028079651, 9007199254740992⟩
def y : QInput := ⟨687903258322587, 1125899906842624⟩
def radius : QInput := ⟨3690033842062272904760115527347695672953396537379006059153081413997691, 5000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000⟩
def z : QInput := ⟨37576012902851399154830185989, 250000000000000000000000000000⟩
def a0 : QInput := ⟨10001, 20000⟩
def a1 : QInput := ⟨30001, 40000⟩
def x0 : QInput := ⟨88471573424730015812668800004660364803, 170141183460469231731687303715884105728⟩
def x1 : QInput := ⟨54471334987363575166937446516365236680666043211, 103845937170696552570609926584401920000000000000⟩
def upper : Int := -400000000000000000000000
theorem check_00 : vertexCheck m c p y radius z a0 x0 50 upper=true := by decide
theorem check_01 : vertexCheck m c p y radius z a0 x1 50 upper=true := by decide
theorem check_10 : vertexCheck m c p y radius z a1 x0 50 upper=true := by decide
theorem check_11 : vertexCheck m c p y radius z a1 x1 50 upper=true := by decide
end Spin.Structured.DenseFourierExact.B214
