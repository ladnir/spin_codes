import SpinCodes.Structured.DenseOccupationFixedVertexDefs

namespace Spin.Structured.DenseFourierExact.B190
open Spin.Structured.DenseOccupationFixed
set_option maxHeartbeats 0
set_option maxRecDepth 100000
def m : QInput := ⟨1, 100⟩
def c : QInput := ⟨58119563056842377037993588567825857659, 170141183460469231731687303715884105728⟩
def p : QInput := ⟨6343963028079651, 9007199254740992⟩
def y : QInput := ⟨687903258322587, 1125899906842624⟩
def radius : QInput := ⟨3690033842062272904760115527347695672953396537379006059153081413997691, 5000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000⟩
def z : QInput := ⟨37576012902851399154830185989, 250000000000000000000000000000⟩
def a0 : QInput := ⟨10001, 20000⟩
def a1 : QInput := ⟨30001, 40000⟩
def x0 : QInput := ⟨84219921256861908840067417679371902625, 170141183460469231731687303715884105728⟩
def x1 : QInput := ⟨129223289418938324344194200630890939964085310021, 259614842926741381426524816461004800000000000000⟩
def upper : Int := -400000000000000000000000
theorem check_00 : vertexCheck m c p y radius z a0 x0 50 upper=true := by decide
theorem check_01 : vertexCheck m c p y radius z a0 x1 50 upper=true := by decide
theorem check_10 : vertexCheck m c p y radius z a1 x0 50 upper=true := by decide
theorem check_11 : vertexCheck m c p y radius z a1 x1 50 upper=true := by decide
end Spin.Structured.DenseFourierExact.B190
