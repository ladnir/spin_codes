import SpinCodes.Structured.DenseOccupationFixedVertexDefs

namespace Spin.Structured.DenseFourierExact.B064
open Spin.Structured.DenseOccupationFixed
set_option maxHeartbeats 0
set_option maxRecDepth 100000
def m : QInput := ⟨6, 5⟩
def c : QInput := ⟨-832718284462267, 10000000000000000⟩
def p : QInput := ⟨7000766952139315, 9007199254740992⟩
def y : QInput := ⟨1214542177824723, 2251799813685248⟩
def radius : QInput := ⟨63774820549165247542636615921141191414111883496455878980235928923140149, 5000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000⟩
def z : QInput := ⟨21267180149540871149313755529, 125000000000000000000000000000⟩
def a0 : QInput := ⟨50003, 80000⟩
def a1 : QInput := ⟨30001, 40000⟩
def x0 : QInput := ⟨428613422794653, 2000000000000000⟩
def x1 : QInput := ⟨499776708257287, 2000000000000000⟩
def upper : Int := -400000000000000000000000
theorem check_00 : vertexCheck m c p y radius z a0 x0 50 upper=true := by decide
theorem check_01 : vertexCheck m c p y radius z a0 x1 50 upper=true := by decide
theorem check_10 : vertexCheck m c p y radius z a1 x0 50 upper=true := by decide
theorem check_11 : vertexCheck m c p y radius z a1 x1 50 upper=true := by decide
end Spin.Structured.DenseFourierExact.B064
