import SpinCodes.Structured.DenseOccupationFixedVertexDefs

namespace Spin.Structured.DenseOccupationFixed.B173
set_option maxHeartbeats 0
set_option maxRecDepth 100000
def m : QInput := ⟨9, 100⟩
def c : QInput := ⟨102967997603379096932443767333991244233, 340282366920938463463374607431768211456⟩
def p : QInput := ⟨6332185026843021, 288230376151711744⟩
def y : QInput := ⟨1745505188924709, 2251799813685248⟩
def radius : QInput := ⟨58303800182826281, 500000000000000000⟩
def z : QInput := ⟨239853066408732671, 250000000000000000⟩
def a0 : QInput := ⟨10127, 1280000⟩
def a1 : QInput := ⟨10063, 640000⟩
def x0 : QInput := ⟨49374602183332977403672480068036683319333956789, 103845937170696552570609926584401920000000000000⟩
def x1 : QInput := ⟨81669610035739215919018503711223740925, 170141183460469231731687303715884105728⟩
def upper : Int := -400000000000000000000000
theorem check_00 : vertexCheck m c p y radius z a0 x0 50 upper = true := by decide
theorem check_01 : vertexCheck m c p y radius z a0 x1 50 upper = true := by decide
theorem check_10 : vertexCheck m c p y radius z a1 x0 50 upper = true := by decide
theorem check_11 : vertexCheck m c p y radius z a1 x1 50 upper = true := by decide
end Spin.Structured.DenseOccupationFixed.B173
