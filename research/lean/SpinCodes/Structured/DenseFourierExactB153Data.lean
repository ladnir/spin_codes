import SpinCodes.Structured.DenseOccupationFixedVertexDefs

namespace Spin.Structured.DenseFourierExact.B153
open Spin.Structured.DenseOccupationFixed
set_option maxHeartbeats 0
set_option maxRecDepth 100000
def m : QInput := ⟨1, 10⟩
def c : QInput := ⟨297841147503783, 1000000000000000⟩
def p : QInput := ⟨8994219463258877, 9007199254740992⟩
def y : QInput := ⟨4327626887123875, 9007199254740992⟩
def radius : QInput := ⟨328175647358845962607755343381275496537438477183591505373918065608707, 10000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000⟩
def z : QInput := ⟨134443316966951039460176022307, 1000000000000000000000000000000⟩
def a0 : QInput := ⟨5110001, 5120000⟩
def a1 : QInput := ⟨10230001, 10240000⟩
def x0 : QInput := ⟨94384754313803, 200000000000000⟩
def x1 : QInput := ⟨49374602183332977403672480068036683319333956789, 103845937170696552570609926584401920000000000000⟩
def upper : Int := -400000000000000000000000
theorem check_00 : vertexCheck m c p y radius z a0 x0 50 upper=true := by decide
theorem check_01 : vertexCheck m c p y radius z a0 x1 50 upper=true := by decide
theorem check_10 : vertexCheck m c p y radius z a1 x0 50 upper=true := by decide
theorem check_11 : vertexCheck m c p y radius z a1 x1 50 upper=true := by decide
end Spin.Structured.DenseFourierExact.B153
