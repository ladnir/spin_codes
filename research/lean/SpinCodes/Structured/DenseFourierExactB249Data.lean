import SpinCodes.Structured.DenseOccupationFixedVertexDefs

namespace Spin.Structured.DenseFourierExact.B249
open Spin.Structured.DenseOccupationFixed
set_option maxHeartbeats 0
set_option maxRecDepth 100000
def m : QInput := ⟨-4, 5⟩
def c : QInput := ⟨1030683842683431, 1250000000000000⟩
def p : QInput := ⟨541560358341605, 2251799813685248⟩
def y : QInput := ⟨6090569745046833, 9007199254740992⟩
def radius : QInput := ⟨4146446175878725442223657684897664622419585541396181090423580831129402209436819222007413571, 2000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000⟩
def z : QInput := ⟨658593850075895085838459574799, 1000000000000000000000000000000⟩
def a0 : QInput := ⟨10007, 80000⟩
def a1 : QInput := ⟨10003, 40000⟩
def x0 : QInput := ⟨1336129721374967, 2000000000000000⟩
def x1 : QInput := ⟨355396920581893, 500000000000000⟩
def upper : Int := -400000000000000000000000
theorem check_00 : vertexCheck m c p y radius z a0 x0 50 upper=true := by decide
theorem check_01 : vertexCheck m c p y radius z a0 x1 50 upper=true := by decide
theorem check_10 : vertexCheck m c p y radius z a1 x0 50 upper=true := by decide
theorem check_11 : vertexCheck m c p y radius z a1 x1 50 upper=true := by decide
end Spin.Structured.DenseFourierExact.B249
