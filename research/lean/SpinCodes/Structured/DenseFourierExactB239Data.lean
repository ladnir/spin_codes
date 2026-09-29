import SpinCodes.Structured.DenseOccupationFixedVertexDefs

namespace Spin.Structured.DenseFourierExact.B239
open Spin.Structured.DenseOccupationFixed
set_option maxHeartbeats 0
set_option maxRecDepth 100000
def m : QInput := ⟨-2, 5⟩
def c : QInput := ⟨14161516955371, 25000000000000⟩
def p : QInput := ⟨2247270629335803, 2251799813685248⟩
def y : QInput := ⟨4887260556201863, 9007199254740992⟩
def radius : QInput := ⟨87904028572394797921961270730054602455215578106741318058487352593811, 1250000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000⟩
def z : QInput := ⟨140735250218194272553824698919, 1000000000000000000000000000000⟩
def a0 : QInput := ⟨1270001, 1280000⟩
def a1 : QInput := ⟨1, 1⟩
def x0 : QInput := ⟨14359431847653, 25000000000000⟩
def x1 : QInput := ⟨1244734237944081, 2000000000000000⟩
def upper : Int := -400000000000000000000000
theorem check_00 : vertexCheck m c p y radius z a0 x0 50 upper=true := by decide
theorem check_01 : vertexCheck m c p y radius z a0 x1 50 upper=true := by decide
theorem check_10 : vertexCheck m c p y radius z a1 x0 50 upper=true := by decide
theorem check_11 : vertexCheck m c p y radius z a1 x1 50 upper=true := by decide
end Spin.Structured.DenseFourierExact.B239
