import SpinCodes.Structured.DenseOccupationFixedVertexDefs

namespace Spin.Structured.DenseFourierExact.B229
open Spin.Structured.DenseOccupationFixed
set_option maxHeartbeats 0
set_option maxRecDepth 100000
def m : QInput := ⟨-1, 5⟩
def c : QInput := ⟨28224076464601, 62500000000000⟩
def p : QInput := ⟨4490181474054469, 4503599627370496⟩
def y : QInput := ⟨1170917191561415, 2251799813685248⟩
def radius : QInput := ⟨310085764870146626903345769205313098718919294765956109701509405329117, 10000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000⟩
def z : QInput := ⟨5359065748611777250582860579, 40000000000000000000000000000⟩
def a0 : QInput := ⟨2550001, 2560000⟩
def a1 : QInput := ⟨5110001, 5120000⟩
def x0 : QInput := ⟨53744075929833, 100000000000000⟩
def x1 : QInput := ⟨14359431847653, 25000000000000⟩
def upper : Int := -400000000000000000000000
theorem check_00 : vertexCheck m c p y radius z a0 x0 50 upper=true := by decide
theorem check_01 : vertexCheck m c p y radius z a0 x1 50 upper=true := by decide
theorem check_10 : vertexCheck m c p y radius z a1 x0 50 upper=true := by decide
theorem check_11 : vertexCheck m c p y radius z a1 x1 50 upper=true := by decide
end Spin.Structured.DenseFourierExact.B229
