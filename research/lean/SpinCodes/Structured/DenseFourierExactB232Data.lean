import SpinCodes.Structured.DenseOccupationFixedVertexDefs

namespace Spin.Structured.DenseFourierExact.B232
open Spin.Structured.DenseOccupationFixed
set_option maxHeartbeats 0
set_option maxRecDepth 100000
def m : QInput := ⟨-2, 5⟩
def c : QInput := ⟨14161516955371, 25000000000000⟩
def p : QInput := ⟨2785485255891209, 18014398509481984⟩
def y : QInput := ⟨1855845974641443, 2251799813685248⟩
def radius : QInput := ⟨166017594244334855580892783173736280460295335334235628978071646254945181957162698507007956981, 1250000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000⟩
def z : QInput := ⟨736947745531436273389874598289, 1000000000000000000000000000000⟩
def a0 : QInput := ⟨2003, 32000⟩
def a1 : QInput := ⟨10007, 80000⟩
def x0 : QInput := ⟨14359431847653, 25000000000000⟩
def x1 : QInput := ⟨1244734237944081, 2000000000000000⟩
def upper : Int := -400000000000000000000000
theorem check_00 : vertexCheck m c p y radius z a0 x0 50 upper=true := by decide
theorem check_01 : vertexCheck m c p y radius z a0 x1 50 upper=true := by decide
theorem check_10 : vertexCheck m c p y radius z a1 x0 50 upper=true := by decide
theorem check_11 : vertexCheck m c p y radius z a1 x1 50 upper=true := by decide
end Spin.Structured.DenseFourierExact.B232
