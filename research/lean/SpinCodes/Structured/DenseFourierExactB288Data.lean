import SpinCodes.Structured.DenseOccupationFixedVertexDefs

namespace Spin.Structured.DenseFourierExact.B288
open Spin.Structured.DenseOccupationFixed
set_option maxHeartbeats 0
set_option maxRecDepth 100000
def m : QInput := ⟨-163, 100⟩
def c : QInput := ⟨29249441311893, 20000000000000⟩
def p : QInput := ⟨4268734832555319, 4503599627370496⟩
def y : QInput := ⟨2720582276689713, 4503599627370496⟩
def radius : QInput := ⟨610515968608144564268048857289423108643377917658375332101417262531693, 5000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000⟩
def z : QInput := ⟨17823056619130430678905657409, 125000000000000000000000000000⟩
def a0 : QInput := ⟨70001, 80000⟩
def a1 : QInput := ⟨1, 1⟩
def x0 : QInput := ⟨12530211486671, 15000000000000⟩
def x1 : QInput := ⟨8390967064323, 10000000000000⟩
def upper : Int := -400000000000000000000000
theorem check_00 : vertexCheck m c p y radius z a0 x0 50 upper=true := by decide
theorem check_01 : vertexCheck m c p y radius z a0 x1 50 upper=true := by decide
theorem check_10 : vertexCheck m c p y radius z a1 x0 50 upper=true := by decide
theorem check_11 : vertexCheck m c p y radius z a1 x1 50 upper=true := by decide
end Spin.Structured.DenseFourierExact.B288
