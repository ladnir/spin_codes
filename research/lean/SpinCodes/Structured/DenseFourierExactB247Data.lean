import SpinCodes.Structured.DenseOccupationFixedVertexDefs

namespace Spin.Structured.DenseFourierExact.B247
open Spin.Structured.DenseOccupationFixed
set_option maxHeartbeats 0
set_option maxRecDepth 100000
def m : QInput := ⟨-3, 5⟩
def c : QInput := ⟨6909341020092481, 10000000000000000⟩
def p : QInput := ⟨1116640388211527, 1125899906842624⟩
def y : QInput := ⟨2598550854512701, 4503599627370496⟩
def radius : QInput := ⟨1355776172936543696140325061912701584986662326726968558139152348456031, 5000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000⟩
def z : QInput := ⟨75941654400256499514385073097, 500000000000000000000000000000⟩
def a0 : QInput := ⟨310001, 320000⟩
def a1 : QInput := ⟨1, 1⟩
def x0 : QInput := ⟨1244734237944081, 2000000000000000⟩
def x1 : QInput := ⟨1336129721374967, 2000000000000000⟩
def upper : Int := -400000000000000000000000
theorem check_00 : vertexCheck m c p y radius z a0 x0 50 upper=true := by decide
theorem check_01 : vertexCheck m c p y radius z a0 x1 50 upper=true := by decide
theorem check_10 : vertexCheck m c p y radius z a1 x0 50 upper=true := by decide
theorem check_11 : vertexCheck m c p y radius z a1 x1 50 upper=true := by decide
end Spin.Structured.DenseFourierExact.B247
