import SpinCodes.Structured.DenseOccupationFixedVertexDefs

namespace Spin.Structured.DenseFourierExact.B248
open Spin.Structured.DenseOccupationFixed
set_option maxHeartbeats 0
set_option maxRecDepth 100000
def m : QInput := ⟨-4, 5⟩
def c : QInput := ⟨1030683842683431, 1250000000000000⟩
def p : QInput := ⟨689225111364677, 4503599627370496⟩
def y : QInput := ⟨7809055730938139, 9007199254740992⟩
def radius : QInput := ⟨35406972011157146750551596374649886646242758189287212423776493445543125352204710961703941353, 500000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000⟩
def z : QInput := ⟨725907287108328496642662119721, 1000000000000000000000000000000⟩
def a0 : QInput := ⟨2003, 32000⟩
def a1 : QInput := ⟨10007, 80000⟩
def x0 : QInput := ⟨1336129721374967, 2000000000000000⟩
def x1 : QInput := ⟨355396920581893, 500000000000000⟩
def upper : Int := -400000000000000000000000
theorem check_00 : vertexCheck m c p y radius z a0 x0 50 upper=true := by decide
theorem check_01 : vertexCheck m c p y radius z a0 x1 50 upper=true := by decide
theorem check_10 : vertexCheck m c p y radius z a1 x0 50 upper=true := by decide
theorem check_11 : vertexCheck m c p y radius z a1 x1 50 upper=true := by decide
end Spin.Structured.DenseFourierExact.B248
