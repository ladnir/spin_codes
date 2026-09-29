import SpinCodes.Structured.DenseOccupationFixedVertexDefs

namespace Spin.Structured.DenseFourierExact.B175
open Spin.Structured.DenseOccupationFixed
set_option maxHeartbeats 0
set_option maxRecDepth 100000
def m : QInput := ⟨1, 100⟩
def c : QInput := ⟨58119563056842377037993588567825857659, 170141183460469231731687303715884105728⟩
def p : QInput := ⟨5600589304269657, 36028797018963968⟩
def y : QInput := ⟨7209227967457949, 9007199254740992⟩
def radius : QInput := ⟨96507269484952492848914406282455620769791399124285105354020406325637678301970747496003567101, 500000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000⟩
def z : QInput := ⟨148668147224178160710734348573, 200000000000000000000000000000⟩
def a0 : QInput := ⟨2003, 32000⟩
def a1 : QInput := ⟨10007, 80000⟩
def x0 : QInput := ⟨84219921256861908840067417679371902625, 170141183460469231731687303715884105728⟩
def x1 : QInput := ⟨129223289418938324344194200630890939964085310021, 259614842926741381426524816461004800000000000000⟩
def upper : Int := -400000000000000000000000
theorem check_00 : vertexCheck m c p y radius z a0 x0 50 upper=true := by decide
theorem check_01 : vertexCheck m c p y radius z a0 x1 50 upper=true := by decide
theorem check_10 : vertexCheck m c p y radius z a1 x0 50 upper=true := by decide
theorem check_11 : vertexCheck m c p y radius z a1 x1 50 upper=true := by decide
end Spin.Structured.DenseFourierExact.B175
