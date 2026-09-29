import SpinCodes.Structured.DenseOccupationFixedVertexDefs

namespace Spin.Structured.DenseFourierExact.B271
open Spin.Structured.DenseOccupationFixed
set_option maxHeartbeats 0
set_option maxRecDepth 100000
def m : QInput := ⟨-7, 5⟩
def c : QInput := ⟨318466707318577, 250000000000000⟩
def p : QInput := ⟨4296353682493619, 4503599627370496⟩
def y : QInput := ⟨5376274426600089, 9007199254740992⟩
def radius : QInput := ⟨392140290789076081042567980375025373936571199034924993285136090480859, 5000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000⟩
def z : QInput := ⟨68814846501743031866678021819, 500000000000000000000000000000⟩
def a0 : QInput := ⟨70001, 80000⟩
def a1 : QInput := ⟨1, 1⟩
def x0 : QInput := ⟨1571386577205347, 2000000000000000⟩
def x1 : QInput := ⟨61043983743719, 75000000000000⟩
def upper : Int := -400000000000000000000000
theorem check_00 : vertexCheck m c p y radius z a0 x0 50 upper=true := by decide
theorem check_01 : vertexCheck m c p y radius z a0 x1 50 upper=true := by decide
theorem check_10 : vertexCheck m c p y radius z a1 x0 50 upper=true := by decide
theorem check_11 : vertexCheck m c p y radius z a1 x1 50 upper=true := by decide
end Spin.Structured.DenseFourierExact.B271
