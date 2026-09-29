import SpinCodes.Structured.DenseOccupationFixedVertexDefs

namespace Spin.Structured.DenseFourierExact.B089
open Spin.Structured.DenseOccupationFixed
set_option maxHeartbeats 0
set_option maxRecDepth 100000
def m : QInput := ⟨4, 5⟩
def c : QInput := ⟨30683842683431, 1250000000000000⟩
def p : QInput := ⟨4946587463656681, 9007199254740992⟩
def y : QInput := ⟨5998491036596071, 9007199254740992⟩
def radius : QInput := ⟨7322597169528805047657608492097560790679607437642299564060995936655516121309, 10000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000⟩
def z : QInput := ⟨363631545815655466409696069, 1600000000000000000000000000⟩
def a0 : QInput := ⟨470081, 1280000⟩
def a1 : QInput := ⟨6001, 16000⟩
def x0 : QInput := ⟨144603079418107, 500000000000000⟩
def x1 : QInput := ⟨663870278625033, 2000000000000000⟩
def upper : Int := -400000000000000000000000
theorem check_00 : vertexCheck m c p y radius z a0 x0 50 upper=true := by decide
theorem check_01 : vertexCheck m c p y radius z a0 x1 50 upper=true := by decide
theorem check_10 : vertexCheck m c p y radius z a1 x0 50 upper=true := by decide
theorem check_11 : vertexCheck m c p y radius z a1 x1 50 upper=true := by decide
end Spin.Structured.DenseFourierExact.B089
