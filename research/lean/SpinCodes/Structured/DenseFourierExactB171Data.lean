import SpinCodes.Structured.DenseOccupationFixedVertexDefs

namespace Spin.Structured.DenseFourierExact.B171
open Spin.Structured.DenseOccupationFixed
set_option maxHeartbeats 0
set_option maxRecDepth 100000
def m : QInput := ⟨7, 100⟩
def c : QInput := ⟨53117391002404332784602253741220096935, 170141183460469231731687303715884105728⟩
def p : QInput := ⟨3829080798801901, 4503599627370496⟩
def y : QInput := ⟨4570397883334507, 9007199254740992⟩
def radius : QInput := ⟨5577789794406192245045080972615375763204150267651171602007760368629933, 10000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000⟩
def z : QInput := ⟨74132833302148096293903502521, 500000000000000000000000000000⟩
def a0 : QInput := ⟨30001, 40000⟩
def a1 : QInput := ⟨70001, 80000⟩
def x0 : QInput := ⟨12032207560909, 25000000000000⟩
def x1 : QInput := ⟨82519260580058937111451454383981630225, 170141183460469231731687303715884105728⟩
def upper : Int := -400000000000000000000000
theorem check_00 : vertexCheck m c p y radius z a0 x0 50 upper=true := by decide
theorem check_01 : vertexCheck m c p y radius z a0 x1 50 upper=true := by decide
theorem check_10 : vertexCheck m c p y radius z a1 x0 50 upper=true := by decide
theorem check_11 : vertexCheck m c p y radius z a1 x1 50 upper=true := by decide
end Spin.Structured.DenseFourierExact.B171
