import SpinCodes.Structured.DenseOccupationFixedVertexDefs

namespace Spin.Structured.DenseOccupationFixed.B163
set_option maxHeartbeats 0
set_option maxRecDepth 100000
def m : QInput := ⟨7, 100⟩
def c : QInput := ⟨53117391002404332784602253741220096935, 170141183460469231731687303715884105728⟩
def p : QInput := ⟨5386065771040615, 2305843009213693952⟩
def y : QInput := ⟨7081728982588543, 9007199254740992⟩
def radius : QInput := ⟨396648920424326979, 500000000000000000⟩
def z : QInput := ⟨996021574517471523, 1000000000000000000⟩
def a0 : QInput := ⟨11023, 10240000⟩
def a1 : QInput := ⟨10511, 5120000⟩
def x0 : QInput := ⟨81669610035739215919018503711223740925, 170141183460469231731687303715884105728⟩
def x1 : QInput := ⟨12032207560909, 25000000000000⟩
def upper : Int := -400000000000000000000000
theorem check_00 : vertexCheck m c p y radius z a0 x0 50 upper = true := by decide
theorem check_01 : vertexCheck m c p y radius z a0 x1 50 upper = true := by decide
theorem check_10 : vertexCheck m c p y radius z a1 x0 50 upper = true := by decide
theorem check_11 : vertexCheck m c p y radius z a1 x1 50 upper = true := by decide
end Spin.Structured.DenseOccupationFixed.B163
