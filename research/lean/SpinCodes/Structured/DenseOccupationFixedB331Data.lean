import SpinCodes.Structured.DenseOccupationFixedVertexDefs

namespace Spin.Structured.DenseOccupationFixed.B331
set_option maxHeartbeats 0
set_option maxRecDepth 100000
def m : QInput := ⟨-2, 5⟩
def c : QInput := ⟨14161516955371, 25000000000000⟩
def p : QInput := ⟨3718134143316723, 4611686018427387904⟩
def y : QInput := ⟨7515472203716963, 9007199254740992⟩
def radius : QInput := ⟨459758014223937883, 500000000000000000⟩
def z : QInput := ⟨998645812946312383, 1000000000000000000⟩
def a0 : QInput := ⟨2819, 8192000⟩
def a1 : QInput := ⟨12047, 20480000⟩
def x0 : QInput := ⟨14359431847653, 25000000000000⟩
def x1 : QInput := ⟨1244734237944081, 2000000000000000⟩
def upper : Int := -400000000000000000000000
theorem check_00 : vertexCheck m c p y radius z a0 x0 50 upper = true := by decide
theorem check_01 : vertexCheck m c p y radius z a0 x1 50 upper = true := by decide
theorem check_10 : vertexCheck m c p y radius z a1 x0 50 upper = true := by decide
theorem check_11 : vertexCheck m c p y radius z a1 x1 50 upper = true := by decide
end Spin.Structured.DenseOccupationFixed.B331
