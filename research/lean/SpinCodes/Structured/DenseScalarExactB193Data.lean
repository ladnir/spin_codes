import SpinCodes.Structured.DenseScalarExactDefs
import SpinCodes.Structured.DenseOccupationFixedVertexDefs

namespace Spin.Structured.DenseScalarExact.B193
open Spin.Structured.DenseOccupationFixed
set_option maxHeartbeats 0
set_option maxRecDepth 100000
def m : QInput := ⟨1, 100⟩
def c : QInput := ⟨58119563056842377037993588567825857659, 170141183460469231731687303715884105728⟩
def p : QInput := ⟨2251513555472729, 2251799813685248⟩
def y : QInput := ⟨2252086108292735, 4503599627370496⟩
def radius : QInput := ⟨88356744200255685402664874480055860043654620815749316566849719008553, 10000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000⟩
def z : QInput := ⟨61797752808988751801433432259, 500000000000000000000000000000⟩
def a0 : QInput := ⟨40950001, 40960000⟩
def a1 : QInput := ⟨81910001, 81920000⟩
def x0 : QInput := ⟨84219921256861908840067417679371902625, 170141183460469231731687303715884105728⟩
def x1 : QInput := ⟨129223289418938324344194200630890939964085310021, 259614842926741381426524816461004800000000000000⟩
def qn : Int := 5070602400912917174410041323815
def qd : Int := 10141204801825835211973625643008
theorem check0 : check qn qd z.num z.den radius.num radius.den 0 = true := by decide
theorem check48 : check qn qd z.num z.den radius.num radius.den 48 = true := by decide
theorem check56 : check qn qd z.num z.den radius.num radius.den 56 = true := by decide
theorem check64 : check qn qd z.num z.den radius.num radius.den 64 = true := by decide
theorem check72 : check qn qd z.num z.den radius.num radius.den 72 = true := by decide
theorem check80 : check qn qd z.num z.den radius.num radius.den 80 = true := by decide
def upper : Int := -400000000000000000000000
theorem check_00 : vertexCheck m c p y radius z a0 x0 50 upper = true := by decide
theorem check_01 : vertexCheck m c p y radius z a0 x1 50 upper = true := by decide
theorem check_10 : vertexCheck m c p y radius z a1 x0 50 upper = true := by decide
theorem check_11 : vertexCheck m c p y radius z a1 x1 50 upper = true := by decide
end Spin.Structured.DenseScalarExact.B193
