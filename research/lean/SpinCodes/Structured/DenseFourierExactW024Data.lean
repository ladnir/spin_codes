import SpinCodes.Structured.DenseFourierExactColumnDefs

namespace Spin.Structured.DenseFourierExact.W024
set_option maxHeartbeats 0
set_option maxRecDepth 1000000
def qn : Int := 38077551626671290060295802870337
def qd : Int := 324518553658426726783156020576256
def zn : Int := 727509670879533099545798794729
def zd : Int := 1000000000000000000000000000000
def rn : Int := 2097128403179466924668705376047236427358820910170529033272571862996891784572294634949427689789
def rd : Int := 10000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
def v : ICoords := ⟨100000000000000000000, 407598478642791, fun i => ([407598478642791, 58461359246583, 8143458835921, 948687996265, 89867546310]).getD i 0⟩
theorem check_Z : checkZ qn qd zn zd rn rd v=true := by decide
theorem check_D : checkD qn qd zn zd rn rd v=true := by decide
theorem check_S : ∀ i:Fin 5,checkS qn qd zn zd rn rd v i=true := by
  intro i
  match i with
  | ⟨0, _⟩ => change checkS qn qd zn zd rn rd v (0:Fin 5)=true; decide
  | ⟨1, _⟩ => change checkS qn qd zn zd rn rd v (1:Fin 5)=true; decide
  | ⟨2, _⟩ => change checkS qn qd zn zd rn rd v (2:Fin 5)=true; decide
  | ⟨3, _⟩ => change checkS qn qd zn zd rn rd v (3:Fin 5)=true; decide
  | ⟨4, _⟩ => change checkS qn qd zn zd rn rd v (4:Fin 5)=true; decide
end Spin.Structured.DenseFourierExact.W024
