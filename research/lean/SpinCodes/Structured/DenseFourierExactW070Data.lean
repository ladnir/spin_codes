import SpinCodes.Structured.DenseFourierExactColumnDefs

namespace Spin.Structured.DenseFourierExact.W070
set_option maxHeartbeats 0
set_option maxRecDepth 1000000
def qn : Int := 5257630681200953454732868713635
def qd : Int := 10141204801825835211973625643008
def zn : Int := 5359065748611777250582860579
def zd : Int := 40000000000000000000000000000
def rn : Int := 310085764870146626903345769205313098718919294765956109701509405329117
def rd : Int := 10000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
def v : ICoords := ⟨100000000000000000000, 4191668454753048700000, fun i => ([677153536676358940000, 1077680580579792700000, 1700481495494691900000, 2670493456803231800000, 4191668454752682000000]).getD i 0⟩
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
end Spin.Structured.DenseFourierExact.W070
