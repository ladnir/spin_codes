import SpinCodes.Structured.DenseFourierExactColumnDefs

namespace Spin.Structured.DenseFourierExact.W071
set_option maxHeartbeats 0
set_option maxRecDepth 1000000
def qn : Int := 10528290430046087248732288126307
def qd : Int := 20282409603651670423947251286016
def zn : Int := 67071607883642631428690298457
def zd : Int := 500000000000000000000000000000
def rn : Int := 158203392309202796562687919547201970848109694137228991731959619586933
def rd : Int := 5000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
def v : ICoords := ⟨100000000000000000000, 4770595265557473400000, fun i => ([723505204625434220000, 1170253191557864600000, 1876009384236520100000, 2992416977346347900000, 4770595265557388100000]).getD i 0⟩
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
end Spin.Structured.DenseFourierExact.W071
