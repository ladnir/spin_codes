import SpinCodes.Structured.DenseFourierExactColumnDefs

namespace Spin.Structured.DenseFourierExact.W080
set_option maxHeartbeats 0
set_option maxRecDepth 1000000
def qn : Int := 10982997105863807271536081200989
def qd : Int := 20282409603651670423947251286016
def zn : Int := 140735250218194272553824698919
def zd : Int := 1000000000000000000000000000000
def rn : Int := 87904028572394797921961270730054602455215578106741318058487352593811
def rd : Int := 1250000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
def v : ICoords := ⟨100000000000000000000, 445111114441560680000000, fun i => ([7584787851536803800000, 21572739092921540000000, 59898466447026294000000, 163458640772159310000000, 445111114441497650000000]).getD i 0⟩
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
end Spin.Structured.DenseFourierExact.W080
