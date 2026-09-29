import SpinCodes.Structured.DenseFourierExactColumnDefs

namespace Spin.Structured.DenseFourierExact.W078
set_option maxHeartbeats 0
set_option maxRecDepth 1000000
def qn : Int := 5749461390938474061979714525231
def qd : Int := 10141204801825835211973625643008
def zn : Int := 134205669425665989180916433867
def zd : Int := 1000000000000000000000000000000
def rn : Int := 58058119024108778004704772379665260268068448606547372694823392420157
def rd : Int := 1000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
def v : ICoords := ⟨100000000000000000000, 194109238874356840000000000, fun i => ([182642602586478800000000, 1218561008073651500000000, 7158717065505018400000000, 37600525794908562000000000, 194109238874307180000000000]).getD i 0⟩
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
end Spin.Structured.DenseFourierExact.W078
