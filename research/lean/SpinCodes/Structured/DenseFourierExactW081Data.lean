import SpinCodes.Structured.DenseFourierExactColumnDefs

namespace Spin.Structured.DenseFourierExact.W081
set_option maxHeartbeats 0
set_option maxRecDepth 1000000
def qn : Int := 21123736265913688168280289983655
def qd : Int := 162259276829213363391578010288128
def zn : Int := 365576028424162091352714126547
def zd : Int := 500000000000000000000000000000
def rn : Int := 29764630106772185940615773792302755677813347449435925180673960280476255944492752246457327561
def rd : Int := 312500000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
def v : ICoords := ⟨100000000000000000000, 548123416028565, fun i => ([548123416028565, 86506693325816, 13503124203679, 1959055708440, 265027298155]).getD i 0⟩
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
end Spin.Structured.DenseFourierExact.W081
