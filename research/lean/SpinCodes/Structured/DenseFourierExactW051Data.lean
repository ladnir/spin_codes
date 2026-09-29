import SpinCodes.Structured.DenseFourierExactColumnDefs

namespace Spin.Structured.DenseFourierExact.W051
set_option maxHeartbeats 0
set_option maxRecDepth 1000000
def qn : Int := 74815525424603770326674984917205
def qd : Int := 162259276829213363391578010288128
def zn : Int := 139947924148734641672489112037
def zd : Int := 1000000000000000000000000000000
def rn : Int := 161692158644144152003292020206135035172141909141910187917564932116019
def rd : Int := 2500000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
def v : ICoords := ⟨100000000000000000000, 77277123422822189, fun i => ([77277123422822178, 30168043241129435, 11755058104470866, 4510271963426033, 1692748818517828]).getD i 0⟩
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
end Spin.Structured.DenseFourierExact.W051
