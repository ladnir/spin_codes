import SpinCodes.Structured.DenseFourierExactColumnDefs

namespace Spin.Structured.DenseFourierExact.W090
set_option maxHeartbeats 0
set_option maxRecDepth 1000000
def qn : Int := 23550546501560155828281522213403
def qd : Int := 40564819207303340847894502572032
def zn : Int := 156167681658821452432281700591
def zd : Int := 1000000000000000000000000000000
def rn : Int := 886106221632040163476806057044589617512891840507445608334194761656801
def rd : Int := 2000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
def v : ICoords := ⟨100000000000000000000, 1165721696129988700000000000, fun i => ([503881421611554380000000, 3736860893702488800000000, 26358918550174584000000000, 176002759510870390000000000, 1165721696128656900000000000]).getD i 0⟩
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
end Spin.Structured.DenseFourierExact.W090
