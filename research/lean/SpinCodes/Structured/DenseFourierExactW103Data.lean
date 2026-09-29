import SpinCodes.Structured.DenseFourierExactColumnDefs

namespace Spin.Structured.DenseFourierExact.W103
set_option maxHeartbeats 0
set_option maxRecDepth 1000000
def qn : Int := 22893076960546663211194156617325
def qd : Int := 40564819207303340847894502572032
def zn : Int := 65594999908555291432260744129
def zd : Int := 500000000000000000000000000000
def rn : Int := 6946212725587510547174942130333153149963649028366927936839931643271
def rd : Int := 156250000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
def v : ICoords := ⟨100000000000000000000, 129961796881153390000000000, fun i => ([148007085247510670000000, 952136008089063760000000, 5333919583815572000000000, 26569332445614133000000000, 129961796881121020000000000]).getD i 0⟩
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
end Spin.Structured.DenseFourierExact.W103
