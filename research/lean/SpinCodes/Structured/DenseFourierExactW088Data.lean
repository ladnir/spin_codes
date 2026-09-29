import SpinCodes.Structured.DenseFourierExactColumnDefs

namespace Spin.Structured.DenseFourierExact.W088
set_option maxHeartbeats 0
set_option maxRecDepth 1000000
def qn : Int := 5382197305808808003199756716103
def qd : Int := 40564819207303340847894502572032
def zn : Int := 725907287108328496642662119721
def zd : Int := 1000000000000000000000000000000
def rn : Int := 35406972011157146750551596374649886646242758189287212423776493445543125352204710961703941353
def rd : Int := 500000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
def v : ICoords := ⟨100000000000000000000, 478186333328335, fun i => ([478186333328335, 73296477947083, 11127619359405, 1584267004409, 212482481196]).getD i 0⟩
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
end Spin.Structured.DenseFourierExact.W088
