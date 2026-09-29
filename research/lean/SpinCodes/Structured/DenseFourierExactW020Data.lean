import SpinCodes.Structured.DenseFourierExactColumnDefs

namespace Spin.Structured.DenseFourierExact.W020
set_option maxHeartbeats 0
set_option maxRecDepth 1000000
def qn : Int := 36312508218925475067862830783361
def qd : Int := 324518553658426726783156020576256
def zn : Int := 369217956000168725636606306593
def zd : Int := 500000000000000000000000000000
def rn : Int := 2150623201060033281000664943609887320214972102241377892890527886725676502451142237479532460507
def rd : Int := 5000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
def v : ICoords := ⟨100000000000000000000, 594747734673326, fun i => ([594747734673326, 90909061092585, 13439489954553, 1625199266354, 152457777195]).getD i 0⟩
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
end Spin.Structured.DenseFourierExact.W020
