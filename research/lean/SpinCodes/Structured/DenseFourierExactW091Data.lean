import SpinCodes.Structured.DenseFourierExactColumnDefs

namespace Spin.Structured.DenseFourierExact.W091
set_option maxHeartbeats 0
set_option maxRecDepth 1000000
def qn : Int := 11624412566171833419637402007811
def qd : Int := 20282409603651670423947251286016
def zn : Int := 28695330300993484403841014761
def zd : Int := 200000000000000000000000000000
def rn : Int := 132427562082508544580943484877309667582931547850404698899097425139761
def rd : Int := 1000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
def v : ICoords := ⟨100000000000000000000, 447601513220611030000000000, fun i => ([289451146865344210000000, 2033258763045898000000000, 13067530440316271000000000, 76986010505214904000000000, 447601513219769760000000000]).getD i 0⟩
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
end Spin.Structured.DenseFourierExact.W091
