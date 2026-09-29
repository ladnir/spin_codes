import SpinCodes.Structured.DenseFourierExactColumnDefs

namespace Spin.Structured.DenseFourierExact.W036
set_option maxHeartbeats 0
set_option maxRecDepth 1000000
def qn : Int := 35402853539990859172220709008001
def qd : Int := 81129638414606681695789005144064
def zn : Int := 139163099236396518512507940269
def zd : Int := 1000000000000000000000000000000
def rn : Int := 1649667967354308707547242832331479311630012919078782378469466335943607
def rd : Int := 10000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
def v : ICoords := ⟨100000000000000000000, 2126130896444520, fun i => ([2126130896444520, 454145592907195, 95972565110948, 18859542835933, 3419396847962]).getD i 0⟩
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
end Spin.Structured.DenseFourierExact.W036
