import SpinCodes.Structured.DenseFourierExactColumnDefs

namespace Spin.Structured.DenseFourierExact.W007
set_option maxHeartbeats 0
set_option maxRecDepth 1000000
def qn : Int := 8981165135689161959402370620679
def qd : Int := 81129638414606681695789005144064
def zn : Int := 740678512027960504523662657083
def zd : Int := 1000000000000000000000000000000
def rn : Int := 503615186152526961321779582643457017781434196959259666720504961470269490137102870567827370477
def rd : Int := 1000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
def v : ICoords := ⟨100000000000000000000, 645845332106979, fun i => ([645845332106979, 99995745269908, 14960999322398, 1823129582317, 170662150800]).getD i 0⟩
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
end Spin.Structured.DenseFourierExact.W007
