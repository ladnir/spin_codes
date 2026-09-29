import SpinCodes.Structured.DenseFourierExactColumnDefs

namespace Spin.Structured.DenseFourierExact.W095
set_option maxHeartbeats 0
set_option maxRecDepth 1000000
def qn : Int := 5528575318792230478225397377575
def qd : Int := 40564819207303340847894502572032
def zn : Int := 718345638703320441893352098433
def zd : Int := 1000000000000000000000000000000
def rn : Int := 92962479656809682902312693283777612838112807447516821886959385629846169771573339515039893297
def rd : Int := 2000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
def v : ICoords := ⟨100000000000000000000, 396348547728115, fun i => ([396348547728115, 58280581413135, 8503409072418, 1176938288052, 155354766786]).getD i 0⟩
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
end Spin.Structured.DenseFourierExact.W095
