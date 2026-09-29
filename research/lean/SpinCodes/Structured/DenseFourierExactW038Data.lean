import SpinCodes.Structured.DenseFourierExactColumnDefs

namespace Spin.Structured.DenseFourierExact.W038
set_option maxHeartbeats 0
set_option maxRecDepth 1000000
def qn : Int := 4630791514299924753805241127603
def qd : Int := 40564819207303340847894502572032
def zn : Int := 734028908017088197304625964779
def zd : Int := 1000000000000000000000000000000
def rn : Int := 797665094007456611009247078822630631653328810946373594041110550717171592422312271415912034773
def rd : Int := 2500000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
def v : ICoords := ⟨100000000000000000000, 508476885282284, fun i => ([508476885282284, 75765973788204, 10937330788318, 1302863221884, 122738707291]).getD i 0⟩
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
end Spin.Structured.DenseFourierExact.W038
