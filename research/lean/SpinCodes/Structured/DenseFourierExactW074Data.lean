import SpinCodes.Structured.DenseFourierExactColumnDefs

namespace Spin.Structured.DenseFourierExact.W074
set_option maxHeartbeats 0
set_option maxRecDepth 1000000
def qn : Int := 5169431599568790023637190774587
def qd : Int := 40564819207303340847894502572032
def zn : Int := 736947745531436273389874598289
def zd : Int := 1000000000000000000000000000000
def rn : Int := 166017594244334855580892783173736280460295335334235628978071646254945181957162698507007956981
def rd : Int := 1250000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
def v : ICoords := ⟨100000000000000000000, 642119574809023, fun i => ([642119574809023, 104706826011434, 16857740657523, 2495805181893, 340092491084]).getD i 0⟩
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
end Spin.Structured.DenseFourierExact.W074
