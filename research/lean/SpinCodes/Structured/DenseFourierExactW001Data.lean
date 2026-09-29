import SpinCodes.Structured.DenseFourierExactColumnDefs

namespace Spin.Structured.DenseFourierExact.W001
set_option maxHeartbeats 0
set_option maxRecDepth 1000000
def qn : Int := 5269567837670373369389975102169
def qd : Int := 10141204801825835211973625643008
def zn : Int := 26856576457192184543084848407
def zd : Int := 200000000000000000000000000000
def rn : Int := 12872242362597962909540481355652147006113511645932736038593337104827
def rd : Int := 400000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
def v : ICoords := ⟨100000000000000000000, 5313190906201121800000, fun i => ([764539050434749790000, 1253412094843392100000, 2035958185374523000000, 3289923300597524800000, 5313190906200753000000]).getD i 0⟩
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
end Spin.Structured.DenseFourierExact.W001
