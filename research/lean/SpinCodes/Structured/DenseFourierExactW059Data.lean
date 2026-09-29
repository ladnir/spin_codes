import SpinCodes.Structured.DenseFourierExactColumnDefs

namespace Spin.Structured.DenseFourierExact.W059
set_option maxHeartbeats 0
set_option maxRecDepth 1000000
def qn : Int := 4782089526781368245411864454747
def qd : Int := 10141204801825835211973625643008
def zn : Int := 136721146514200557757045665589
def zd : Int := 1000000000000000000000000000000
def rn : Int := 27057525647649150063702156134911806760153586215022924037816735921857
def rd : Int := 625000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
def v : ICoords := ⟨100000000000000000000, 492786216649927780, fun i => ([492786216649927950, 246714675605427510, 123387128962449920, 61142651709812656, 29852598913596607]).getD i 0⟩
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
end Spin.Structured.DenseFourierExact.W059
