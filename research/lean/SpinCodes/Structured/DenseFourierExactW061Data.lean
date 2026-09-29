import SpinCodes.Structured.DenseFourierExactColumnDefs

namespace Spin.Structured.DenseFourierExact.W061
set_option maxHeartbeats 0
set_option maxRecDepth 1000000
def qn : Int := 9738525602863367032332033138631
def qd : Int := 20282409603651670423947251286016
def zn : Int := 134344346016041137816724063577
def zd : Int := 1000000000000000000000000000000
def rn : Int := 40532412814737869077750591183413660747546602719725177978649935623573
def rd : Int := 1250000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
def v : ICoords := ⟨100000000000000000000, 2376222954891032300, fun i => ([2376222954891017700, 1463047332263222800, 900279914614946980, 551059829847845920, 334149389229523860]).getD i 0⟩
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
end Spin.Structured.DenseFourierExact.W061
