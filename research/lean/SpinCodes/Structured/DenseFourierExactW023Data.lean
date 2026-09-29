import SpinCodes.Structured.DenseFourierExactColumnDefs

namespace Spin.Structured.DenseFourierExact.W023
set_option maxHeartbeats 0
set_option maxRecDepth 1000000
def qn : Int := 15369252938610814672879401204045
def qd : Int := 40564819207303340847894502572032
def zn : Int := 26455571036350080901981216427
def zd : Int := 125000000000000000000000000000
def rn : Int := 434766015880792877315826077418935534229729254002174689659500179198945385473
def rd : Int := 10000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
def v : ICoords := ⟨100000000000000000000, 20619921138546, fun i => ([20619921138546, 1621316717289, 127272031846, 9871797281, 759393289]).getD i 0⟩
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
end Spin.Structured.DenseFourierExact.W023
