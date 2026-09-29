import SpinCodes.Structured.DenseFourierExactColumnDefs

namespace Spin.Structured.DenseFourierExact.W079
set_option maxHeartbeats 0
set_option maxRecDepth 1000000
def qn : Int := 2901646834970457662170496104427
def qd : Int := 5070602400912917605986812821504
def zn : Int := 75941654400256499514385073097
def zd : Int := 500000000000000000000000000000
def rn : Int := 1355776172936543696140325061912701584986662326726968558139152348456031
def rd : Int := 5000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
def v : ICoords := ⟨100000000000000000000, 228004837970916230000000000, fun i => ([208279381088566970000000, 1273273074517810200000000, 7424198002773716900000000, 41283223703997728000000000, 228004837971066260000000000]).getD i 0⟩
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
end Spin.Structured.DenseFourierExact.W079
