import SpinCodes.Numeric.BAEvalDefs
set_option maxRecDepth 4000000
open Spin.Numeric Spin.Numeric.Fix

def show2 (name : String) (x : Fix) : String := s!"{name} {x.lo} {x.hi}"
def showO (name : String) : Option Fix → String
  | none => s!"{name} NONE"
  | some x => show2 name x

def A1 : Fix := ⟨32000000000000000000000000000, 35000000000000000000000000000⟩
def B1 : Fix := ⟨190000000000000000000000000000, 195000000000000000000000000000⟩
def W1 : Fix := ⟨99000000000000000000000000000, 102000000000000000000000000000⟩
def A2 : Fix := ⟨300000000000000000000000000000, 400000000000000000000000000000⟩
def C2 : Fix := ⟨400000000000000000000000000000, 500000000000000000000000000000⟩
def A3 : Fix := ⟨0, 200000000000000000000000000000⟩
def A4 : Fix := ⟨400000000000000000000000000000, 400000000000000000000000000000⟩
def C4 : Fix := ⟨200000000000000000000000000000, 210000000000000000000000000000⟩
def U1 : Int := 330000000000000000000000000000

#eval IO.println (show2 "log2" log2)
#eval IO.println (show2 "flogA_9_10" (flogA 9 10 9))
#eval IO.println (show2 "flogA_5_4" (flogA 5 4 9))
#eval IO.println (show2 "flogA_3_2" (flogA 3 2 9))
#eval IO.println (showO "flogQ_1_3" (flogQ 1 3 10))
#eval IO.println (showO "flogQ_big" (flogQ 4096000000000000000000000000000000 scale 10))
#eval IO.println (show2 "golay1" (fGolayG (ofInt 1)))
#eval IO.println (showO "fxlogx_01_02" (fxlogxL (directLog 10) ⟨100000000000000000000000000000, 200000000000000000000000000000⟩))
#eval IO.println (showO "fxlogx_0_02" (fxlogxL (directLog 10) A3))
#eval IO.println (showO "fpiEval_A2C2" (fpiEvalL (directLog 10) A2 C2))
#eval IO.println (showO "fpiCentered_A2C2" (fpiCenteredG (directLog 10) A2 C2))
#eval IO.println (showO "fpiBest_A2C2" (fpiBestL (directLog 10) A2 C2))
#eval IO.println (showO "fpiBest_A4C4" (fpiBestL (directLog 10) A4 C4))
#eval IO.println (showO "fDpa_A2C2" (fDpaL (directLog 10) A2 C2))
#eval IO.println (showO "fDpc_A2C2" (fDpcL (directLog 10) A2 C2))
#eval IO.println (showO "fgObj_U1A1" (fgObjL (directLog 10) (sc U1) A1))
#eval IO.println (showO "fpiBest_A1B1" (fpiBestL (directLog 10) A1 B1))
#eval IO.println (showO "fpiBest_B1W1" (fpiBestL (directLog 10) B1 W1))
#eval IO.println (showO "fpiClamp_A2C2" (fpiEvalClampL (directLog 10) A2 C2))
#eval IO.println (showO "fpiBestClamp_straddle" (fpiBestClampL (directLog 10) ⟨8000000000000000000000000000, 500000000000000000000000000000⟩ ⟨0, 500000000000000000000000000000⟩))
#eval IO.println s!"infeas_A1B1W1 {infeasible A1 B1 W1}"
#eval IO.println s!"leafOK_A1B1W1 {leafOK (directLog 10) (-76800000000000000000000) U1 A1 B1 W1}"
