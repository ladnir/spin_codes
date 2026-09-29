import SpinCodes.Structured.ConcreteMarkedUniformSubset
open Spin Spin.Structured.ConcreteMarked Finset
example {n Q:ℕ} (marks : Fin Q ↪ Fin n) :
    fairOnMarks (univ.map marks)=(FinPMF.uniform (Finset (Fin Q))).map (fun A => A.map marks) :=
  fairOnMarks_uniform_map marks
#print axioms fairOnMarks_uniform_map
#print axioms fairOnMarks_expect_uniform_sum
