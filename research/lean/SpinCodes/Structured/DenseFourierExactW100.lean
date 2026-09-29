import SpinCodes.Structured.DenseFourierExactW100Data
import SpinCodes.Structured.DenseFourierExactCheck

noncomputable section
namespace Spin.Structured.DenseFourierExact.W100
open Spin.Imt
def wmin : ℝ := 3222932884605
def prefactor : ℝ := 31027640
theorem collatz : ((ConcreteFourier.matrix ((qn:ℝ)/qd) ((zn:ℝ)/zd)).applyCol v.real).le
    (Coords.smul ((rn:ℝ)/rd) v.real) :=
  checked_collatz (by decide) (by decide) (by decide) (by decide) (by decide)
    v (by decide) check_Z check_D check_S
theorem witness_floor : 0 < wmin ∧ wmin ≤ v.real.Z ∧ wmin ≤ v.real.D ∧ ∀ i,wmin ≤ v.real.S i := by
  refine ⟨by norm_num [wmin], by norm_num [wmin,v,ICoords.real], by norm_num [wmin,v,ICoords.real], ?_⟩
  intro i
  fin_cases i <;> norm_num [wmin,v,ICoords.real]
theorem witness_prefactor : v.real.Z/wmin ≤ prefactor := by norm_num [v,ICoords.real,wmin,prefactor]
#print axioms collatz
#print axioms witness_floor
#print axioms witness_prefactor
end Spin.Structured.DenseFourierExact.W100
