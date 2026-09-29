import SpinCodes.Structured.MapSpectrumBridge
import SpinCodes.Structured.FiberNumericsAll
import SpinCodes.Structured.ConcreteShells

/-! Identification of the frozen numerical data with the actual table maps. -/
noncomputable section
namespace Spin.Structured.ConcreteMaps

theorem transpose_spectrum (w : Fin 129) :
    weightCounts CtransposeSet w = FiberNumerics.Data.spectrum.getD w 0 := by
  rw [CtransposeSet_weightCounts w.isLt]
  rfl

theorem actual_shell_counts (i : Fin 5) :
    weightCounts Aset (shellWeight i) = shellCount i := by
  fin_cases i
  · change weightCounts Aset 48 = 5166
    exact Aset_weightCounts (i := 48) (by decide)
  · change weightCounts Aset 56 = 110288
    exact Aset_weightCounts (i := 56) (by decide)
  · change weightCounts Aset 64 = 293455
    exact Aset_weightCounts (i := 64) (by decide)
  · change weightCounts Aset 72 = 110128
    exact Aset_weightCounts (i := 72) (by decide)
  · change weightCounts Aset 80 = 5250
    exact Aset_weightCounts (i := 80) (by decide)

def actualShellSystem : Spin.Imt.ShellSystem 19 5 := shellSystem actual_shell_counts

theorem kernel_card_frozen (j : Fin 129) :
    (kernelLayer j).card = (SparsePolynomial.weightN j).kernel :=
  FiberNumerics.kernel_counts_frozen transpose_spectrum j

theorem syndromeFiber_cap_frozen (j : Fin 129) {q : Finset (Fin 19)} (hq : q ≠ ∅) :
    (syndromeFiber j q).card ≤ (SparsePolynomial.weightN j).cap :=
  FiberNumerics.fiber_caps_frozen transpose_spectrum j hq

theorem pair_card_frozen (j : Fin 129) :
    (equalSyndromePairs j).card = FiberNumerics.pairTotals.getD j 0 :=
  FiberNumerics.pair_counts transpose_spectrum j

theorem actualShellSystem_card (i : Fin 5) :
    ((actualShellSystem.shell i).card : ℝ) = Spin.Imt.Occupation.Sparse.count i :=
  shellSystem_card actual_shell_counts i

end Spin.Structured.ConcreteMaps
