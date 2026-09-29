import SpinCodes.Structured.ConcreteRouteTranspose
import SpinCodes.FiniteLaw

namespace Spin.Structured.ConcreteRoute

open Finset

/-- Transposing the independent row law gives independent region laws. -/
theorem transpose_bernoulli {L b : ℕ} (p : Fin L → ℝ)
    (hp0 : ∀ i, 0 ≤ p i) (hp1 : ∀ i, p i ≤ 1) :
    (piPMF (fun i => poissonBinom (fun _ : Fin b => p i)
      (fun _ => hp0 i) (fun _ => hp1 i))).map transpose =
    piPMF (fun _ : Fin b => poissonBinom p hp0 hp1) := by
  apply FinPMF.ext
  funext regions
  rw [FinPMF.map_p]
  have he : ∀ rows : Fin L → Finset (Fin b),
      transpose rows = regions ↔ rows = transpose regions :=
    fun rows => (transposeEquiv L b).apply_eq_iff_eq_symm_apply
  simp only [he, Finset.sum_ite_eq', Finset.mem_univ, if_true]
  exact transpose_density p hp0 hp1 regions

end Spin.Structured.ConcreteRoute

