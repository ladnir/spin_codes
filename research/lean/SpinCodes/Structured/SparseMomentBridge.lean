import SpinCodes.Structured.LiveInduction
import SpinCodes.Structured.SparseIteration

/-! The numerical sparse iterate bound applied to a state-measure evolution.
Only one-step domination on nonnegative budgets with live diffuse mass is required. The actual
encoder's one-step law is a separate obligation. -/

noncomputable section
namespace Spin.Imt.Occupation.Sparse

theorem sparse_live_moment_of_step {sys : ShellSystem 19 5} {α : ℝ}
    (h0 : 0 < α) (h1 : α ≤ 1 / 10000)
    (step : (Finset (Fin 19) → ℝ) → (Finset (Fin 19) → ℝ))
    (hstep : ∀ (c : Coords 5) μ, c.Nonneg → LiveDominates sys c μ →
      LiveDominates sys ((numericalMatrix ((4 / 5) * α) (1 - (8 / 5) * α)).apply c) (step μ))
    (R : ℕ) :
    ∑ q, (step^[R] (fun q => if q = ∅ then (1 : ℝ) else 0)) q ≤
      2048 * (1 - 96 * α) ^ R := by
  have hT := numericalMatrix_nonneg (β := (4 / 5) * α) (z := 1 - (8 / 5) * α)
    (by positivity) (by linarith) (by linarith)
  exact (live_moment_eZ _ hT step hstep R).trans (sparse_iterate_bound h0 h1 R)

end Spin.Imt.Occupation.Sparse
