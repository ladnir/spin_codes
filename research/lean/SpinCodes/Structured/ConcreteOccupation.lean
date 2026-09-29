import SpinCodes.Structured.ConcreteBernoulli
import SpinCodes.Structured.ConcreteTransferAll
import SpinCodes.Structured.SparseMomentBridge

/-! The actual Bernoulli-input occupation step and its sparse finite-round
moment bound. No one-step or low-table hypothesis remains in these statements. -/
noncomputable section
namespace Spin.Structured.ConcreteMaps
open Spin.Imt

theorem bernoulliStep_liveDominates {β z : ℝ}
    (hβ0 : 0 ≤ β) (hβ1 : β ≤ 1) (hz : 0 ≤ z) (hz1 : z ≤ 1)
    (c : Coords 5) (μ : Finset (Fin 19) → ℝ) (hc : c.Nonneg)
    (hμ : LiveDominates actualShellSystem c μ) :
    LiveDominates actualShellSystem ((Occupation.Sparse.numericalMatrix β z).apply c)
      (bernoulliStep β z μ) := by
  rw [bernoulliStep_eq_mixture]
  unfold mixtureStep
  rw [Occupation.Sparse.numericalMatrix, Occupation.matrix_apply]
  apply liveDominates_sum
  intro j
  exact liveDominates_smul (fixedStep_liveDominates hz hz1 j c μ hc hμ)
    (Occupation.probability_nonneg 128 j hβ0 hβ1)

theorem bernoulli_moment_bound {β z : ℝ}
    (hβ0 : 0 ≤ β) (hβ1 : β ≤ 1) (hz : 0 ≤ z) (hz1 : z ≤ 1) (R : Nat) :
    ∑ q, ((bernoulliStep β z)^[R] (fun q => if q = ∅ then (1 : ℝ) else 0)) q ≤
      (((Occupation.Sparse.numericalMatrix β z).apply)^[R] (Coords.eZ 5)).total :=
  live_moment_eZ _ (Occupation.Sparse.numericalMatrix_nonneg hβ0 hβ1 hz)
    (bernoulliStep β z) (bernoulliStep_liveDominates hβ0 hβ1 hz hz1) R

theorem sparse_bernoulli_moment {α : ℝ} (h0 : 0 < α) (h1 : α ≤ 1 / 10000) (R : Nat) :
    ∑ q, ((bernoulliStep ((4 / 5) * α) (1 - (8 / 5) * α))^[R]
      (fun q => if q = ∅ then (1 : ℝ) else 0)) q ≤ 2048 * (1 - 96 * α) ^ R := by
  apply Occupation.Sparse.sparse_live_moment_of_step (sys := actualShellSystem) h0 h1
  exact bernoulliStep_liveDominates (by positivity) (by linarith) (by linarith) (by linarith)

end Spin.Structured.ConcreteMaps
