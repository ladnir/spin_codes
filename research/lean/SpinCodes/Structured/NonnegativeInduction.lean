import SpinCodes.Structured.Induction

/-! Transfer induction on nonnegative coordinate budgets. A one-step bound
only needs to hold on the budgets actually propagated by a nonnegative
matrix. The older unrestricted induction remains unchanged. -/

namespace Spin.Imt
open Finset
variable {s k : ℕ}

theorem Dominates.D_nonneg {sys : ShellSystem s k} {c : Coords k}
    {μ : Finset (Fin s) → ℝ} (h : Dominates sys c μ) : 0 ≤ c.D := by
  obtain ⟨ν, hν, hsum, _⟩ := h
  exact (Finset.sum_nonneg (fun q _ => hν q)).trans hsum

theorem Coords.eZ_nonneg : (Coords.eZ k).Nonneg := by
  simp [Coords.Nonneg, Coords.eZ]

theorem imt_moment_nonneg {sys : ShellSystem s k} (T : Transfer k) (hT : T.Nonneg)
    (step : (Finset (Fin s) → ℝ) → (Finset (Fin s) → ℝ))
    (hstep : ∀ c μ, c.Nonneg → Dominates sys c μ → Dominates sys (T.apply c) (step μ))
    {c₀ : Coords k} {μ₀ : Finset (Fin s) → ℝ}
    (hc₀ : c₀.Nonneg) (h₀ : Dominates sys c₀ μ₀) (R : ℕ) :
    ∑ q, (step^[R] μ₀) q ≤ ((T.apply)^[R] c₀).total := by
  have key : ∀ n, ((T.apply)^[n] c₀).Nonneg ∧
      Dominates sys ((T.apply)^[n] c₀) (step^[n] μ₀) := by
    intro n
    induction n with
    | zero => exact ⟨hc₀, h₀⟩
    | succ n ih =>
      rw [Function.iterate_succ_apply', Function.iterate_succ_apply']
      exact ⟨Transfer.apply_nonneg hT ih.1, hstep _ _ ih.1 ih.2⟩
  exact (key R).2.total_le

theorem imt_moment_eZ_nonneg {sys : ShellSystem s k} (T : Transfer k) (hT : T.Nonneg)
    (step : (Finset (Fin s) → ℝ) → (Finset (Fin s) → ℝ))
    (hstep : ∀ c μ, c.Nonneg → Dominates sys c μ → Dominates sys (T.apply c) (step μ))
    (R : ℕ) :
    ∑ q, (step^[R] (fun q => if q = ∅ then (1 : ℝ) else 0)) q ≤
      ((T.apply)^[R] (Coords.eZ k)).total :=
  imt_moment_nonneg T hT step hstep Coords.eZ_nonneg (dominates_eZ sys) R

/-- The unrestricted one-step interface cannot hold for a matrix with a
positive zero-to-diffuse entry: the allowed input budget `(-1,0,0)` would
produce a negative diffuse budget. This explains why nonnegativity must be
part of the propagated invariant. -/
theorem unrestricted_step_forces_rowZ_D_nonpos {sys : ShellSystem s k} (T : Transfer k)
    (step : (Finset (Fin s) → ℝ) → (Finset (Fin s) → ℝ))
    (hstep : ∀ c μ, Dominates sys c μ → Dominates sys (T.apply c) (step μ)) :
    T.rowZ.D ≤ 0 := by
  classical
  let c : Coords k := ⟨-1, 0, fun _ => 0⟩
  let μ : Finset (Fin s) → ℝ := fun q => if q = ∅ then -1 else 0
  have hdom : Dominates sys c μ := by
    refine ⟨fun _ => 0, fun _ => le_rfl, ?_, ?_⟩
    · simp [c]
    · intro q
      simp [c, μ]
  have h := (hstep c μ hdom).D_nonneg
  simp only [Transfer.apply, c, neg_one_mul, zero_mul, Finset.sum_const_zero,
    add_zero] at h
  linarith

end Spin.Imt
