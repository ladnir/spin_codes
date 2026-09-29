import SpinCodes.Structured.NonnegativeInduction

/-! Domination with diffuse mass restricted to nonzero states. The original
`Dominates` predicate alone allows its diffuse witness to charge zero.
Adding the zero-mass bound is equivalent to choosing a witness that vanishes
at zero, and retains the already proved total-mass inequality. -/

namespace Spin.Imt
open Finset
variable {s k : ℕ}

def LiveDominates (sys : ShellSystem s k) (c : Coords k)
    (μ : Finset (Fin s) → ℝ) : Prop := Dominates sys c μ ∧ μ ∅ ≤ c.Z

theorem LiveDominates.total_le {sys : ShellSystem s k} {c : Coords k}
    {μ : Finset (Fin s) → ℝ} (h : LiveDominates sys c μ) : ∑ q, μ q ≤ c.total := h.1.total_le

theorem liveDominates_iff_witness (sys : ShellSystem s k) (c : Coords k)
    (μ : Finset (Fin s) → ℝ) :
    LiveDominates sys c μ ↔
      ∃ ν : Finset (Fin s) → ℝ, (∀ q, 0 ≤ ν q) ∧ ν ∅ = 0 ∧
        (∑ q, ν q ≤ c.D) ∧ ∀ q, μ q ≤ (if q = ∅ then c.Z else 0) + ν q +
          ∑ i, (if q ∈ sys.shell i then c.S i / ((sys.shell i).card : ℝ) else 0) := by
  classical
  have hempty : ∀ i, (∅ : Finset (Fin s)) ∉ sys.shell i :=
    fun _ hi => sys.ne_empty_of_mem hi rfl
  constructor
  · rintro ⟨⟨ν, hν, hsum, hbound⟩, hzero⟩
    refine ⟨fun q => if q = ∅ then 0 else ν q, ?_, by simp, ?_, ?_⟩
    · intro q
      change 0 ≤ if q = ∅ then 0 else ν q
      split_ifs <;> first | exact le_rfl | exact hν q
    · refine (Finset.sum_le_sum (fun q _ => ?_)).trans hsum
      split_ifs <;> first | exact hν q | exact le_rfl
    · intro q
      by_cases hq : q = ∅
      · subst q
        simpa [hempty] using hzero
      · simpa only [if_neg hq] using hbound q
  · rintro ⟨ν, hν, hzero, hsum, hbound⟩
    refine ⟨⟨ν, hν, hsum, hbound⟩, ?_⟩
    simpa [hempty, hzero] using hbound ∅

theorem liveDominates_eZ (sys : ShellSystem s k) :
    LiveDominates sys (Coords.eZ k) (fun q => if q = ∅ then (1 : ℝ) else 0) := by
  exact ⟨dominates_eZ sys, by simp [Coords.eZ]⟩

theorem LiveDominates.mono {sys : ShellSystem s k} {c c' : Coords k}
    {μ : Finset (Fin s) → ℝ} (h : LiveDominates sys c μ)
    (hZ : c.Z ≤ c'.Z) (hD : c.D ≤ c'.D) (hS : ∀ i, c.S i ≤ c'.S i) :
    LiveDominates sys c' μ := ⟨h.1.mono hZ hD hS, h.2.trans hZ⟩

theorem liveDominates_smul {sys : ShellSystem s k} {c : Coords k}
    {μ : Finset (Fin s) → ℝ} (h : LiveDominates sys c μ) {a : ℝ} (ha : 0 ≤ a) :
    LiveDominates sys (Coords.smul a c) (fun q => a * μ q) :=
  ⟨dominates_smul h.1 ha, mul_le_mul_of_nonneg_left h.2 ha⟩

theorem liveDominates_sum {sys : ShellSystem s k} {ι : Type*} [Fintype ι]
    {cs : ι → Coords k} {μs : ι → Finset (Fin s) → ℝ}
    (h : ∀ i, LiveDominates sys (cs i) (μs i)) :
    LiveDominates sys (Coords.sum cs) (fun q => ∑ i, μs i q) :=
  ⟨dominates_sum (fun i => (h i).1), Finset.sum_le_sum (fun i _ => (h i).2)⟩

theorem live_moment {sys : ShellSystem s k} (T : Transfer k) (hT : T.Nonneg)
    (step : (Finset (Fin s) → ℝ) → (Finset (Fin s) → ℝ))
    (hstep : ∀ (c : Coords k) μ, c.Nonneg → LiveDominates sys c μ →
      LiveDominates sys (T.apply c) (step μ))
    {c₀ : Coords k} {μ₀ : Finset (Fin s) → ℝ}
    (hc₀ : c₀.Nonneg) (h₀ : LiveDominates sys c₀ μ₀) (R : ℕ) :
    ∑ q, (step^[R] μ₀) q ≤ ((T.apply)^[R] c₀).total := by
  have key : ∀ n, ((T.apply)^[n] c₀).Nonneg ∧
      LiveDominates sys ((T.apply)^[n] c₀) (step^[n] μ₀) := by
    intro n
    induction n with
    | zero => exact ⟨hc₀, h₀⟩
    | succ n ih =>
      rw [Function.iterate_succ_apply', Function.iterate_succ_apply']
      exact ⟨Transfer.apply_nonneg hT ih.1, hstep _ _ ih.1 ih.2⟩
  exact (key R).2.total_le

theorem live_moment_eZ {sys : ShellSystem s k} (T : Transfer k) (hT : T.Nonneg)
    (step : (Finset (Fin s) → ℝ) → (Finset (Fin s) → ℝ))
    (hstep : ∀ (c : Coords k) μ, c.Nonneg → LiveDominates sys c μ →
      LiveDominates sys (T.apply c) (step μ)) (R : ℕ) :
    ∑ q, (step^[R] (fun q => if q = ∅ then (1 : ℝ) else 0)) q ≤
      ((T.apply)^[R] (Coords.eZ k)).total :=
  live_moment T hT step hstep Coords.eZ_nonneg (liveDominates_eZ sys) R

end Spin.Imt
