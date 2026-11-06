import SpinCodes.UniformMixer.RegionalContinuum

/-!
# Stability of the continuum regional product

The two-state product is substochastic. Moving its gap lengths changes
each row by at most `γ` times the `ℓ₁` change in the gaps. This is the
analytic estimate needed to compare rounded bit sites with the exact
finite encoder gaps.
-/

noncomputable section

namespace Spin.UniformMixer

open Spin.Structured.FiniteKernel (RowError rowAbs)
open Spin.Structured.ConcreteEncoder (exp_neg_lipschitz)

theorem coarseP_substochastic :
    Spin.Structured.FiniteKernel.Substochastic coarseP := by
  have hD : (0 : ℝ) < liveCount State22 := by
    exact_mod_cast (liveCount_pos (S := State22))
  have hD1 : (1 : ℝ) ≤ liveCount State22 := by
    exact_mod_cast (liveCount_pos (S := State22))
  have hfrac0 : (0 : ℝ) ≤ 1 / liveCount State22 := by positivity
  have hfrac1 : (1 : ℝ) / liveCount State22 ≤ 1 := by
    apply (div_le_iff₀ hD).mpr
    nlinarith
  constructor
  · intro b d
    cases b <;> cases d
    · change (0 : ℝ) ≤ 0; norm_num
    · change (0 : ℝ) ≤ 1; norm_num
    · change (0 : ℝ) ≤ 1 / liveCount State22; exact hfrac0
    · change (0 : ℝ) ≤ 1 - 1 / liveCount State22
      linarith
  · intro b
    cases b
    · rw [Fintype.sum_bool]
      change (1 : ℝ) + 0 ≤ 1
      norm_num
    · rw [Fintype.sum_bool]
      change (1 : ℝ) - 1 / liveCount State22 +
        1 / liveCount State22 ≤ 1
      linarith

theorem timeRegionGap_substochastic {γ u : ℝ}
    (hγ : 0 ≤ γ) (hu : 0 ≤ u) :
    Spin.Structured.FiniteKernel.Substochastic (timeRegionGap γ u) := by
  have he1 : Real.exp (-γ * u) ≤ 1 :=
    Real.exp_le_one_iff.mpr (by nlinarith)
  constructor
  · intro b d
    cases b <;> cases d <;>
      simp [timeRegionGap, coarseD, (Real.exp_pos _).le]
  · intro b
    cases b <;>
      simpa [timeRegionGap, coarseD, Fintype.sum_bool] using he1

theorem timeRegionGap_rowError {γ s t : ℝ}
    (hγ : 0 ≤ γ) (hs : 0 ≤ s) (ht : 0 ≤ t) :
    RowError (timeRegionGap γ s) (timeRegionGap γ t)
      (γ * |s - t|) := by
  have he := exp_neg_lipschitz
    (mul_nonneg hγ hs) (mul_nonneg hγ ht)
  have hh : |Real.exp (-γ * s) - Real.exp (-γ * t)| ≤
      γ * |s - t| := by
    simpa only [neg_mul, ← mul_sub, abs_mul, abs_of_nonneg hγ]
      using he
  intro b
  cases b
  · simp [rowAbs, timeRegionGap, coarseD, Matrix.sub_apply,
      Fintype.sum_bool, mul_nonneg hγ (abs_nonneg _)]
  · simpa [rowAbs, timeRegionGap, coarseD, Matrix.sub_apply,
      Fintype.sum_bool] using hh

theorem timeRegionProduct_substochastic {a : ℕ} {γ : ℝ}
    (hγ : 0 ≤ γ) (gaps : Fin a → ℝ) (last : ℝ)
    (hg : ∀ i, 0 ≤ gaps i) (hl : 0 ≤ last) :
    Spin.Structured.FiniteKernel.Substochastic
      (timeRegionProduct γ gaps last) := by
  induction a with
  | zero => exact timeRegionGap_substochastic hγ hl
  | succ a ih =>
      exact ((timeRegionGap_substochastic hγ (hg 0)).mul
        coarseP_substochastic).mul
        (ih (Fin.tail gaps) (fun i => hg i.succ))

theorem timeRegionProduct_rowError {a : ℕ} {γ : ℝ}
    (hγ : 0 ≤ γ) (gaps gaps' : Fin a → ℝ)
    (last last' : ℝ)
    (hg : ∀ i, 0 ≤ gaps i) (hg' : ∀ i, 0 ≤ gaps' i)
    (hl : 0 ≤ last) (hl' : 0 ≤ last') :
    RowError (timeRegionProduct γ gaps last)
      (timeRegionProduct γ gaps' last')
      (γ * (|last - last'| + ∑ i, |gaps i - gaps' i|)) := by
  induction a with
  | zero => simpa [timeRegionProduct] using
      timeRegionGap_rowError hγ hl hl'
  | succ a ih =>
      have hfirst := (timeRegionGap_rowError hγ (hg 0) (hg' 0)).mul_right
        coarseP_substochastic
      have htail := ih (Fin.tail gaps) (Fin.tail gaps')
        (fun i => hg i.succ) (fun i => hg' i.succ)
      have h := hfirst.mul htail
        (timeRegionProduct_substochastic hγ _ last
          (fun i => hg i.succ) hl)
        ((timeRegionGap_substochastic hγ (hg' 0)).mul
          coarseP_substochastic)
        (by positivity)
      have he : γ * |gaps 0 - gaps' 0| +
          γ * (|last - last'| +
            ∑ i, |Fin.tail gaps i - Fin.tail gaps' i|) =
          γ * (|last - last'| +
            ∑ i, |gaps i - gaps' i|) := by
        rw [Fin.sum_univ_succ]
        change γ * |gaps 0 - gaps' 0| +
          γ * (|last - last'| +
            ∑ i : Fin a, |gaps i.succ - gaps' i.succ|) = _
        ring
      rw [he] at h
      exact h

end Spin.UniformMixer
