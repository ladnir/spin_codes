import SpinCodes.Structured.ConcretePlacementSimplexGood
import SpinCodes.Structured.ConcreteImpulseCoarseProductUniform

/-! A quantitative continuity bound for the actual two-state product integrand. -/
noncomputable section
namespace Spin.Structured.Placement
open Finset ConcreteEncoder FiniteKernel

def timeEmpty (γ t : ℝ) : Matrix (Fin 2) (Fin 2) ℝ := !![1, 0; 0, Real.exp (-γ * t)]

def timeProduct (γ : ℝ) : {a : Nat} → (Fin a → ℝ) → ℝ → Matrix (Fin 2) (Fin 2) ℝ
  | 0, _, last => timeEmpty γ last
  | _ + 1, gaps, last => timeEmpty γ (gaps 0) * impulseMatrix * timeProduct γ (Fin.tail gaps) last

theorem impulseMatrix_substochastic : Substochastic impulseMatrix := by
  constructor
  · intro i j; fin_cases i <;> fin_cases j <;> norm_num [impulseMatrix]
  · intro i; fin_cases i <;> norm_num [impulseMatrix, Fin.sum_univ_two]

theorem timeEmpty_substochastic {γ t : ℝ} (hγ : 0 ≤ γ) (ht : 0 ≤ t) :
    Substochastic (timeEmpty γ t) := by
  have he : Real.exp (-γ * t) ≤ 1 := Real.exp_le_one_iff.mpr (by nlinarith)
  constructor
  · intro i j; fin_cases i <;> fin_cases j <;> simp [timeEmpty, (Real.exp_pos _).le]
  · intro i; fin_cases i <;> simpa [timeEmpty, Fin.sum_univ_two] using he

theorem timeEmpty_rowError {γ s t : ℝ} (hγ : 0 ≤ γ) (hs : 0 ≤ s) (ht : 0 ≤ t) :
    RowError (timeEmpty γ s) (timeEmpty γ t) (γ * |s - t|) := by
  have he := exp_neg_lipschitz (mul_nonneg hγ hs) (mul_nonneg hγ ht)
  have hh : |Real.exp (-γ * s) - Real.exp (-γ * t)| ≤ γ * |s - t| := by
    simpa only [neg_mul, ← mul_sub, abs_mul, abs_of_nonneg hγ] using he
  intro i
  fin_cases i
  · simp [rowAbs, timeEmpty, Matrix.sub_apply, Fin.sum_univ_two, mul_nonneg hγ (abs_nonneg _)]
  · simpa [rowAbs, timeEmpty, Matrix.sub_apply, Fin.sum_univ_two] using hh

theorem timeProduct_substochastic {a : Nat} {γ : ℝ} (hγ : 0 ≤ γ)
    (gaps : Fin a → ℝ) (last : ℝ) (hg : ∀ i, 0 ≤ gaps i) (hl : 0 ≤ last) :
    Substochastic (timeProduct γ gaps last) := by
  induction a with
  | zero => exact timeEmpty_substochastic hγ hl
  | succ a ih =>
    exact ((timeEmpty_substochastic hγ (hg 0)).mul impulseMatrix_substochastic).mul
      (ih (Fin.tail gaps) (fun i => hg i.succ))

theorem timeProduct_rowError {a : Nat} {γ : ℝ} (hγ : 0 ≤ γ)
    (gaps gaps' : Fin a → ℝ) (last last' : ℝ)
    (hg : ∀ i, 0 ≤ gaps i) (hg' : ∀ i, 0 ≤ gaps' i) (hl : 0 ≤ last) (hl' : 0 ≤ last') :
    RowError (timeProduct γ gaps last) (timeProduct γ gaps' last')
      (γ * (|last - last'| + ∑ i, |gaps i - gaps' i|)) := by
  induction a with
  | zero => simpa [timeProduct] using timeEmpty_rowError hγ hl hl'
  | succ a ih =>
    have hfirst := (timeEmpty_rowError hγ (hg 0) (hg' 0)).mul_right impulseMatrix_substochastic
    have htail := ih (Fin.tail gaps) (Fin.tail gaps') (fun i => hg i.succ) (fun i => hg' i.succ)
    have h := hfirst.mul htail
      (timeProduct_substochastic hγ _ last (fun i => hg i.succ) hl)
      ((timeEmpty_substochastic hγ (hg' 0)).mul impulseMatrix_substochastic)
      (by positivity)
    have he : γ * |gaps 0 - gaps' 0| + γ * (|last - last'| + ∑ i, |Fin.tail gaps i - Fin.tail gaps' i|) =
        γ * (|last - last'| + ∑ i, |gaps i - gaps' i|) := by
      rw [Fin.sum_univ_succ]
      change γ * |gaps 0 - gaps' 0| + γ * (|last - last'| + ∑ i : Fin a, |gaps i.succ - gaps' i.succ|) = _
      ring
    rw [he] at h
    exact h

theorem coarseImpulseProduct_timeProduct {a : Nat} (c : ℝ) (gaps : Fin a → Nat) (last : Nat) :
    coarseImpulseProduct c gaps last = timeProduct (c * epochMean) (fun i => (gaps i : ℝ)) last := by
  have he (g : Nat) : coarseEmpty c g = timeEmpty (c * epochMean) g := by
    unfold coarseEmpty timeEmpty
    congr 2
    ring
  induction a with
  | zero => exact he last
  | succ a ih =>
    simp only [coarseImpulseProduct, timeProduct, he, ih]
    rfl

end Spin.Structured.Placement

