import SpinCodes.Structured.ConcreteFixedInsertionExpansion

/-! Unused potential impulse sites can be deleted exactly, merging their empty gaps. -/
noncomputable section
namespace Spin.Structured.Placement
open Finset ConcreteEncoder

def sitePath (γ : ℝ) (P : Matrix (Fin 2) (Fin 2) ℝ) (start : ℝ) : List ℝ → ℝ → Matrix (Fin 2) (Fin 2) ℝ
  | [], last => timeEmpty γ (last-start)
  | x :: xs, last => timeEmpty γ (x-start) * P * sitePath γ P x xs last

def markedSitePath (γ : ℝ) (P : Matrix (Fin 2) (Fin 2) ℝ) (start : ℝ) :
    List (ℝ × Bool) → ℝ → Matrix (Fin 2) (Fin 2) ℝ
  | [], last => timeEmpty γ (last-start)
  | (x,bit) :: xs, last => timeEmpty γ (x-start) * (if bit then P else 1) * markedSitePath γ P x xs last

def selectedSites (xs : List (ℝ × Bool)) : List ℝ :=
  xs.filterMap (fun p => if p.2 then some p.1 else none)

theorem sitePath_prepend_empty (γ : ℝ) (P : Matrix (Fin 2) (Fin 2) ℝ)
    (start mid last : ℝ) (xs : List ℝ) :
    timeEmpty γ (mid-start) * sitePath γ P mid xs last = sitePath γ P start xs last := by
  cases xs with
  | nil =>
    simp only [sitePath, ← timeEmpty_add]
    congr 1
    ring
  | cons x xs =>
    simp only [sitePath]
    rw [← Matrix.mul_assoc, ← Matrix.mul_assoc, ← timeEmpty_add]
    rw [show mid-start+(x-mid)=x-start by ring]

theorem markedSitePath_eq_selected (γ : ℝ) (P : Matrix (Fin 2) (Fin 2) ℝ)
    (start last : ℝ) (xs : List (ℝ × Bool)) :
    markedSitePath γ P start xs last = sitePath γ P start (selectedSites xs) last := by
  induction xs generalizing start with
  | nil => rfl
  | cons p xs ih =>
    rcases p with ⟨x,bit⟩
    cases bit
    · simp only [markedSitePath, Bool.false_eq_true, ↓reduceIte, Matrix.mul_one, selectedSites,
        List.filterMap_cons, List.filterMap_nil]
      rw [ih]
      exact sitePath_prepend_empty γ P start x last (selectedSites xs)
    · simp only [markedSitePath, ↓reduceIte, selectedSites, List.filterMap_cons, sitePath]
      rw [ih]
      rfl

end Spin.Structured.Placement

