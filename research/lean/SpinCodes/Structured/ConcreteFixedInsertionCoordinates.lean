import SpinCodes.Structured.ConcreteFixedInsertionDeletion

/-! Site coordinates and gap products agree for arbitrary inserted matrices. -/
noncomputable section
namespace Spin.Structured.Placement
open Finset ConcreteEncoder

def relativeGaps {a : Nat} (start last : ℝ) (x : Fin a → ℝ) (i : Fin (a+1)) : ℝ :=
  (Fin.snoc x last : Fin (a+1) → ℝ) i - (Fin.cons start x : Fin (a+1) → ℝ) i

theorem relativeGaps_zero {a : Nat} (start last : ℝ) (x : Fin (a+1) → ℝ) :
    relativeGaps start last x 0 = x 0-start := by
  simp [relativeGaps, Fin.snoc_apply_zero]

theorem relativeGaps_succ {a : Nat} (start last : ℝ) (x : Fin (a+1) → ℝ) (i : Fin (a+1)) :
    relativeGaps start last x i.succ = relativeGaps (x 0) last (Fin.tail x) i := by
  unfold relativeGaps
  rw [Fin.cons_succ, Fin.cons_self_tail]
  congr 1
  refine Fin.lastCases ?_ (fun j => ?_) i
  · simp only [Fin.succ_last, Fin.snoc_last]
  · simp only [← Fin.castSucc_succ, Fin.snoc_castSucc, Fin.tail]

def factorSitePath (γ start : ℝ) : List (ℝ × Matrix (Fin 2) (Fin 2) ℝ) → ℝ → Matrix (Fin 2) (Fin 2) ℝ
  | [], last => timeEmpty γ (last-start)
  | (x,P) :: xs, last => timeEmpty γ (x-start) * P * factorSitePath γ x xs last

theorem factorSitePath_ofFn {a : Nat} (γ start last : ℝ) (x : Fin a → ℝ)
    (P : Fin a → Matrix (Fin 2) (Fin 2) ℝ) :
    factorSitePath γ start (List.ofFn (fun i => (x i,P i))) last =
      (List.ofFn (fun i => timeEmpty γ (relativeGaps start last x i.castSucc) * P i)).prod *
        timeEmpty γ (relativeGaps start last x (Fin.last a)) := by
  induction a generalizing start with
  | zero => simp [factorSitePath, relativeGaps, Fin.snoc, Fin.cons]
  | succ a ih =>
    simp only [List.ofFn_succ, factorSitePath, List.prod_cons, Fin.castSucc_zero, relativeGaps_zero]
    rw [ih]
    simp only [Fin.castSucc_succ, relativeGaps_succ, ← Fin.succ_last, Fin.tail]
    simp only [Matrix.mul_assoc]
    rfl

theorem factorSitePath_constant (γ start last : ℝ) (P : Matrix (Fin 2) (Fin 2) ℝ) (xs : List ℝ) :
    factorSitePath γ start (xs.map (fun x => (x,P))) last = sitePath γ P start xs last := by
  induction xs generalizing start with
  | nil => rfl
  | cons x xs ih => simp only [List.map_cons, factorSitePath, sitePath, ih]

theorem factorSitePath_marked (γ start last : ℝ) (P : Matrix (Fin 2) (Fin 2) ℝ) (xs : List (ℝ × Bool)) :
    factorSitePath γ start (xs.map (fun p => (p.1,if p.2 then P else 1))) last = markedSitePath γ P start xs last := by
  induction xs generalizing start with
  | nil => rfl
  | cons p xs ih =>
    rcases p with ⟨x,b⟩
    simp only [List.map_cons, factorSitePath, markedSitePath, ih]

theorem selectedGapProduct_eq_selectedPath {Q : Nat} (γ : ℝ) (P : Matrix (Fin 2) (Fin 2) ℝ)
    (x : Fin Q → ℝ) (A : Finset (Fin Q)) :
    selectedGapProduct γ P (fun i => positionGaps x i.castSucc) (positionGaps x (Fin.last Q)) A =
      sitePath γ P 0 (selectedSites (List.ofFn (fun i => (x i, decide (i ∈ A))))) 1 := by
  classical
  rw [← markedSitePath_eq_selected]
  rw [← factorSitePath_marked, List.map_ofFn]
  simp only [Function.comp_def, decide_eq_true_eq]
  rw [factorSitePath_ofFn]
  rfl

end Spin.Structured.Placement


