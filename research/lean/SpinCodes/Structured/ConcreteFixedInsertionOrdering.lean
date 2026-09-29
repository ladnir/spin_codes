import SpinCodes.Structured.ConcreteFixedInsertionCoordinates

/-! Selected-site paths use exactly the increasing subtuple of marked positions. -/
noncomputable section
namespace Spin.Structured.Placement
open Finset ConcreteEncoder
attribute [local instance] Classical.propDecidable

theorem filterMap_choice {α β : Type*} (p : α → Bool) (f : α → β) (xs : List α) :
    xs.filterMap (fun i => if p i then some (f i) else none) = (xs.filter p).map f := by
  induction xs with
  | nil => rfl
  | cons x xs ih => cases h : p x <;> simp [h,ih]

theorem selectedSites_ofFn {Q : Nat} (x : Fin Q → ℝ) (A : Finset (Fin Q)) :
    selectedSites (List.ofFn (fun i => (x i,decide (i ∈ A)))) =
      List.ofFn (fun i : Fin A.card => x (A.orderEmbOfFin rfl i)) := by
  unfold selectedSites
  rw [List.ofFn_eq_map, List.filterMap_map]
  change (List.finRange Q).filterMap (fun i => if decide (i ∈ A) then some (x i) else none) = _
  rw [filterMap_choice]
  have hs : (List.finRange Q).filter (fun i => decide (i ∈ A)) = A.sort := by
    apply ((List.sortedLT_finRange Q).pairwise.filter _).sortedLT.eq_of_mem_iff A.sortedLT_sort
    intro i
    simp
  rw [hs, ← A.listMap_orderEmbOfFin_finRange rfl, List.map_map, ← List.ofFn_eq_map]
  rfl

theorem timeProduct_list_prod {a : Nat} (γ : ℝ) (gaps : Fin a → ℝ) (last : ℝ) :
    timeProduct γ gaps last =
      (List.ofFn (fun i => timeEmpty γ (gaps i) * impulseMatrix)).prod * timeEmpty γ last := by
  induction a with
  | zero => simp [timeProduct]
  | succ a ih => simp only [timeProduct, List.ofFn_succ, List.prod_cons, ih, Fin.tail, Matrix.mul_assoc]

theorem sitePath_ofFn {a : Nat} (γ : ℝ) (x : Fin a → ℝ) :
    sitePath γ impulseMatrix 0 (List.ofFn x) 1 = siteProduct γ x := by
  rw [← factorSitePath_constant, List.map_ofFn, factorSitePath_ofFn]
  symm
  exact timeProduct_list_prod _ _ _

theorem selectedGapProduct_eq_subtuple {Q : Nat} (γ : ℝ) (x : Fin Q → ℝ) (A : Finset (Fin Q)) :
    selectedGapProduct γ impulseMatrix (fun i => positionGaps x i.castSucc) (positionGaps x (Fin.last Q)) A =
      siteProduct γ (fun i : Fin A.card => x (A.orderEmbOfFin rfl i)) := by
  rw [selectedGapProduct_eq_selectedPath, selectedSites_ofFn, sitePath_ofFn]

end Spin.Structured.Placement

