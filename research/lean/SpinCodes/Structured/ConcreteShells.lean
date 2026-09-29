import SpinCodes.Structured.ConcreteWeightSums
import SpinCodes.Structured.Transfer
import SpinCodes.Structured.SparsePolynomialDefs
import SpinCodes.Structured.SparseWitness

/-! The five actual weight shells of the expansion map. Once their checked
cardinalities are supplied, disjointness and the total count prove that they
partition every nonzero state. No separate coverage assumption is needed. -/

noncomputable section
namespace Spin.Structured.ConcreteMaps
open Finset

def shellWeight (i : Fin 5) : ℕ := SparsePolynomial.levels.getD i 0
def shellCount (i : Fin 5) : ℕ := SparsePolynomial.counts.getD i 0
def weightShell (i : Fin 5) : Finset (Finset (Fin 19)) :=
  univ.filter (fun q => (Aset q).card = shellWeight i)

theorem shellWeight_pos (i : Fin 5) : 0 < shellWeight i := by
  fin_cases i <;> decide

theorem shellWeight_injective : Function.Injective shellWeight := by decide

theorem shellCount_pos (i : Fin 5) : 0 < shellCount i := by
  fin_cases i <;> decide

theorem sum_shellCount : ∑ i : Fin 5, shellCount i = 524287 := by decide

theorem shellCount_matrix (i : Fin 5) : (shellCount i : ℝ) = Spin.Imt.Occupation.Sparse.count i := by
  fin_cases i <;> norm_num [shellCount, SparsePolynomial.counts, Spin.Imt.Occupation.Sparse.count]

theorem shellWeight_matrix (i : Fin 5) : shellWeight i = Spin.Imt.Occupation.Sparse.level i := by
  fin_cases i <;> rfl

theorem weightShell_card (i : Fin 5) :
    (weightShell i).card = weightCounts Aset (shellWeight i) := rfl

theorem weightShell_weight {i : Fin 5} {q : Finset (Fin 19)} (hq : q ∈ weightShell i) :
    (Aset q).card = shellWeight i := (Finset.mem_filter.mp hq).2

theorem weightShell_nonzero {i : Fin 5} {q : Finset (Fin 19)} (hq : q ∈ weightShell i) : q ≠ ∅ := by
  intro hzero
  have hw := weightShell_weight hq
  rw [hzero, Aset_empty, Finset.card_empty] at hw
  have hp := shellWeight_pos i
  omega

theorem weightShell_disjoint (i j : Fin 5) (hij : i ≠ j) :
    Disjoint (weightShell i) (weightShell j) := by
  apply Finset.disjoint_left.mpr
  intro q hi hj
  exact hij (shellWeight_injective ((weightShell_weight hi).symm.trans (weightShell_weight hj)))

theorem weightShell_covers
    (hc : ∀ i : Fin 5, weightCounts Aset (shellWeight i) = shellCount i) :
    (univ : Finset (Fin 5)).biUnion weightShell = Spin.nonzeroStates 19 := by
  apply Finset.eq_of_subset_of_card_le
  · intro q hq
    obtain ⟨i, _, hi⟩ := Finset.mem_biUnion.mp hq
    exact Spin.mem_nonzeroStates (weightShell_nonzero hi)
  · rw [Finset.card_biUnion (fun i _ j _ hij => weightShell_disjoint i j hij)]
    simp only [weightShell_card, hc, sum_shellCount]
    rw [Spin.nonzeroStates_eq_erase]
    simp

def shellSystem
    (hc : ∀ i : Fin 5, weightCounts Aset (shellWeight i) = shellCount i) : Spin.Imt.ShellSystem 19 5 where
  shell := weightShell
  nonempty := fun i => Finset.card_pos.mp (by rw [weightShell_card, hc]; exact shellCount_pos i)
  disjoint := weightShell_disjoint
  covers := weightShell_covers hc

theorem shellSystem_card
    (hc : ∀ i : Fin 5, weightCounts Aset (shellWeight i) = shellCount i) (i : Fin 5) :
    (((shellSystem hc).shell i).card : ℝ) = Spin.Imt.Occupation.Sparse.count i := by
  change ((weightShell i).card : ℝ) = _
  rw [weightShell_card, hc, shellCount_matrix]

end Spin.Structured.ConcreteMaps
