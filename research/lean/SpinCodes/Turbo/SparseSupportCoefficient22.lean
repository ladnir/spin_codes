import SpinCodes.Turbo.SparseShellCoefficient22
import SpinCodes.Structured.ConcreteMarkedBits

/-! Convert exact row-support choices into chronological shell coefficients.
The old `ConcreteMarkedBits.supportEquiv` is the exact support/Boolean-word
bijection. We reuse it for each of the two separately constrained groups. -/
noncomputable section
namespace Spin.Turbo
open Finset
namespace SparseSupportCoefficient22
open Spin.Structured.ConcreteMarked

theorem support_card (bits : Fin 64 → Bool) :
    (support bits).card = ∑ i, sparseBoolDegree22 (bits i) := by
  simp only [support, Finset.card_filter, sparseBoolDegree22]

def twoSupportEquiv : (Fin 64 → Bool × Bool) ≃
    (Finset (Fin 64) × Finset (Fin 64)) where
  toFun bits := (support (fun i => (bits i).1), support (fun i => (bits i).2))
  invFun AB i := (decide (i ∈ AB.1), decide (i ∈ AB.2))
  left_inv bits := by
    funext i
    apply Prod.ext <;> simp [support]
  right_inv AB := by
    apply Prod.ext <;> ext i <;> simp [support]

def oneSupportCoefficient (W : Bool → Matrix I I ℝ) [Fintype I] [DecidableEq I]
    (v : ℕ) : Matrix I I ℝ :=
  ∑ T : Finset (Fin 64), if T.card = v then
    (List.ofFn (fun i : Fin 64 => W (decide (i ∈ T)))).prod else 0

def twoSupportCoefficient (W : Bool × Bool → Matrix I I ℝ) [Fintype I] [DecidableEq I]
    (u v : ℕ) : Matrix I I ℝ :=
  ∑ AB : Finset (Fin 64) × Finset (Fin 64),
    if AB.1.card = u ∧ AB.2.card = v then
      (List.ofFn (fun i : Fin 64 => W (decide (i ∈ AB.1), decide (i ∈ AB.2)))).prod else 0

theorem oneSupportCoefficient_eq_ordered {I : Type*} [Fintype I] [DecidableEq I]
    (W : Bool → Matrix I I ℝ) (v : ℕ) :
    oneSupportCoefficient W v = orderedMatrixCoefficient sparseBoolDegree22 W 64 v := by
  unfold oneSupportCoefficient orderedMatrixCoefficient
  rw [← Equiv.sum_comp (supportEquiv 64)]
  apply Finset.sum_congr rfl
  intro bits _
  have hinv := (supportEquiv 64).left_inv bits
  change (fun i => decide (i ∈ support bits)) = bits at hinv
  change (if (support bits).card = v then
    (List.ofFn (fun i : Fin 64 => W ((fun j => decide (j ∈ support bits)) i))).prod else 0) = _
  rw [hinv, support_card]

theorem twoSupportCoefficient_eq_shell {I : Type*} [Fintype I] [DecidableEq I]
    (W : Bool × Bool → Matrix I I ℝ) (u v : ℕ) :
    twoSupportCoefficient W u v = sparseTwoShellCoefficient22 W u v := by
  unfold twoSupportCoefficient sparseTwoShellCoefficient22
  rw [← Equiv.sum_comp twoSupportEquiv]
  apply Finset.sum_congr rfl
  intro bits _
  have hfirst := (supportEquiv 64).left_inv (fun i => (bits i).1)
  have hsecond := (supportEquiv 64).left_inv (fun i => (bits i).2)
  change (fun i => decide (i ∈ support (fun j => (bits j).1))) = (fun i => (bits i).1) at hfirst
  change (fun i => decide (i ∈ support (fun j => (bits j).2))) = (fun i => (bits i).2) at hsecond
  change (if (support (fun i => (bits i).1)).card = u ∧
      (support (fun i => (bits i).2)).card = v then
    (List.ofFn (fun i : Fin 64 => W
      ((fun j => decide (j ∈ support (fun k => (bits k).1))) i,
       (fun j => decide (j ∈ support (fun k => (bits k).2))) i))).prod else 0) = _
  rw [hfirst, hsecond, support_card, support_card]

#print axioms oneSupportCoefficient_eq_ordered
#print axioms twoSupportCoefficient_eq_shell
end SparseSupportCoefficient22
end Spin.Turbo

