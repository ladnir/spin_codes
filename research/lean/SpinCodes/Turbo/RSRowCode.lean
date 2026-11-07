import SpinCodes.Turbo.RSRow
import Mathlib.RingTheory.Polynomial.Basic
import Mathlib.Data.Finset.Card

/-!
# Four-symbol messages as degree-below-four RS polynomials

The coefficient equivalence makes the row message space exactly `F^4`.
The information-set lemma from `RSRow` then applies to the actual row
encoding, with no random-code assumption.
-/

noncomputable section

namespace Spin.Turbo

variable {F : Type*} [Field F]

abbrev RSRow (F : Type*) [Field F] := Polynomial.degreeLT F 4

noncomputable instance [Fintype F] : Fintype (RSRow F) :=
  Fintype.ofEquiv (Fin 4 → F)
    ((Polynomial.degreeLTEquiv F 4).symm.toEquiv)

/-- A row has four field symbols of message dimension. -/
theorem rsRow_card [Fintype F] :
    Fintype.card (RSRow F) = (Fintype.card F) ^ 4 := by
  rw [Fintype.card_congr (Polynomial.degreeLTEquiv F 4).toEquiv]
  simp

theorem rsRow_natDegree_lt (p : RSRow F) :
    (p.1 : Polynomial F).natDegree < 4 := by
  by_cases hp : (p.1 : Polynomial F) = 0
  · simp [hp]
  · exact (Polynomial.natDegree_lt_iff_degree_lt hp).2
      (Polynomial.mem_degreeLT.mp p.2)

/-- Evaluation at any four distinct positions is injective on row messages. -/
theorem rsRow_eq_of_four_evals_code
    (a : Fin 4 → F) (ha : Function.Injective a)
    (p q : RSRow F)
    (heval : ∀ i : Fin 4, p.1.eval (a i) = q.1.eval (a i)) : p = q := by
  apply Subtype.ext
  exact rsRow_eq_of_four_evals a ha p.1 q.1
    (rsRow_natDegree_lt p) (rsRow_natDegree_lt q) heval

/-- Evaluation at any four distinct positions is also onto the four-symbol
space. This identifies the systematic interpolating encoder with our
polynomial message model. -/
noncomputable def rsRow_four_eval_equiv [Fintype F]
    (a : Fin 4 → F) (ha : Function.Injective a) :
    RSRow F ≃ (Fin 4 → F) := by
  let f : RSRow F → Fin 4 → F := fun p i => p.1.eval (a i)
  have hinj : Function.Injective f := by
    intro p q hpq
    exact rsRow_eq_of_four_evals_code a ha p q (fun i => congrFun hpq i)
  exact Equiv.ofBijective f ⟨hinj,
    hinj.surjective_of_finite (Polynomial.degreeLTEquiv F 4).toEquiv⟩

/-- A nonzero RS row has at most three zeros among eight distinct points. -/
theorem rsRow_zero_positions_le_three
    [DecidableEq F]
    (points : Fin 8 → F) (hpoints : Function.Injective points)
    (p : RSRow F) (hp : p ≠ 0) :
    (Finset.univ.filter (fun i : Fin 8 => p.1.eval (points i) = 0)).card ≤ 3 := by
  classical
  let zeros : Finset (Fin 8) :=
    Finset.univ.filter (fun i => p.1.eval (points i) = 0)
  change zeros.card ≤ 3
  by_contra h
  have hfour : 4 ≤ zeros.card := by omega
  obtain ⟨s, hs, hsc⟩ := Finset.exists_subset_card_eq hfour
  have hinj : Function.Injective (fun i : s => points i.1) := by
    intro a b hab
    apply Subtype.ext
    exact hpoints hab
  have heval : ∀ i : s, p.1.eval (points i.1) = 0 := by
    intro i
    have hi : i.1 ∈ zeros := hs i.2
    simpa only [zeros, Finset.mem_filter, Finset.mem_univ, true_and] using hi
  have hdeg : p.1.natDegree < Fintype.card s := by
    simpa [Fintype.card_coe, hsc] using (rsRow_natDegree_lt p)
  have hpzero : p.1 = 0 :=
    Polynomial.eq_zero_of_natDegree_lt_card_of_eval_eq_zero p.1 hinj heval hdeg
  exact hp (Subtype.ext hpzero)

/-- Consequently a nonzero RS[8,4] row has at least five active columns. -/
theorem rsRow_support_ge_five
    [DecidableEq F]
    (points : Fin 8 → F) (hpoints : Function.Injective points)
    (p : RSRow F) (hp : p ≠ 0) :
    5 ≤ (Finset.univ.filter (fun i : Fin 8 => p.1.eval (points i) ≠ 0)).card := by
  classical
  have hpart := Finset.card_filter_add_card_filter_not
    (s := (Finset.univ : Finset (Fin 8)))
    (p := fun i : Fin 8 => p.1.eval (points i) = 0)
  have hsum :
      (Finset.univ.filter (fun i : Fin 8 => p.1.eval (points i) = 0)).card +
      (Finset.univ.filter (fun i : Fin 8 => p.1.eval (points i) ≠ 0)).card = 8 := by
    simpa using hpart
  have hzero := rsRow_zero_positions_le_three points hpoints p hp
  omega

#print axioms rsRow_card
#print axioms rsRow_natDegree_lt
#print axioms rsRow_eq_of_four_evals_code
#print axioms rsRow_four_eval_equiv
#print axioms rsRow_zero_positions_le_three
#print axioms rsRow_support_ge_five

end Spin.Turbo
