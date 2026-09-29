import SpinCodes.Structured.ConcreteShellRows
import SpinCodes.Structured.LiveKernel

/-! When every input has zero syndrome, the lazy half preserves the entering
uniform shell instead of consuming diffuse mass. -/

noncomputable section
namespace Spin.Structured.ConcreteMaps
open Finset
open Spin.Imt

def closedShellBudget (i : Fin 5) (m : ℝ) : Coords 5 :=
  ⟨0, 0, fun h => ((actualShellSystem.shell h).card : ℝ) * m / (2 * 524287) +
    if h = i then m / 2 else 0⟩

theorem closedShell_density (i : Fin 5) (m : ℝ) (r : Finset (Fin 19)) :
    (∑ h, if r ∈ actualShellSystem.shell h then
      (closedShellBudget i m).S h / ((actualShellSystem.shell h).card : ℝ) else 0) =
      (if r = ∅ then 0 else m / (2 * 524287)) +
      (if r ∈ weightShell i then m / (2 * shellCount i) else 0) := by
  classical
  have hbase : (∑ h, if r ∈ actualShellSystem.shell h then
      ((((actualShellSystem.shell h).card : ℝ) * m / (2 * 524287)) /
        (actualShellSystem.shell h).card) else 0) =
      if r = ∅ then 0 else m / (2 * 524287) := by
    rw [← actualShellSystem.sum_mem r (m / (2 * 524287))]
    apply Finset.sum_congr rfl
    intro h _
    by_cases hr : r ∈ actualShellSystem.shell h
    · simp only [if_pos hr]
      have hc := (actualShellSystem.card_pos h).ne'
      field_simp
    · simp only [if_neg hr]
  have hlazy : (∑ h, if r ∈ actualShellSystem.shell h then
      (if h = i then m / 2 else 0) / ((actualShellSystem.shell h).card : ℝ) else 0) =
      if r ∈ weightShell i then m / (2 * shellCount i) else 0 := by
    rw [Finset.sum_eq_single i]
    · change (if r ∈ weightShell i then (if i = i then m / 2 else 0) / ((weightShell i).card : ℝ) else 0) = _
      rw [weightShell_card_actual]
      by_cases hr : r ∈ weightShell i
      · simp only [if_pos hr, ite_true]; ring
      · simp only [if_neg hr, ite_true]
    · intro h _ hne
      simp only [if_neg hne, zero_div, ite_self]
    · intro hi
      exact (hi (Finset.mem_univ i)).elim
  have he (h : Fin 5) :
      (if r ∈ actualShellSystem.shell h then (closedShellBudget i m).S h / ((actualShellSystem.shell h).card : ℝ) else 0) =
      (if r ∈ actualShellSystem.shell h then ((((actualShellSystem.shell h).card : ℝ) * m / (2 * 524287)) /
        (actualShellSystem.shell h).card) else 0) +
      (if r ∈ actualShellSystem.shell h then (if h = i then m / 2 else 0) /
        ((actualShellSystem.shell h).card : ℝ) else 0) := by
    dsimp only [closedShellBudget]
    split_ifs <;> ring
  rw [Finset.sum_congr rfl (fun h _ => he h), Finset.sum_add_distrib, hbase, hlazy]

theorem liveDominates_closedShell (i : Fin 5) (m : ℝ) :
    LiveDominates actualShellSystem (closedShellBudget i m)
      (fun r => (if r = ∅ then 0 else m / (2 * 524287)) +
        (if r ∈ weightShell i then m / (2 * shellCount i) else 0)) := by
  classical
  rw [liveDominates_iff_witness]
  refine ⟨fun _ => 0, fun _ => le_rfl, rfl, by simp [closedShellBudget], ?_⟩
  intro r
  change _ ≤ (if r = ∅ then (0 : ℝ) else 0) + 0 + _
  simp only [ite_self, add_zero, zero_add]
  rw [closedShell_density]

theorem shellStepRow_closed (j : Fin 129) (z : ℝ) (i : Fin 5)
    (hl : Occupation.Sparse.live j (SparsePolynomial.weightN j) = 0) (r : Finset (Fin 19)) :
    shellStepRow j z i r =
      (if r = ∅ then 0 else Occupation.hyperMoment (shellWeight i) j z / (2 * 524287)) +
      (if r ∈ weightShell i then Occupation.hyperMoment (shellWeight i) j z / (2 * shellCount i) else 0) := by
  classical
  let m := Occupation.hyperMoment (shellWeight i) j z
  have hc : (shellCount i : ℝ) ≠ 0 := by exact_mod_cast (shellCount_pos i).ne'
  have he (q : Finset (Fin 19)) (hq : q ∈ weightShell i) :
      stepRow j z q r = (if r = q then m / 2 else 0) + (if r = ∅ then 0 else m / (2 * 524287)) := by
    rw [stepRow_all_kernel j z hl, emittedMoment_shell j z hq, Spin.Imt.actLaw,
      Spin.prob_act (weightShell_nonzero hq), nonzeroStates_card_actual]
    change m * ((if r = q then 1 / 2 else 0) + (if r = ∅ then 0 else 1 / (2 * 524287))) = _
    split_ifs <;> ring
  unfold shellStepRow
  rw [Finset.sum_congr rfl he, Finset.sum_add_distrib, Finset.sum_ite_eq,
    Finset.sum_const, nsmul_eq_mul, weightShell_card_actual]
  change ((if r ∈ weightShell i then m / 2 else 0) +
    (shellCount i : ℝ) * (if r = ∅ then 0 else m / (2 * 524287))) / shellCount i =
    (if r = ∅ then 0 else m / (2 * 524287)) + (if r ∈ weightShell i then m / (2 * shellCount i) else 0)
  split_ifs <;> field_simp <;> ring

theorem shellStepRow_dominates_closed {z : ℝ} (hz : 0 ≤ z) (j : Fin 129) (i : Fin 5)
    (hl : Occupation.Sparse.live j (SparsePolynomial.weightN j) = 0)
    (hc : shellCancellationMoment j z i ≤ (Occupation.Sparse.row j (SparsePolynomial.weightN j) z).cancelS i) :
    LiveDominates actualShellSystem
      ((Occupation.fixed Occupation.Sparse.count 524287
        (Occupation.Sparse.row j (SparsePolynomial.weightN j) z)).rowS i) (shellStepRow j z i) := by
  have he := funext (shellStepRow_closed j z i hl)
  rw [he]
  apply (liveDominates_closedShell i (Occupation.hyperMoment (shellWeight i) j z)).mono
  · have hm : 0 ≤ Occupation.hyperMoment (shellWeight i) j z := by
      obtain ⟨q, hq⟩ := actualShellSystem.nonempty i
      rw [← emittedMoment_shell j z hq]
      exact emittedMoment_nonneg hz j q
    have hc0 : 0 ≤ (Occupation.Sparse.row j (SparsePolynomial.weightN j) z).cancelS i :=
      (div_nonneg (Finset.sum_nonneg fun q _ => cancellationMoment_nonneg hz j q)
        (Nat.cast_nonneg (shellCount i))).trans hc
    change 0 ≤ _ / 2 + min (Occupation.hyperMoment (shellWeight i) j z)
      (Occupation.Sparse.live j (SparsePolynomial.weightN j)) / (2 * 524287)
    rw [hl, min_eq_right hm]
    positivity
  · change 0 ≤ if Occupation.Sparse.live j (SparsePolynomial.weightN j) = 0 then 0 else _
    rw [if_pos hl]
  · intro h
    simp only [closedShellBudget, Occupation.fixed, Occupation.Sparse.row, hl, true_and,
      actualShellSystem_card]
    exact le_rfl

theorem shellStepRow_dominates_generic {z : ℝ} (hz : 0 ≤ z) (hz1 : z ≤ 1)
    (j : Fin 129) (hj1 : j.val ≠ 1) (hj2 : j.val ≠ 2) (i : Fin 5) :
    LiveDominates actualShellSystem
      ((Occupation.fixed Occupation.Sparse.count 524287
        (Occupation.Sparse.row j (SparsePolynomial.weightN j) z)).rowS i) (shellStepRow j z i) := by
  have hc := shellCancellationMoment_le_rowS hz hz1 j hj1 hj2 i
  by_cases hl : Occupation.Sparse.live j (SparsePolynomial.weightN j) = 0
  · exact shellStepRow_dominates_closed hz j i hl hc
  · exact shellStepRow_dominates_live hz hz1 j i hl hc

end Spin.Structured.ConcreteMaps
