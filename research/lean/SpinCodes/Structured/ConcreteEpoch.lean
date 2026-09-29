import SpinCodes.Structured.ConcreteEncoderMoment
import SpinCodes.Structured.EmptyEpoch

/-! Exact finite empty-epoch laws for the concrete encoder experiment. -/
noncomputable section
namespace Spin.Structured.ConcreteEncoder
open ConcreteMaps
open scoped symmDiff

def terminalState : {R : Nat} → (Fin R → Input) → (Fin R → Transvection) → State → State
  | 0, _, _, q => q
  | _ + 1, xs, ts, q =>
      terminalState (Fin.tail xs) (Fin.tail ts) (nextState q (xs 0) (ts 0))

def emptyState {g : Nat} (ts : Fin g → Transvection) (q : State) : State :=
  terminalState (fun _ => ∅) ts q

def liveMean (f : State → ℝ) : ℝ :=
  (∑ q ∈ Spin.nonzeroStates 19, f q) / 524287

def emptyExpect (g : Nat) (f : State → ℝ) (q : State) : ℝ :=
  (Spin.piPMF (fun _ : Fin g => transvectionLaw)).expect (fun ts => f (emptyState ts q))

theorem nextState_empty (q : State) (t : Transvection) :
    nextState q ∅ t = Spin.act t.val.1 t.val.2 q := by
  rw [nextState, Cset_empty, Spin.symmDiff_empty_right]

theorem act_nonzero_of_mem {s : Nat} (p : Finset (Fin s) × Finset (Fin s))
    (hp : p ∈ Spin.transPairs s) {q : Finset (Fin s)} (hq : q ≠ ∅) :
    Spin.act p.1 p.2 q ≠ ∅ := by
  have hn := Finset.filter_eq_empty_iff.mp
    (Finset.card_eq_zero.mp (Spin.count_target_zero hq))
  exact hn hp

theorem transvection_nonzero_generic {s : Nat}
    (t : {p : Finset (Fin s) × Finset (Fin s) // p ∈ Spin.transPairs s})
    {q : Finset (Fin s)} (hq : q ≠ ∅) : Spin.act t.val.1 t.val.2 q ≠ ∅ :=
  act_nonzero_of_mem t.val t.property hq

theorem transvection_nonzero (t : Transvection) {q : State} (hq : q ≠ ∅) :
    Spin.act t.val.1 t.val.2 q ≠ ∅ :=
  @transvection_nonzero_generic 19 t q hq

theorem emptyState_succ {g : Nat} (t : Transvection) (ts : Fin g → Transvection)
    (q : State) : emptyState (Fin.cons t ts) q =
      emptyState ts (Spin.act t.val.1 t.val.2 q) := by
  simp only [emptyState, terminalState, Fin.cons_zero, Fin.tail_cons, nextState_empty]
  rfl

theorem emptyState_zero {g : Nat} (ts : Fin g → Transvection) : emptyState ts ∅ = ∅ := by
  induction g with
  | zero => rfl
  | succ g ih =>
    rw [← Fin.cons_self_tail ts, emptyState_succ, Spin.act_empty]
    exact ih _

theorem emptyState_nonzero {g : Nat} (ts : Fin g → Transvection) {q : State}
    (hq : q ≠ ∅) : emptyState ts q ≠ ∅ := by
  induction g generalizing q with
  | zero => exact hq
  | succ g ih =>
    rw [← Fin.cons_self_tail ts, emptyState_succ]
    exact ih _ (transvection_nonzero _ hq)

theorem transvection_expect_live {q : State} (hq : q ≠ ∅) (f : State → ℝ) :
    transvectionLaw.expect (fun t => f (Spin.act t.val.1 t.val.2 q)) =
      (f q + liveMean f) / 2 := by
  have he : transvectionLaw.expect (fun t => f (Spin.act t.val.1 t.val.2 q)) =
      (∑ p ∈ Spin.transPairs 19, f (Spin.act p.1 p.2 q)) /
        ((Spin.transPairs 19).card : ℝ) := by
    unfold transvectionLaw Spin.FinPMF.expect Spin.FinPMF.uniform
    simp only [Fintype.card_coe, div_mul_eq_mul_div, one_mul, ← Finset.sum_div]
    apply congrArg (fun a : ℝ => a / ((Spin.transPairs 19).card : ℝ))
    exact Finset.sum_coe_sort (Spin.transPairs 19)
      (fun p : State × State => f (Spin.act p.1 p.2 q))
  rw [he, Spin.expect_act hq, nonzeroStates_card_actual, liveMean]
  push_cast
  ring

theorem expect_affine {α : Type*} [Fintype α] (P : Spin.FinPMF α)
    (a b : ℝ) (f : α → ℝ) :
    P.expect (fun x => a + b * (f x - a)) = a + b * (P.expect f - a) := by
  unfold Spin.FinPMF.expect
  simp only [mul_add, mul_sub, Finset.sum_add_distrib, Finset.sum_sub_distrib]
  simp_rw [← mul_assoc, mul_comm (P.p _) b, mul_assoc]
  rw [← Finset.sum_mul, P.total, one_mul, ← Finset.mul_sum,
    ← Finset.mul_sum, ← Finset.sum_mul, P.total, one_mul]

theorem emptyExpect_zero (f : State → ℝ) (q : State) : emptyExpect 0 f q = f q := by
  exact Spin.FinPMF.expect_const (Spin.piPMF (fun _ : Fin 0 => transvectionLaw)) (f q)

theorem emptyExpect_succ (g : Nat) (f : State → ℝ) (q : State) :
    emptyExpect (g + 1) f q =
      transvectionLaw.expect (fun t => emptyExpect g f (Spin.act t.val.1 t.val.2 q)) := by
  unfold emptyExpect
  rw [expect_iid_succ]
  simp only [emptyState_succ]

/-- Actual finite empty epochs have the exact lazy/uniform mixing law. -/
theorem emptyExpect_live (g : Nat) (f : State → ℝ) {q : State} (hq : q ≠ ∅) :
    emptyExpect g f q = liveMean f + (1 / 2 : ℝ) ^ g * (f q - liveMean f) := by
  induction g generalizing q with
  | zero => rw [emptyExpect_zero]; ring
  | succ g ih =>
    rw [emptyExpect_succ]
    have he : (fun t : Transvection => emptyExpect g f (Spin.act t.val.1 t.val.2 q)) =
        fun t => liveMean f + (1 / 2 : ℝ) ^ g *
          (f (Spin.act t.val.1 t.val.2 q) - liveMean f) := by
      funext t
      exact ih (transvection_nonzero t hq)
    rw [he, expect_affine, transvection_expect_live hq]
    ring

theorem emptyExpect_dead (g : Nat) (f : State → ℝ) : emptyExpect g f ∅ = f ∅ := by
  unfold emptyExpect
  simp only [emptyState_zero]
  exact Spin.FinPMF.expect_const _ _

theorem emitted_liveMean : liveMean (fun q => ((Aset q).card : ℝ)) =
    128 * (262144 / 524287 : ℝ) := by
  unfold liveMean
  rw [← weightShell_covers actual_shell_counts,
    Finset.sum_biUnion (fun i _ j _ hij => weightShell_disjoint i j hij)]
  have he (i : Fin 5) : (∑ q ∈ weightShell i, ((Aset q).card : ℝ)) =
      (shellCount i : ℝ) * shellWeight i := by
    rw [Finset.sum_congr rfl (fun q hq => by rw [weightShell_weight hq])]
    rw [Finset.sum_const, nsmul_eq_mul, weightShell_card_actual]
  simp only [he]
  norm_num [shellCount, shellWeight, SparsePolynomial.counts, SparsePolynomial.levels,
    Fin.sum_univ_succ]

end Spin.Structured.ConcreteEncoder
