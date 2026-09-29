import SpinCodes.Structured.ConcreteEpoch

/-! The mean bound of EmptyEpoch applied to the actual emitted output weight. -/
noncomputable section
namespace Spin.Structured.ConcreteEncoder
open ConcreteMaps
open scoped symmDiff

def emptyOutputMean (g : Nat) (q : State) : ℝ :=
  (Spin.piPMF (fun _ : Fin g => transvectionLaw)).expect
    (fun ts => (outputWeight (fun _ => ∅) ts q : ℝ))

theorem emptyOutputMean_zero (q : State) : emptyOutputMean 0 q = 0 := by
  simp only [emptyOutputMean, outputWeight, Nat.cast_zero, Spin.FinPMF.expect_const]

theorem emptyOutputMean_succ (g : Nat) (q : State) :
    emptyOutputMean (g + 1) q = ((Aset q).card : ℝ) +
      transvectionLaw.expect (fun t => emptyOutputMean g (Spin.act t.val.1 t.val.2 q)) := by
  unfold emptyOutputMean
  rw [expect_iid_succ]
  have he (t : Transvection) (ts : Fin g → Transvection) :
      (outputWeight (fun _ => ∅) (Fin.cons t ts) q : ℝ) =
        ((Aset q).card : ℝ) +
          (outputWeight (fun _ => ∅) ts (Spin.act t.val.1 t.val.2 q) : ℝ) := by
    simp only [outputWeight, Fin.cons_zero, Fin.tail_cons, nextState_empty, Nat.cast_add]
    rw [show (∅ : Input) ∆ Aset q = Aset q from bot_symmDiff _]
    rfl
  simp only [he, Spin.FinPMF.expect_add, Spin.FinPMF.expect_const]

theorem emptyOutputMean_eq_sum (g : Nat) (q : State) :
    emptyOutputMean g q = ∑ i ∈ Finset.range g,
      emptyExpect i (fun r => ((Aset r).card : ℝ)) q := by
  induction g generalizing q with
  | zero => simp only [emptyOutputMean_zero, Finset.range_zero, Finset.sum_empty]
  | succ g ih =>
    rw [emptyOutputMean_succ]
    simp_rw [ih]
    rw [Spin.FinPMF.expect_sum]
    simp_rw [← emptyExpect_succ]
    rw [Finset.sum_range_succ']
    rw [emptyExpect_zero]
    ring

/-- The paper's |E S_g - μg| ≤ 2t, for the actual finite encoder with t=128. -/
theorem emptyOutputMean_bound (g : Nat) {q : State} (hq : q ≠ ∅) :
    |emptyOutputMean g q - (128 * (262144 / 524287 : ℝ)) * g| ≤ 256 := by
  let f : State → ℝ := fun r => ((Aset r).card : ℝ)
  have hrec (i : Nat) : emptyExpect (i + 1) f q =
      (emptyExpect i f q + liveMean f) / 2 := by
    rw [emptyExpect_live _ _ hq, emptyExpect_live _ _ hq]
    ring
  have hf0 : 0 ≤ f q := Nat.cast_nonneg _
  have hf1 : f q ≤ 128 := by
    change ((Aset q).card : ℝ) ≤ 128
    have hn : (Aset q).card ≤ 128 := by
      simpa only [Fintype.card_fin] using Finset.card_le_univ (Aset q)
    exact_mod_cast hn
  have hmean : liveMean f = 128 * (262144 / 524287 : ℝ) := emitted_liveMean
  have ht : |emptyExpect 0 f q - liveMean f| ≤ 128 := by
    rw [emptyExpect_zero, hmean, abs_le]
    constructor <;> linarith
  have hb := Spin.mean_sum_bound (m := fun i => emptyExpect i f q)
    (mu := liveMean f) (t := 128) hrec ht g
  rw [hmean] at hb
  simpa only [emptyOutputMean_eq_sum, f, show (2 : ℝ) * 128 = 256 by norm_num] using hb

end Spin.Structured.ConcreteEncoder
