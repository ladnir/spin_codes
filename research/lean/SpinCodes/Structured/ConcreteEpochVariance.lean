import SpinCodes.Structured.ConcreteEpochMean

/-! A direct second-moment estimate for actual empty epochs. Centering at
the stationary total gives an O(g) bound without assuming covariance laws. -/
noncomputable section
namespace Spin.Structured.ConcreteEncoder
open ConcreteMaps
open scoped symmDiff

def epochMean : ℝ := 128 * (262144 / 524287 : ℝ)
def epochGeom (g : Nat) : ℝ := ∑ i ∈ Finset.range g, (1 / 2 : ℝ) ^ i

theorem epochGeom_nonneg (g : Nat) : 0 ≤ epochGeom g :=
  Finset.sum_nonneg (fun _ _ => by positivity)

theorem epochGeom_le (g : Nat) : epochGeom g ≤ 2 := by
  unfold epochGeom
  rw [geom_sum_eq (by norm_num : (1 / 2 : ℝ) ≠ 1)]
  have hp : 0 ≤ (1 / 2 : ℝ) ^ g := by positivity
  norm_num
  linarith

theorem expect_mul_const {α : Type*} [Fintype α] (P : Spin.FinPMF α)
    (f : α → ℝ) (c : ℝ) : P.expect (fun x => f x * c) = P.expect f * c := by
  simp only [Spin.FinPMF.expect, Finset.sum_mul, mul_assoc]

theorem expect_const_mul {α : Type*} [Fintype α] (P : Spin.FinPMF α)
    (c : ℝ) (f : α → ℝ) : P.expect (fun x => c * f x) = c * P.expect f := by
  simp only [mul_comm c, expect_mul_const]

theorem expect_sub_const {α : Type*} [Fintype α] (P : Spin.FinPMF α)
    (f : α → ℝ) (c : ℝ) : P.expect (fun x => f x - c) = P.expect f - c := by
  simp only [Spin.FinPMF.expect, mul_sub, Finset.sum_sub_distrib,
    ← Finset.sum_mul, P.total, one_mul]

theorem expect_sq_shift {α : Type*} [Fintype α] (P : Spin.FinPMF α)
    (f : α → ℝ) (c : ℝ) :
    P.expect (fun x => (c + f x) ^ 2) =
      c ^ 2 + 2 * c * P.expect f + P.expect (fun x => f x ^ 2) := by
  have he (x) : (c + f x) ^ 2 = c ^ 2 + 2 * c * f x + f x ^ 2 := by ring
  simp only [he, Spin.FinPMF.expect_add, Spin.FinPMF.expect_const, expect_const_mul]

theorem emptyOutputMean_formula (g : Nat) {q : State} (hq : q ≠ ∅) :
    emptyOutputMean g q = epochMean * g + epochGeom g * (((Aset q).card : ℝ) - epochMean) := by
  rw [emptyOutputMean_eq_sum]
  simp only [emptyExpect_live _ _ hq, emitted_liveMean, Finset.sum_add_distrib,
    Finset.sum_const, Finset.card_range, nsmul_eq_mul, ← Finset.sum_mul,
    epochGeom, epochMean]
  ring

theorem emitted_deviation_bound (q : State) : |((Aset q).card : ℝ) - epochMean| ≤ 128 := by
  have h0 : (0 : ℝ) ≤ (Aset q).card := Nat.cast_nonneg _
  have hn : (Aset q).card ≤ 128 := by
    simpa only [Fintype.card_fin] using Finset.card_le_univ (Aset q)
  have h1 : ((Aset q).card : ℝ) ≤ 128 := by exact_mod_cast hn
  rw [abs_le]
  unfold epochMean
  constructor <;> linarith

theorem tailMean_deviation (g : Nat) {q : State} (hq : q ≠ ∅) :
    transvectionLaw.expect (fun t => emptyOutputMean g (Spin.act t.val.1 t.val.2 q) -
      epochMean * g) = epochGeom g * (((Aset q).card : ℝ) - epochMean) / 2 := by
  have he (t : Transvection) :
      emptyOutputMean g (Spin.act t.val.1 t.val.2 q) - epochMean * g =
        epochGeom g * (((Aset (Spin.act t.val.1 t.val.2 q)).card : ℝ) - epochMean) := by
    rw [emptyOutputMean_formula g (transvection_nonzero t hq)]
    ring
  simp only [he, expect_const_mul, expect_sub_const]
  rw [transvection_expect_live hq (fun r => ((Aset r).card : ℝ)), emitted_liveMean]
  unfold epochMean
  ring

def emptyCenteredSecond (g : Nat) (q : State) : ℝ :=
  (Spin.piPMF (fun _ : Fin g => transvectionLaw)).expect
    (fun ts => ((outputWeight (fun _ => ∅) ts q : ℝ) - epochMean * g) ^ 2)

theorem emptyCenteredSecond_zero (q : State) : emptyCenteredSecond 0 q = 0 := by
  simp only [emptyCenteredSecond, outputWeight, Nat.cast_zero, mul_zero, sub_self,
    zero_pow (by decide : 2 ≠ 0), Spin.FinPMF.expect_const]

theorem emptyCenteredSecond_succ (g : Nat) {q : State} (hq : q ≠ ∅) :
    emptyCenteredSecond (g + 1) q =
      transvectionLaw.expect (fun t => emptyCenteredSecond g (Spin.act t.val.1 t.val.2 q)) +
        (1 + epochGeom g) * (((Aset q).card : ℝ) - epochMean) ^ 2 := by
  unfold emptyCenteredSecond
  rw [expect_iid_succ]
  have he (t : Transvection) (ts : Fin g → Transvection) :
      (outputWeight (fun _ => ∅) (Fin.cons t ts) q : ℝ) - epochMean * (g + 1 : Nat) =
        (((Aset q).card : ℝ) - epochMean) +
          ((outputWeight (fun _ => ∅) ts (Spin.act t.val.1 t.val.2 q) : ℝ) - epochMean * g) := by
    simp only [outputWeight, Fin.cons_zero, Fin.tail_cons, nextState_empty, Nat.cast_add,
      Nat.cast_one]
    rw [show (∅ : Input) ∆ Aset q = Aset q from bot_symmDiff _]
    change ((Aset q).card : ℝ) +
      (outputWeight (fun _ => ∅) ts (Spin.act t.val.1 t.val.2 q) : ℝ) -
      epochMean * ((g : ℝ) + 1) = _
    ring
  simp only [he, expect_sq_shift, expect_sub_const, Spin.FinPMF.expect_add,
    Spin.FinPMF.expect_const, expect_const_mul]
  change (((Aset q).card : ℝ) - epochMean) ^ 2 +
    2 * (((Aset q).card : ℝ) - epochMean) *
      (transvectionLaw.expect (fun t => emptyOutputMean g (Spin.act t.val.1 t.val.2 q)) - epochMean * g) +
    transvectionLaw.expect (fun t => emptyCenteredSecond g (Spin.act t.val.1 t.val.2 q)) =
      transvectionLaw.expect (fun t => emptyCenteredSecond g (Spin.act t.val.1 t.val.2 q)) +
        (1 + epochGeom g) * (((Aset q).card : ℝ) - epochMean) ^ 2
  rw [← expect_sub_const, tailMean_deviation g hq]
  ring

/-- A centered second moment for the actual emitted sum; no covariance
hypothesis and no stationary initial-state assumption. -/
theorem emptyCenteredSecond_bound (g : Nat) {q : State} (hq : q ≠ ∅) :
    emptyCenteredSecond g q ≤ 3 * 128 ^ 2 * g := by
  induction g generalizing q with
  | zero => simp only [emptyCenteredSecond_zero, Nat.cast_zero, mul_zero, le_refl]
  | succ g ih =>
    rw [emptyCenteredSecond_succ g hq]
    have htail : transvectionLaw.expect
        (fun t => emptyCenteredSecond g (Spin.act t.val.1 t.val.2 q)) ≤ 3 * 128 ^ 2 * g := by
      calc _ ≤ transvectionLaw.expect (fun _ => 3 * 128 ^ 2 * (g : ℝ)) :=
        Spin.FinPMF.expect_mono _ (fun t => ih (transvection_nonzero t hq))
        _ = _ := Spin.FinPMF.expect_const _ _
    have hd : (((Aset q).card : ℝ) - epochMean) ^ 2 ≤ 128 ^ 2 := by
      have hh := emitted_deviation_bound q
      nlinarith [sq_abs (((Aset q).card : ℝ) - epochMean),
        abs_nonneg (((Aset q).card : ℝ) - epochMean)]
    have hg := epochGeom_le g
    have hterm := mul_le_mul_of_nonneg_right (show 1 + epochGeom g ≤ 3 by linarith)
      (sq_nonneg (((Aset q).card : ℝ) - epochMean))
    push_cast
    nlinarith

end Spin.Structured.ConcreteEncoder
