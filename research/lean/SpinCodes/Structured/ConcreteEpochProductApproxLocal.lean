import SpinCodes.Structured.ConcreteEpochProductApprox

/-! Substochastic bounds and local perturbations for actual empty and impulse kernels. -/
noncomputable section
namespace Spin.Structured.ConcreteEncoder
open ConcreteMaps FiniteKernel
open scoped symmDiff

theorem endpointKernel_substochastic {R : Nat} {z : ℝ} (hz : 0 ≤ z) (hz1 : z ≤ 1)
    (xs : Fin R → Input) : Substochastic (endpointKernel z xs) := by
  constructor
  · exact endpointKernel_nonneg hz xs
  · intro q
    rw [endpointKernel_total]
    calc inputMoment z xs q ≤ (Spin.piPMF (fun _ : Fin R => transvectionLaw)).expect (fun _ => 1) :=
        Spin.FinPMF.expect_mono _ (fun _ => pow_le_one₀ hz hz1)
      _ = 1 := Spin.FinPMF.expect_const _ _

theorem pathKernel_substochastic {z : ℝ} (hz : 0 ≤ z) (hz1 : z ≤ 1) (xs : List Input) :
    Substochastic (pathKernel z xs) := by
  have h := endpointKernel_substochastic hz hz1 xs.get
  simpa only [endpointKernel_eq_path, List.ofFn_get] using h

theorem roundKernel_substochastic {z : ℝ} (hz : 0 ≤ z) (hz1 : z ≤ 1) (x : Input) :
    Substochastic (roundKernel z x) := by
  simpa only [pathKernel, mul_one] using pathKernel_substochastic hz hz1 [x]

theorem matrixExpect_substochastic {α : Type*} [Fintype α] (P : Spin.FinPMF α)
    (K : α → Matrix State State ℝ) (h : ∀ x, Substochastic (K x)) :
    Substochastic (matrixExpect P K) := by
  constructor
  · intro q r
    exact P.expect_nonneg (fun x => (h x).nonneg q r)
  · intro q
    change (∑ r, P.expect (fun x => K x q r)) ≤ 1
    rw [← Spin.FinPMF.expect_sum]
    calc _ ≤ P.expect (fun _ => 1) := P.expect_mono (fun x => (h x).row_le q)
      _ = 1 := P.expect_const 1

theorem uniformImpulseKernel_substochastic {z : ℝ} (hz : 0 ≤ z) (hz1 : z ≤ 1) :
    Substochastic (uniformImpulseKernel z) :=
  matrixExpect_substochastic _ _ (fun i => roundKernel_substochastic hz hz1 {i})

theorem emptyKernel_substochastic {z : ℝ} (hz : 0 ≤ z) (hz1 : z ≤ 1) (g : Nat) :
    Substochastic (Matrix.of (emptyKernel z g)) :=
  endpointKernel_substochastic hz hz1 (fun _ : Fin g => ∅)

theorem liftedEmpty_substochastic {c : ℝ} (hc : 0 ≤ c) (g : Nat) :
    Substochastic (liftedEmpty c g) := by
  have he0 : 0 ≤ Real.exp (-c * (epochMean * g)) := (Real.exp_pos _).le
  have he1 : Real.exp (-c * (epochMean * g)) ≤ 1 := by
    apply Real.exp_le_one_iff.mpr
    have hm : 0 ≤ epochMean * (g : ℝ) := by unfold epochMean; positivity
    nlinarith
  constructor
  · intro q r
    rw [liftedEmpty_apply]
    split_ifs <;> positivity
  · intro q
    simp only [liftedEmpty_apply]
    by_cases hq : q = ∅
    · simp only [hq, ite_true, Finset.sum_ite_eq', Finset.mem_univ, le_refl]
    · simp only [hq, ite_false]
      have hh := liveMean_as_sum (fun _ : State => Real.exp (-c * (epochMean * g)))
      have hm : liveMean (fun _ : State => Real.exp (-c * (epochMean * g))) =
          Real.exp (-c * (epochMean * g)) := by
        unfold liveMean
        rw [Finset.sum_const, nsmul_eq_mul, nonzeroStates_card_actual]
        push_cast
        ring
      rw [hm] at hh
      calc
        _ = ∑ r : State, (if r = ∅ then 0 else (1 / 524287 : ℝ)) *
            Real.exp (-c * (epochMean * g)) := by
          apply Finset.sum_congr rfl
          intro r _
          split_ifs <;> ring
        _ = Real.exp (-c * (epochMean * g)) := hh.symm
        _ ≤ 1 := he1

theorem emptyKernel_rowError {c : ℝ} (hc : 0 ≤ c) (g : Nat) :
    RowError (Matrix.of (emptyKernel (Real.exp (-c)) g)) (liftedEmpty c g)
      (524288 * (c * (128 * Real.sqrt (3 * g)) + (1 / 2 : ℝ) ^ g)) := by
  have h := rowError_of_entry (K := Matrix.of (emptyKernel (Real.exp (-c)) g))
    (L := liftedEmpty c g) (emptyKernel_lifted_error hc g)
  simpa only [State, Fintype.card_finset, Fintype.card_fin, Nat.cast_pow, Nat.cast_ofNat,
    show (2 : ℝ) ^ 19 = 524288 by norm_num] using h

theorem roundKernel_exp_rowError {c : ℝ} (hc : 0 ≤ c) (x : Input) :
    RowError (roundKernel (Real.exp (-c)) x) (roundKernel 1 x) (128 * c) := by
  intro q
  let w := (x ∆ Aset q).card
  have hw : w ≤ 128 := by simpa only [Fintype.card_fin] using Finset.card_le_univ (x ∆ Aset q)
  have hwR : (w : ℝ) ≤ 128 := by exact_mod_cast hw
  have he0 : 0 ≤ Real.exp (-c * w) := (Real.exp_pos _).le
  have he1 : Real.exp (-c * w) ≤ 1 := by
    apply Real.exp_le_one_iff.mpr
    have hw0 : (0 : ℝ) ≤ w := Nat.cast_nonneg _
    nlinarith
  have hsum : (∑ r, Spin.Imt.actLaw q (r ∆ Cset x)) = 1 := by
    have h := nextState_expect q x (fun _ => 1)
    simpa only [Spin.FinPMF.expect_const, mul_one] using h.symm
  have hnonneg (r : State) : 0 ≤ Spin.Imt.actLaw q (r ∆ Cset x) :=
    div_nonneg (Nat.cast_nonneg _) (Nat.cast_nonneg _)
  have hd := exp_neg_lipschitz (mul_nonneg hc (Nat.cast_nonneg w)) (by norm_num : (0 : ℝ) ≤ 0)
  simp only [neg_zero, Real.exp_zero, sub_zero, abs_of_nonneg (mul_nonneg hc (Nat.cast_nonneg w))] at hd
  rw [← neg_mul] at hd
  have hb : 1 - Real.exp (-c * w) ≤ 128 * c := by
    rw [abs_of_nonpos (sub_nonpos.mpr he1)] at hd
    have hmul := mul_le_mul_of_nonneg_left hwR hc
    nlinarith
  calc rowAbs (roundKernel (Real.exp (-c)) x - roundKernel 1 x) q =
      ∑ r, (1 - Real.exp (-c * w)) * Spin.Imt.actLaw q (r ∆ Cset x) := by
        apply Finset.sum_congr rfl
        intro r _
        simp only [Matrix.sub_apply, roundKernel, emitted, exp_neg_pow, one_pow, one_mul]
        change |Real.exp (-c * w) * Spin.Imt.actLaw q (r ∆ Cset x) - Spin.Imt.actLaw q (r ∆ Cset x)| = _
        rw [abs_of_nonpos (by nlinarith [hnonneg r])]
        ring
    _ = 1 - Real.exp (-c * w) := by rw [← Finset.mul_sum, hsum, mul_one]
    _ ≤ 128 * c := hb

theorem matrixExpect_rowError {α : Type*} [Fintype α] (P : Spin.FinPMF α)
    (K L : α → Matrix State State ℝ) (ε : ℝ) (h : ∀ x, RowError (K x) (L x) ε) :
    RowError (matrixExpect P K) (matrixExpect P L) ε := by
  intro q
  calc rowAbs (matrixExpect P K - matrixExpect P L) q ≤
      ∑ r, P.expect (fun x => |K x q r - L x q r|) :=
        Finset.sum_le_sum (fun r _ => expect_abs_sub_le P _ _)
    _ = P.expect (fun x => rowAbs (K x - L x) q) := by
      rw [← Spin.FinPMF.expect_sum]
      rfl
    _ ≤ P.expect (fun _ => ε) := P.expect_mono (fun x => h x q)
    _ = ε := P.expect_const ε

theorem uniformImpulseKernel_exp_rowError {c : ℝ} (hc : 0 ≤ c) :
    RowError (uniformImpulseKernel (Real.exp (-c))) (uniformImpulseKernel 1) (128 * c) :=
  matrixExpect_rowError _ _ _ _ (fun i => roundKernel_exp_rowError hc {i})

end Spin.Structured.ConcreteEncoder
