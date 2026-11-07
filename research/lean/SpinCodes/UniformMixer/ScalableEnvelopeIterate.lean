import SpinCodes.UniformMixer.ScalableEnvelopeRows

/-!
# Iterating the weighted-state envelope

The initial state is zero. A checked row bound propagates the three
coordinates through any number of steps, and their terminal sum
dominates the full-state weighted moment.
-/

noncomputable section

namespace Spin.UniformMixer

theorem weightedEnvelope_initial :
    WeightedEnvelope (fun q : State22 => if q = 0 then 1 else 0)
      1 0 0 := by
  refine ⟨fun _ => 0, ?_, ?_, ?_, ?_⟩
  · intro q
    exact le_refl _
  · simp only [Finset.sum_const_zero, le_refl]
  · simp only [ite_true, le_refl]
  · intro q hq
    simp only [hq, ite_false, zero_add, zero_div, le_refl]

theorem weightedEnvelope_total_le (μ : State22 → ℝ)
    (a f u : ℝ) (h : WeightedEnvelope μ a f u) :
    ∑ q : State22, μ q ≤ a + f + u := by
  obtain ⟨ν, _hν, hmass, hzero, hlive⟩ := h
  have hsplit := Finset.add_sum_erase (Finset.univ : Finset State22)
    μ (Finset.mem_univ (0 : State22))
  change μ 0 + ∑ q ∈ liveStateSet, μ q = ∑ q : State22, μ q at hsplit
  rw [← hsplit]
  calc
    _ ≤ a + ∑ q ∈ liveStateSet,
        (ν q + u / (liveCount State22 : ℝ)) := by
      apply add_le_add hzero
      apply Finset.sum_le_sum
      intro q hq
      exact hlive q (Finset.ne_of_mem_erase hq)
    _ = a + (∑ q ∈ liveStateSet, ν q) + u := by
      simp only [Finset.sum_add_distrib, liveStateSet_sum_constant]
      ring
    _ ≤ a + f + u := by linarith

def weightedBudgetStep (H G F U V W : ℝ)
    (v : Fin 3 → ℝ) : Fin 3 → ℝ :=
  ![v 0 * H + v 1 * V / (liveCount State22 : ℝ) +
      v 2 * W / (liveCount State22 : ℝ),
    v 0 * G,
    v 1 * F + v 2 * U]

def weightedTransfer (H G F U V W : ℝ) :
    Matrix (Fin 3) (Fin 3) ℝ :=
  fun i j =>
    if i = 0 then
      if j = 0 then H else if j = 1 then
        V / (liveCount State22 : ℝ) else
        W / (liveCount State22 : ℝ)
    else if i = 1 then
      if j = 0 then G else 0
    else if j = 0 then 0 else if j = 1 then F else U

theorem weightedBudgetStep_matrix (H G F U V W : ℝ)
    (v : Fin 3 → ℝ) :
    weightedBudgetStep H G F U V W v =
      (weightedTransfer H G F U V W).mulVec v := by
  funext i
  fin_cases i
  · change v 0 * H + v 1 * V / (liveCount State22 : ℝ) +
      v 2 * W / (liveCount State22 : ℝ) =
      (weightedTransfer H G F U V W).mulVec v 0
    simp only [Matrix.mulVec, dotProduct, Fin.sum_univ_three]
    change v 0 * H + v 1 * V / (liveCount State22 : ℝ) +
      v 2 * W / (liveCount State22 : ℝ) =
      H * v 0 + (V / (liveCount State22 : ℝ)) * v 1 +
        (W / (liveCount State22 : ℝ)) * v 2
    ring
  · change v 0 * G = (weightedTransfer H G F U V W).mulVec v 1
    simp [weightedTransfer, Matrix.mulVec, dotProduct]
    ring
  · change v 1 * F + v 2 * U =
      (weightedTransfer H G F U V W).mulVec v 2
    simp [weightedTransfer, Matrix.mulVec, dotProduct,
      Fin.sum_univ_three]
    ring

def weightedBudget (H G F U V W : ℝ) : ℕ → Fin 3 → ℝ
  | 0 => ![1, 0, 0]
  | n + 1 => weightedBudgetStep H G F U V W
      (weightedBudget H G F U V W n)

theorem weightedBudget_matrix_power (H G F U V W : ℝ)
    (n : ℕ) :
    weightedBudget H G F U V W n =
      ((weightedTransfer H G F U V W) ^ n).mulVec ![1, 0, 0] := by
  induction n with
  | zero =>
      simp only [weightedBudget, pow_zero, Matrix.one_mulVec]
  | succ n ih =>
      calc
        weightedBudget H G F U V W (n + 1) =
            (weightedTransfer H G F U V W).mulVec
              (weightedBudget H G F U V W n) := by
          exact weightedBudgetStep_matrix _ _ _ _ _ _ _
        _ = (weightedTransfer H G F U V W).mulVec
            (((weightedTransfer H G F U V W) ^ n).mulVec
              ![1, 0, 0]) := by rw [ih]
        _ = (((weightedTransfer H G F U V W) *
              (weightedTransfer H G F U V W) ^ n)).mulVec
              ![1, 0, 0] :=
          Matrix.mulVec_mulVec _ _ _
        _ = _ := by rw [pow_succ']

theorem weightedBudget_nonneg
    (H G F U V W : ℝ)
    (hH : 0 ≤ H) (hG : 0 ≤ G) (hF : 0 ≤ F)
    (hU : 0 ≤ U) (hV : 0 ≤ V) (hW : 0 ≤ W)
    (n : ℕ) (i : Fin 3) :
    0 ≤ weightedBudget H G F U V W n i := by
  induction n generalizing i with
  | zero =>
      fin_cases i <;> norm_num [weightedBudget]
  | succ n ih =>
      have hD : (0 : ℝ) ≤ (liveCount State22 : ℝ) :=
        state22_liveCount_real_pos.le
      fin_cases i
      · change 0 ≤ weightedBudget H G F U V W n 0 * H +
          weightedBudget H G F U V W n 1 * V /
            (liveCount State22 : ℝ) +
          weightedBudget H G F U V W n 2 * W /
            (liveCount State22 : ℝ)
        exact add_nonneg
          (add_nonneg (mul_nonneg (ih 0) hH)
            (div_nonneg (mul_nonneg (ih 1) hV) hD))
          (div_nonneg (mul_nonneg (ih 2) hW) hD)
      · change 0 ≤ weightedBudget H G F U V W n 0 * G
        exact mul_nonneg (ih 0) hG
      · change 0 ≤ weightedBudget H G F U V W n 1 * F +
          weightedBudget H G F U V W n 2 * U
        exact add_nonneg (mul_nonneg (ih 1) hF)
          (mul_nonneg (ih 2) hU)

theorem weightedEnvelope_iterate
    (K : Matrix State22 State22 ℝ) (H G F U V W : ℝ)
    (B : WeightedRowBounds K H G F U V W) (n : ℕ) :
    WeightedEnvelope (fun r => (K ^ n) 0 r)
      (weightedBudget H G F U V W n 0)
      (weightedBudget H G F U V W n 1)
      (weightedBudget H G F U V W n 2) := by
  induction n with
  | zero =>
      change WeightedEnvelope
        (fun r : State22 => if (0 : State22) = r then 1 else 0)
        1 0 0
      simpa only [eq_comm] using weightedEnvelope_initial
  | succ n ih =>
      let v := weightedBudget H G F U V W n
      have ha : 0 ≤ v 0 :=
        weightedBudget_nonneg H G F U V W
          B.H_nonneg B.G_nonneg B.F_nonneg B.U_nonneg
          B.V_nonneg B.W_nonneg n 0
      have hf : 0 ≤ v 1 :=
        weightedBudget_nonneg H G F U V W
          B.H_nonneg B.G_nonneg B.F_nonneg B.U_nonneg
          B.V_nonneg B.W_nonneg n 1
      have hu : 0 ≤ v 2 :=
        weightedBudget_nonneg H G F U V W
          B.H_nonneg B.G_nonneg B.F_nonneg B.U_nonneg
          B.V_nonneg B.W_nonneg n 2
      have hstep := weightedEnvelope_step K H G F U V W B
        (fun r => (K ^ n) 0 r) (v 0) (v 1) (v 2)
        ha hf hu ih
      have hpow :
          weightedKernelStep K (fun r => (K ^ n) 0 r) =
            (fun r => (K ^ (n + 1)) 0 r) := by
        funext r
        rw [pow_succ, Matrix.mul_apply]
        rfl
      rw [hpow] at hstep
      exact hstep

theorem weighted_kernel_power_total_bound
    (K : Matrix State22 State22 ℝ) (H G F U V W : ℝ)
    (B : WeightedRowBounds K H G F U V W) (n : ℕ) :
    ∑ r : State22, (K ^ n) 0 r ≤
      weightedBudget H G F U V W n 0 +
      weightedBudget H G F U V W n 1 +
      weightedBudget H G F U V W n 2 :=
  weightedEnvelope_total_le _ _ _ _
    (weightedEnvelope_iterate K H G F U V W B n)

theorem weighted_kernel_power_transfer_bound
    (K : Matrix State22 State22 ℝ) (H G F U V W : ℝ)
    (B : WeightedRowBounds K H G F U V W) (n : ℕ) :
    ∑ r : State22, (K ^ n) 0 r ≤
      ∑ i : Fin 3,
        (((weightedTransfer H G F U V W) ^ n).mulVec
          ![1, 0, 0]) i := by
  calc
    _ ≤ weightedBudget H G F U V W n 0 +
        weightedBudget H G F U V W n 1 +
        weightedBudget H G F U V W n 2 :=
      weighted_kernel_power_total_bound K H G F U V W B n
    _ = _ := by
      rw [weightedBudget_matrix_power]
      simp only [Fin.sum_univ_three]

theorem scalable_bernoulli_power_three_bound (β z : ℝ)
    (hβ0 : 0 ≤ β) (hβ1 : β ≤ 1)
    (hz0 : 0 ≤ z) (hz1 : z ≤ 1) (n : ℕ) :
    ∑ r : State22,
      (bernoulliStepKernel β z hβ0 hβ1 ^ n) 0 r ≤
    let H := scalableZeroHold β z hβ0 hβ1
    let G := scalableLiveMoment β z hβ0 hβ1 0 - H
    let F := scalableMaxLiveMoment β z hβ0 hβ1
    let U := scalableMeanLiveMoment β z hβ0 hβ1
    let u := scalableFeedbackActive β hβ0 hβ1
    weightedBudget H G F U (min F u) (min U u) n 0 +
      weightedBudget H G F U (min F u) (min U u) n 1 +
      weightedBudget H G F U (min F u) (min U u) n 2 := by
  exact weighted_kernel_power_total_bound _ _ _ _ _ _ _
    (scalable_actual_weighted_rows β z hβ0 hβ1 hz0 hz1) n

theorem scalable_bernoulli_transfer_bound (β z : ℝ)
    (hβ0 : 0 ≤ β) (hβ1 : β ≤ 1)
    (hz0 : 0 ≤ z) (hz1 : z ≤ 1) (n : ℕ) :
    ∑ r : State22,
      (bernoulliStepKernel β z hβ0 hβ1 ^ n) 0 r ≤
    let H := scalableZeroHold β z hβ0 hβ1
    let G := scalableLiveMoment β z hβ0 hβ1 0 - H
    let F := scalableMaxLiveMoment β z hβ0 hβ1
    let U := scalableMeanLiveMoment β z hβ0 hβ1
    let u := scalableFeedbackActive β hβ0 hβ1
    ∑ i : Fin 3,
      (((weightedTransfer H G F U (min F u) (min U u)) ^ n).mulVec
        ![1, 0, 0]) i := by
  exact weighted_kernel_power_transfer_bound _ _ _ _ _ _ _
    (scalable_actual_weighted_rows β z hβ0 hβ1 hz0 hz1) n

end Spin.UniformMixer
