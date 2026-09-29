import SpinCodes.Structured.ConcreteEpochProductApproxLocal

/-! The actual finite gap/impulse experiment approximates its two-state
product, with an explicit additive error summed over all gaps and impulses. -/
noncomputable section
namespace Spin.Structured.ConcreteEncoder
open FiniteKernel

def approximateImpulseProduct (c : ℝ) : {a : Nat} → (Fin a → Nat) → Nat → Matrix State State ℝ
  | 0, _, last => liftedEmpty c last
  | _ + 1, gaps, last => liftedEmpty c (gaps 0) * uniformImpulseKernel 1 *
      approximateImpulseProduct c (Fin.tail gaps) last

def coarseImpulseProduct (c : ℝ) : {a : Nat} → (Fin a → Nat) → Nat → Matrix (Fin 2) (Fin 2) ℝ
  | 0, _, last => coarseEmpty c last
  | _ + 1, gaps, last => coarseEmpty c (gaps 0) * impulseMatrix *
      coarseImpulseProduct c (Fin.tail gaps) last

theorem approximateImpulseProduct_eq_lift {a : Nat} (c : ℝ) (gaps : Fin a → Nat) (last : Nat) :
    approximateImpulseProduct c gaps last =
      liveLift * coarseImpulseProduct c gaps last * liveProjection := by
  induction a with
  | zero => rfl
  | succ a ih =>
    rw [approximateImpulseProduct, coarseImpulseProduct, ih, liftedEmpty]
    calc
      _ = liveLift * coarseEmpty c (gaps 0) *
          (liveProjection * uniformImpulseKernel 1 * liveLift) *
            coarseImpulseProduct c (Fin.tail gaps) last * liveProjection := by
        simp only [Matrix.mul_assoc]
      _ = _ := by
        rw [uniformImpulseKernel_one_coarse]
        simp only [Matrix.mul_assoc]

def emptyError (c : ℝ) (g : Nat) : ℝ :=
  524288 * (c * (128 * Real.sqrt (3 * g)) + (1 / 2 : ℝ) ^ g)

def impulseErrorBudget (c : ℝ) : {a : Nat} → (Fin a → Nat) → Nat → ℝ
  | 0, _, last => emptyError c last
  | _ + 1, gaps, last => emptyError c (gaps 0) + 128 * c + impulseErrorBudget c (Fin.tail gaps) last

theorem emptyError_nonneg {c : ℝ} (hc : 0 ≤ c) (g : Nat) : 0 ≤ emptyError c g := by
  unfold emptyError
  positivity

theorem impulseErrorBudget_nonneg {a : Nat} {c : ℝ} (hc : 0 ≤ c)
    (gaps : Fin a → Nat) (last : Nat) : 0 ≤ impulseErrorBudget c gaps last := by
  induction a with
  | zero => exact emptyError_nonneg hc last
  | succ a ih =>
    change 0 ≤ emptyError c (gaps 0) + 128 * c + impulseErrorBudget c (Fin.tail gaps) last
    exact add_nonneg (add_nonneg (emptyError_nonneg hc _) (by positivity)) (ih _)

theorem impulseErrorBudget_sum {a : Nat} (c : ℝ) (gaps : Fin a → Nat) (last : Nat) :
    impulseErrorBudget c gaps last = emptyError c last + ∑ i, (emptyError c (gaps i) + 128 * c) := by
  induction a with
  | zero => simp only [impulseErrorBudget, Finset.univ_eq_empty, Finset.sum_empty, add_zero]
  | succ a ih =>
    rw [impulseErrorBudget, ih, Fin.sum_univ_succ]
    have ht : (∑ i : Fin a, (emptyError c (Fin.tail gaps i) + 128 * c)) =
        ∑ i : Fin a, (emptyError c (gaps i.succ) + 128 * c) := rfl
    rw [ht]
    ring

theorem uniformImpulseProduct_substochastic {a : Nat} {z : ℝ} (hz : 0 ≤ z) (hz1 : z ≤ 1)
    (gaps : Fin a → Nat) (last : Nat) : Substochastic (uniformImpulseProduct z gaps last) := by
  induction a with
  | zero => exact emptyKernel_substochastic hz hz1 last
  | succ a ih =>
    exact ((emptyKernel_substochastic hz hz1 (gaps 0)).mul
      (uniformImpulseKernel_substochastic hz hz1)).mul (ih (Fin.tail gaps))

theorem approximateImpulseProduct_substochastic {a : Nat} {c : ℝ} (hc : 0 ≤ c)
    (gaps : Fin a → Nat) (last : Nat) : Substochastic (approximateImpulseProduct c gaps last) := by
  induction a with
  | zero => exact liftedEmpty_substochastic hc last
  | succ a ih =>
    exact ((liftedEmpty_substochastic hc (gaps 0)).mul
      (uniformImpulseKernel_substochastic (by norm_num) (by norm_num))).mul (ih (Fin.tail gaps))

/-- Each empty-gap error and each impulse error occurs once; the estimate
does not amplify exponentially with the number of factors. -/
theorem uniformImpulseProduct_rowError {a : Nat} {c : ℝ} (hc : 0 ≤ c)
    (gaps : Fin a → Nat) (last : Nat) :
    RowError (uniformImpulseProduct (Real.exp (-c)) gaps last)
      (approximateImpulseProduct c gaps last) (impulseErrorBudget c gaps last) := by
  have hz : 0 ≤ Real.exp (-c) := (Real.exp_pos _).le
  have hz1 : Real.exp (-c) ≤ 1 := Real.exp_le_one_iff.mpr (by linarith)
  induction a with
  | zero => exact emptyKernel_rowError hc last
  | succ a ih =>
    have hlocal := RowError.mul (emptyKernel_rowError hc (gaps 0))
      (uniformImpulseKernel_exp_rowError hc)
      (uniformImpulseKernel_substochastic hz hz1)
      (liftedEmpty_substochastic hc (gaps 0)) (by positivity : 0 ≤ 128 * c)
    exact RowError.mul hlocal (ih (Fin.tail gaps))
      (uniformImpulseProduct_substochastic hz hz1 (Fin.tail gaps) last)
      ((liftedEmpty_substochastic hc (gaps 0)).mul
        (uniformImpulseKernel_substochastic (by norm_num) (by norm_num)))
      (impulseErrorBudget_nonneg hc (Fin.tail gaps) last)

theorem rowError_entry {K L : Matrix State State ℝ} {ε : ℝ} (h : RowError K L ε) (q r : State) :
    |K q r - L q r| ≤ ε := by
  calc
    |K q r - L q r| ≤ rowAbs (K - L) q := by
      exact Finset.single_le_sum (f := fun s => |K q s - L q s|)
        (fun _ _ => abs_nonneg _) (Finset.mem_univ r)
    _ ≤ ε := h q

theorem uniformImpulseProduct_coarse_error {a : Nat} {c : ℝ} (hc : 0 ≤ c)
    (gaps : Fin a → Nat) (last : Nat) (q r : State) :
    |uniformImpulseProduct (Real.exp (-c)) gaps last q r -
      (liveLift * coarseImpulseProduct c gaps last * liveProjection) q r| ≤
        impulseErrorBudget c gaps last := by
  rw [← approximateImpulseProduct_eq_lift]
  exact rowError_entry (uniformImpulseProduct_rowError hc gaps last) q r

/-- Actual random-coordinate, random-transvection experiment, with no
unproved kernel-identification premise. -/
theorem endpointKernel_uniform_impulses_coarse_error {a : Nat} {θ L : ℝ}
    (hθ : 0 ≤ θ) (hL : 0 < L) (gaps : Fin a → Nat) (last : Nat) (q r : State) :
    |(matrixExpect (Spin.piPMF (fun _ : Fin a => coordinateLaw))
        (fun coords => endpointKernel (Real.exp (-(θ / L))) (impulseInputs gaps coords last).get)) q r -
      (liveLift * coarseImpulseProduct (θ / L) gaps last * liveProjection) q r| ≤
        impulseErrorBudget (θ / L) gaps last := by
  rw [endpointKernel_uniform_impulses]
  exact uniformImpulseProduct_coarse_error (div_nonneg hθ hL.le) gaps last q r

end Spin.Structured.ConcreteEncoder
