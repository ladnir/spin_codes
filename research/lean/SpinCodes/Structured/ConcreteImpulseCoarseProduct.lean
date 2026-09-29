import SpinCodes.Structured.ConcreteImpulseComposition
import SpinCodes.Structured.ConcreteEpochTilt

/-! Exact lift/projection algebra for the full-state impulse kernel. -/
noncomputable section
namespace Spin.Structured.ConcreteEncoder
open ConcreteMaps

def liveLift : Matrix State (Fin 2) ℝ := fun q b =>
  coarseValue (if b = 0 then 1 else 0) (if b = 1 then 1 else 0) q

def liveProjection : Matrix (Fin 2) State ℝ := fun a q =>
  if a = 0 then (if q = ∅ then 1 else 0) else (if q = ∅ then 0 else 1 / 524287)

theorem liveMean_as_sum (f : State → ℝ) :
    liveMean f = ∑ q, (if q = ∅ then 0 else (1 / 524287 : ℝ)) * f q := by
  unfold liveMean Spin.nonzeroStates
  rw [Finset.sum_filter, Finset.sum_div]
  apply Finset.sum_congr rfl
  intro q _
  by_cases hq : q = ∅ <;> simp [hq] <;> ring

theorem liveProjection_action (f : State → ℝ) (a : Fin 2) :
    (∑ q, liveProjection a q * f q) = if a = 0 then f ∅ else liveMean f := by
  by_cases ha : a = 0
  · simp only [liveProjection, ha, ite_true, ite_mul, one_mul, zero_mul,
      Finset.sum_ite_eq', Finset.mem_univ]
  · simp only [liveProjection, ha, ite_false]
    exact (liveMean_as_sum f).symm

theorem liveMean_coarseValue (a b : ℝ) : liveMean (coarseValue a b) = b := by
  unfold liveMean coarseValue
  rw [Finset.sum_congr rfl (fun q hq =>
    if_neg (Spin.ne_empty_of_mem_nonzeroStates hq))]
  rw [Finset.sum_const, nsmul_eq_mul, nonzeroStates_card_actual]
  push_cast
  ring

theorem liveProjection_lift : liveProjection * liveLift = 1 := by
  ext a b
  rw [Matrix.mul_apply, liveProjection_action]
  have hl : liveMean (fun q => liveLift q b) = (if b = 1 then 1 else 0) :=
    liveMean_coarseValue _ _
  rw [hl]
  fin_cases a <;> fin_cases b <;>
    norm_num [liveLift, coarseValue, Matrix.one_apply]

theorem coarse_product_apply (K : Matrix State State ℝ) (a b : Fin 2) :
    (liveProjection * K * liveLift) a b =
      if a = 0 then (∑ r, K ∅ r * liveLift r b)
      else liveMean (fun q => ∑ r, K q r * liveLift r b) := by
  rw [Matrix.mul_assoc, Matrix.mul_apply, liveProjection_action]
  rfl

theorem roundKernel_one_action (i : Fin 128) (q : State) (f : State → ℝ) :
    (∑ r, roundKernel 1 {i} q r * f r) =
      transvectionLaw.expect (fun t => f (nextState q {i} t)) := by
  rw [nextState_expect q {i} f]
  simp only [roundKernel, emitted, one_pow, one_mul]

theorem roundKernel_one_coarse (i : Fin 128) :
    liveProjection * roundKernel 1 {i} * liveLift = impulseMatrix := by
  ext a b
  rw [coarse_product_apply]
  simp only [roundKernel_one_action]
  change impulseCoarse i a b = impulseMatrix a b
  rw [impulseCoarse_eq]

theorem matrixExpect_rect_left {α ι κ : Type*} [Fintype α] [Fintype ι] [Fintype κ]
    (P : Spin.FinPMF α) (K : Matrix ι State ℝ) (f : α → Matrix State κ ℝ)
    (a : ι) (b : κ) :
    P.expect (fun x => (K * f x) a b) = ∑ q, K a q * P.expect (fun x => f x q b) := by
  simp only [Matrix.mul_apply, Spin.FinPMF.expect_sum, expect_const_mul]

theorem uniformImpulseKernel_one_coarse :
    liveProjection * uniformImpulseKernel 1 * liveLift = impulseMatrix := by
  ext a b
  have he : (liveProjection * uniformImpulseKernel 1 * liveLift) a b =
      coordinateLaw.expect (fun i => (liveProjection * roundKernel 1 {i} * liveLift) a b) := by
    simp only [Matrix.mul_apply, uniformImpulseKernel, matrixExpect,
      Spin.FinPMF.expect_sum, expect_const_mul, expect_mul_const]
  rw [he]
  simp only [roundKernel_one_coarse, Spin.FinPMF.expect_const]

def coarseEmpty (c : ℝ) (g : Nat) : Matrix (Fin 2) (Fin 2) ℝ :=
  !![1, 0; 0, Real.exp (-c * (epochMean * g))]

def liftedEmpty (c : ℝ) (g : Nat) : Matrix State State ℝ :=
  liveLift * coarseEmpty c g * liveProjection

theorem liftedEmpty_apply (c : ℝ) (g : Nat) (q r : State) :
    liftedEmpty c g q r = if q = ∅ then (if r = ∅ then 1 else 0)
      else (if r = ∅ then 0 else Real.exp (-c * (epochMean * g)) / 524287) := by
  by_cases hq : q = ∅ <;> by_cases hr : r = ∅ <;>
    simp [liftedEmpty, liveLift, liveProjection, coarseValue, coarseEmpty,
      Matrix.mul_apply, Fin.sum_univ_two, hq, hr] <;> ring

theorem emptyKernel_lifted_error {c : ℝ} (hc : 0 ≤ c) (g : Nat) (q r : State) :
    |emptyKernel (Real.exp (-c)) g q r - liftedEmpty c g q r| ≤
      c * (128 * Real.sqrt (3 * g)) + (1 / 2 : ℝ) ^ g := by
  rw [liftedEmpty_apply]
  by_cases hq : q = ∅
  · subst q
    rw [emptyKernel_zero]
    simp only [ite_true, sub_self, abs_zero]
    positivity
  · rw [if_neg hq]
    by_cases hr : r = ∅
    · subst r
      rw [emptyKernel_live_zero _ _ hq]
      simp only [ite_true, sub_self, abs_zero]
      positivity
    · rw [if_neg hr]
      exact emptyKernel_tilt_live hc g hq hr

end Spin.Structured.ConcreteEncoder
