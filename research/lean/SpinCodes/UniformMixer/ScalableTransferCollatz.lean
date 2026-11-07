import SpinCodes.UniformMixer.ScalableTransferFirstMoment

/-!
# Positive-vector bounds for the three-coordinate transfer

A positive Collatz witness bounds every power of a nonnegative transfer
matrix. The dense and sparse numerical gates can supply their own
parameters and witnesses to this lemma.
-/

noncomputable section

namespace Spin.UniformMixer

open Spin.Structured Spin.Structured.ConcreteOuter
  Spin.Structured.ConcreteNativeFamily Spin.Structured.ConcreteRoute

theorem nonnegative_mulVec_mono
    (T : Matrix (Fin 3) (Fin 3) ℝ)
    (hT : ∀ i j, 0 ≤ T i j)
    (x y : Fin 3 → ℝ)
    (hxy : ∀ j, x j ≤ y j) (i : Fin 3) :
    T.mulVec x i ≤ T.mulVec y i := by
  change ∑ j : Fin 3, T i j * x j ≤
    ∑ j : Fin 3, T i j * y j
  apply Finset.sum_le_sum
  intro j _
  exact mul_le_mul_of_nonneg_left (hxy j) (hT i j)

theorem weightedTransfer_nonneg
    (K : Matrix State22 State22 ℝ) (H G F U V W : ℝ)
    (B : WeightedRowBounds K H G F U V W)
    (i j : Fin 3) :
    0 ≤ weightedTransfer H G F U V W i j := by
  have hD : (0 : ℝ) ≤ (liveCount State22 : ℝ) :=
    state22_liveCount_real_pos.le
  fin_cases i <;> fin_cases j <;>
    simp only [weightedTransfer] <;>
    first
    | exact B.H_nonneg
    | exact B.G_nonneg
    | exact B.F_nonneg
    | exact B.U_nonneg
    | exact div_nonneg B.V_nonneg hD
    | exact div_nonneg B.W_nonneg hD
    | exact le_refl _

theorem transfer_collatz_pointwise
    (T : Matrix (Fin 3) (Fin 3) ℝ)
    (hT : ∀ i j, 0 ≤ T i j)
    (v e : Fin 3 → ℝ) (c ρ : ℝ)
    (hc : 0 ≤ c) (hρ : 0 ≤ ρ)
    (he : ∀ i, e i ≤ c * v i)
    (hv : ∀ i, T.mulVec v i ≤ ρ * v i)
    (n : ℕ) (i : Fin 3) :
    (T ^ n).mulVec e i ≤ c * ρ ^ n * v i := by
  induction n generalizing i with
  | zero =>
      simpa only [pow_zero, Matrix.one_mulVec, mul_one] using he i
  | succ n ih =>
      calc
        _ = T.mulVec ((T ^ n).mulVec e) i := by
          rw [pow_succ', Matrix.mulVec_mulVec]
        _ ≤ T.mulVec (fun j => (c * ρ ^ n) * v j) i :=
          nonnegative_mulVec_mono T hT _ _ ih i
        _ = (c * ρ ^ n) * T.mulVec v i := by
          change T.mulVec ((c * ρ ^ n) • v) i = _
          rw [Matrix.mulVec_smul]
          rfl
        _ ≤ (c * ρ ^ n) * (ρ * v i) :=
          mul_le_mul_of_nonneg_left (hv i)
            (mul_nonneg hc (pow_nonneg hρ _))
        _ = c * ρ ^ (n + 1) * v i := by
          rw [pow_succ]
          ring

theorem transfer_collatz_total
    (T : Matrix (Fin 3) (Fin 3) ℝ)
    (hT : ∀ i j, 0 ≤ T i j)
    (v e : Fin 3 → ℝ) (c ρ : ℝ)
    (hc : 0 ≤ c) (hρ : 0 ≤ ρ)
    (he : ∀ i, e i ≤ c * v i)
    (hv : ∀ i, T.mulVec v i ≤ ρ * v i)
    (n : ℕ) :
    ∑ i : Fin 3, (T ^ n).mulVec e i ≤
      c * ρ ^ n * ∑ i : Fin 3, v i := by
  calc
    _ ≤ ∑ i : Fin 3, c * ρ ^ n * v i := by
      exact Finset.sum_le_sum (fun i _ =>
        transfer_collatz_pointwise T hT v e c ρ hc hρ he hv n i)
    _ = _ := by rw [Finset.mul_sum]

theorem scalableThreeTransfer_nonneg (β z : ℝ)
    (hβ0 : 0 ≤ β) (hβ1 : β ≤ 1)
    (hz0 : 0 ≤ z) (hz1 : z ≤ 1)
    (i j : Fin 3) :
    0 ≤ scalableThreeTransfer β z hβ0 hβ1 i j := by
  exact weightedTransfer_nonneg _ _ _ _ _ _ _
    (scalable_actual_weighted_rows β z hβ0 hβ1 hz0 hz1) i j

theorem scalable_bernoulli_collatz_bound (β z : ℝ)
    (hβ0 : 0 ≤ β) (hβ1 : β ≤ 1)
    (hz0 : 0 ≤ z) (hz1 : z ≤ 1)
    (v : Fin 3 → ℝ) (c ρ : ℝ)
    (hc : 0 ≤ c) (hρ : 0 ≤ ρ)
    (he : ∀ i, (![1, 0, 0] : Fin 3 → ℝ) i ≤ c * v i)
    (hv : ∀ i, (scalableThreeTransfer β z hβ0 hβ1).mulVec v i ≤
      ρ * v i) (n : ℕ) :
    ∑ r : State22, (bernoulliStepKernel β z hβ0 hβ1 ^ n) 0 r ≤
      c * ρ ^ n * ∑ i : Fin 3, v i := by
  calc
    _ ≤ scalableTransferMoment β z hβ0 hβ1 n :=
      scalable_bernoulli_transfer_bound β z hβ0 hβ1 hz0 hz1 n
    _ ≤ c * ρ ^ n * ∑ i : Fin 3, v i := by
      exact transfer_collatz_total _
        (scalableThreeTransfer_nonneg β z hβ0 hβ1 hz0 hz1)
        v ![1, 0, 0] c ρ hc hρ he hv n

/-- A positive witness for the actual three-coordinate transfer bounds
the routed moment of each message under the selected shared outer. -/
theorem scalable_selected_collatz_moment_bound (m : ℕ)
    {z : ℝ} (hz0 : 0 ≤ z) (hz1 : z ≤ 1)
    {out : NativeSeed m} (hg : nativeGood m out)
    {message : Fin (Nsched m / 2) → Bool} (hm : message ≠ 0)
    (v : Fin 3 → ℝ) (c ρ : ℝ)
    (hc : 0 ≤ c) (hρ : 0 ≤ ρ)
    (he : ∀ i, (![1, 0, 0] : Fin 3 → ℝ) i ≤ c * v i)
    (hv : ∀ i,
      (scalableThreeTransfer
        (density (nativeRows m out message)) z
        (density_nonneg (nativeRows m out message))
        (density_le_one (Lsched_pos m) (bsched_pos m)
          (nativeRows m out message))).mulVec v i ≤ ρ * v i) :
    (law (nativeRows m out message)).expect
      (scalableRoutedMoment m z) ≤
    scalableRouteFactor m (nativeRows m out message) *
      (c * ρ ^ (bsched m * scalableStepsPerRegion m) *
        ∑ i : Fin 3, v i) := by
  let rows := nativeRows m out message
  let β := density rows
  let hβ0 := density_nonneg rows
  let hβ1 := density_le_one (Lsched_pos m) (bsched_pos m) rows
  let n := bsched m * scalableStepsPerRegion m
  have hfactor : 0 ≤ scalableRouteFactor m rows := by
    unfold scalableRouteFactor
    positivity
  calc
    _ ≤ scalableRouteFactor m rows *
        scalableTransferMoment β z hβ0 hβ1 n :=
      scalable_selected_transfer_moment_bound m hz0 hz1 hg hm
    _ ≤ scalableRouteFactor m rows *
        (c * ρ ^ n * ∑ i : Fin 3, v i) := by
      apply mul_le_mul_of_nonneg_left _ hfactor
      exact transfer_collatz_total _
        (scalableThreeTransfer_nonneg β z hβ0 hβ1 hz0 hz1)
        v ![1, 0, 0] c ρ hc hρ he hv n

end Spin.UniformMixer
