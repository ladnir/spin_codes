import SpinCodes.UniformMixer.ScalableEnvelopeIterate

/-!
# The Scalable routed moment and the three-coordinate transfer

This is the algebraic bridge from the sampled encoder to the paper's
three-coordinate bound. The remaining distance proof must evaluate
the transfer parameters, sum over the selected outer spectrum, and
close the dense, sparse, and fixed-occupation numerical gates.
-/

noncomputable section

namespace Spin.UniformMixer

open Spin.Structured Spin.Structured.ConcreteOuter
  Spin.Structured.ConcreteNativeFamily Spin.Structured.ConcreteRoute

def scalableRouteFactor (m : ℕ)
    (rows : Fin (Lsched m) → Finset (Fin (bsched m))) : ℝ :=
  ((bsched m : ℝ) + 1) ^ activeRows rows *
    ((Lsched m : ℝ) + 1) ^ bsched m

def scalableTransferMoment (β z : ℝ)
    (hβ0 : 0 ≤ β) (hβ1 : β ≤ 1) (n : ℕ) : ℝ :=
  let H := scalableZeroHold β z hβ0 hβ1
  let G := scalableLiveMoment β z hβ0 hβ1 0 - H
  let F := scalableMaxLiveMoment β z hβ0 hβ1
  let U := scalableMeanLiveMoment β z hβ0 hβ1
  let u := scalableFeedbackActive β hβ0 hβ1
  ∑ i : Fin 3,
    (((weightedTransfer H G F U (min F u) (min U u)) ^ n).mulVec
      ![1, 0, 0]) i

theorem scalable_selected_transfer_moment_bound (m : ℕ)
    {z : ℝ} (hz0 : 0 ≤ z) (hz1 : z ≤ 1)
    {out : NativeSeed m} (hg : nativeGood m out)
    {message : Fin (Nsched m / 2) → Bool} (hm : message ≠ 0) :
    (law (nativeRows m out message)).expect
      (scalableRoutedMoment m z) ≤
    scalableRouteFactor m (nativeRows m out message) *
      scalableTransferMoment
        (density (nativeRows m out message)) z
        (density_nonneg (nativeRows m out message))
        (density_le_one (Lsched_pos m) (bsched_pos m)
          (nativeRows m out message))
        (bsched m * scalableStepsPerRegion m) := by
  let rows := nativeRows m out message
  let β := density rows
  let hβ0 := density_nonneg rows
  let hβ1 := density_le_one (Lsched_pos m) (bsched_pos m) rows
  let n := bsched m * scalableStepsPerRegion m
  change (law rows).expect (scalableRoutedMoment m z) ≤
    scalableRouteFactor m rows * scalableTransferMoment β z hβ0 hβ1 n
  calc
    _ ≤ scalableRouteFactor m rows *
        ∑ r : State22, (bernoulliStepKernel β z hβ0 hβ1 ^ n) 0 r :=
      scalable_selected_outer_moment_bound m hz0 hg hm
    _ ≤ scalableRouteFactor m rows *
        scalableTransferMoment β z hβ0 hβ1 n := by
      exact mul_le_mul_of_nonneg_left
        (scalable_bernoulli_transfer_bound β z hβ0 hβ1 hz0 hz1 n)
        (by unfold scalableRouteFactor; positivity)

end Spin.UniformMixer
