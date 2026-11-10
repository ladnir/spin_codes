import SpinCodes.Turbo.PhysicalPacketFourier22

/-!
# Zero feedback under the actual complete packet-input law

This is the same finite character-orthogonality argument as
`scalableZeroHold_fourier_sum`. It holds for any complete packet-word law;
no independence between output weight and feedback is assumed.
-/
noncomputable section
namespace Spin.Turbo
open Spin.UniformMixer Spin.Structured.ConcreteFourier

def physicalZeroHold22 (P : FinPMF PacketStep16) (z : ℝ) : ℝ :=
  P.expect (fun x => if stepFeedback (packetStepInput22 x) = 0 then
    z ^ stepWeight (packetStepInput22 x) 0 else 0)

theorem physicalZeroHold22_fourier_sum (P : FinPMF PacketStep16) (z : ℝ) :
    (2^22 : ℝ) * physicalZeroHold22 P z =
      ∑ q : State22, P.expect (fun x => z ^ stepWeight (packetStepInput22 x) 0 *
        (binarySign (statePairing q (stepFeedback (packetStepInput22 x))) : ℝ)) := by
  unfold physicalZeroHold22 FinPMF.expect
  rw [Finset.mul_sum]
  calc
    _ = ∑ x : PacketStep16, P.p x * (z ^ stepWeight (packetStepInput22 x) 0 *
        ∑ q : State22,
          (binarySign (statePairing q (stepFeedback (packetStepInput22 x))) : ℝ)) := by
      apply Finset.sum_congr rfl
      intro x _
      rw [← feedback_indicator_character (packetStepInput22 x)]
      by_cases hx : stepFeedback (packetStepInput22 x) = 0
      · simp only [hx, ite_true]; ring
      · simp only [hx, ite_false]; ring
    _ = _ := by
      simp only [Finset.mul_sum]
      rw [Finset.sum_comm]

theorem physicalZeroHold22_at_one (P : FinPMF PacketStep16) :
    physicalZeroHold22 P 1 = P.prob (fun x => stepFeedback (packetStepInput22 x) = 0) := by
  unfold physicalZeroHold22
  rw [FinPMF.prob_eq_expect_indicator]
  apply congrArg P.expect
  funext x
  simp only [one_pow]

theorem physicalFeedbackActive22_complement (P : FinPMF PacketStep16) :
    P.prob (fun x => stepFeedback (packetStepInput22 x) ≠ 0) =
      1 - physicalZeroHold22 P 1 := by
  rw [FinPMF.prob_compl, physicalZeroHold22_at_one]

private theorem stepKernel22_zero_indicator (z : ℝ) (x : Input64) :
    stepKernel z x 0 0 = if stepFeedback x = 0 then z ^ stepWeight x 0 else 0 := by
  simp [stepKernel, update_zero, eq_comm]

theorem physicalInputKernel22_zero_zero (P : FinPMF PacketStep16) (z : ℝ) :
    physicalInputKernel22 P z 0 0 = physicalZeroHold22 P z := by
  simp only [physicalInputKernel22, Matrix.of_apply, physicalZeroHold22,
    stepKernel22_zero_indicator]

#print axioms physicalZeroHold22_fourier_sum
#print axioms physicalInputKernel22_zero_zero
end Spin.Turbo
