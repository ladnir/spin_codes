import SpinCodes.UniformMixer.ScalableEncoderMoment

/-!
# Realized encoder and the chronological regional kernel

The kernel used by the regional comparison is the exact averaged
weighted moment of the sampled encoder, for every deterministic input
sequence and every entering state.
-/

noncomputable section

namespace Spin.UniformMixer

theorem sampledInputMoment_kernel {R : ℕ} (z : ℝ)
    (xs : Fin R → Input64) (q : State22) :
    sampledInputMoment z xs q =
      ∑ r : State22, inputListKernel z (List.ofFn xs) q r := by
  induction R generalizing q with
  | zero =>
      rw [sampledInputMoment_zero]
      simp [inputListKernel, identity]
  | succ R ih =>
      rw [sampledInputMoment_succ]
      simp only [List.ofFn_succ, inputListKernel]
      simp_rw [ih]
      simp only [compose, Finset.mul_sum]
      rw [Finset.sum_comm]
      rfl

theorem sampled_region_moment {R : ℕ} (z : ℝ)
    (T : Finset (Fin (64 * R))) (q : State22) :
    sampledInputMoment z (fun i => regionInput T i) q =
      ∑ r : State22, chronologicalRegionKernel z T q r := by
  simpa only [chronologicalRegionKernel, chronologicalInputList] using
    sampledInputMoment_kernel z (fun i => regionInput T i) q

end Spin.UniformMixer
