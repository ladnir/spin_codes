import SpinCodes.UniformMixer.ScalableStepLaw

/-!
# Moment of the sampled regional encoder

This is the realized 64-bit-step encoder, with one independently sampled
invertible state map per step. The weighted moment obeys the same
`stepKernel` recurrence as the regional proof.
-/

noncomputable section

namespace Spin.UniformMixer

abbrev Input64 := Fin 64 → ZMod 2

/-- Output weight of the actual deterministic encoder for fixed state maps. -/
def sampledOutputWeight : {R : ℕ} →
    (Fin R → Input64) → (Fin R → StateMap) → State22 → ℕ
  | 0, _, _, _ => 0
  | _R + 1, xs, gs, q =>
      stepWeight (xs 0) q +
        sampledOutputWeight (Fin.tail xs) (Fin.tail gs)
          ((gs 0) q + stepFeedback (xs 0))

def sampledInputMoment {R : ℕ} (z : ℝ)
    (xs : Fin R → Input64) (q : State22) : ℝ :=
  (piPMF (fun _ : Fin R => FinPMF.uniform StateMap)).expect
    (fun gs => z ^ sampledOutputWeight xs gs q)

theorem sampledInputMoment_zero (z : ℝ)
    (xs : Fin 0 → Input64) (q : State22) :
    sampledInputMoment z xs q = 1 := by
  exact FinPMF.expect_const _ 1

theorem sampledInputMoment_succ {R : ℕ} (z : ℝ)
    (xs : Fin (R + 1) → Input64) (q : State22) :
    sampledInputMoment z xs q =
      ∑ r : State22,
        stepKernel z (xs 0) q r *
          sampledInputMoment z (Fin.tail xs) r := by
  have hsplit : sampledInputMoment z xs q =
      z ^ stepWeight (xs 0) q *
        (FinPMF.uniform StateMap).expect (fun g =>
          sampledInputMoment z (Fin.tail xs)
            (g q + stepFeedback (xs 0))) := by
    unfold sampledInputMoment
    rw [Spin.Structured.ConcreteEncoder.expect_iid_succ]
    simp only [sampledOutputWeight, Fin.cons_zero, Fin.tail_cons, pow_add]
    simp only [FinPMF.expect, Finset.mul_sum]
    apply Finset.sum_congr rfl
    intro g _
    apply Finset.sum_congr rfl
    intro gs _
    ring
  rw [hsplit, uniform_state_map_expect]
  simp only [stepKernel, Finset.mul_sum]
  apply Finset.sum_congr rfl
  intro r _
  ring

end Spin.UniformMixer
