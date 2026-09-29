import SpinCodes.Structured.ConcreteShufflePin
import SpinCodes.Structured.ConcreteRoutePin
import SpinCodes.Structured.ConcreteEncoderPin
import SpinCodes.Structured.ConcreteRoutedMoment

open Spin Spin.Imt Spin.Structured
open ConcreteRoute ConcreteRoutedEncoder

example {L b : ℕ} (hL : 0 < L) (hb : 0 < b) (hdiv : 128 ∣ b * L)
    (rows : Fin L → Finset (Fin b)) (hw0 : 0 < totalWeight rows)
    (hw1 : totalWeight rows < L * b) {z : ℝ} (hz0 : 0 ≤ z) (hz1 : z ≤ 1) :
    ((seedLaw L b).prod (piPMF (fun _ : Fin (b * L / 128) =>
      ConcreteEncoder.transvectionLaw))).expect
        (fun ω => z ^ ConcreteEncoder.outputWeight
          (reshape (streamWiring hdiv) (routeEval rows ω.1)) ω.2 ∅) ≤
      (((b : ℝ) + 1) ^ activeRows rows * ((L : ℝ) + 1) ^ b) *
        (((Occupation.Sparse.numericalMatrix
          ((totalWeight rows : ℝ) / ((L : ℝ) * b)) z).apply)^[b * L / 128]
          (Coords.eZ 5)).total :=
  stream_moment_bound hL hb hdiv rows hw0 hw1 hz0 hz1

example {L b R : ℕ} (hL : 0 < L) (hb : 0 < b)
    (e : (Fin b × Fin L) ≃ (Fin R × Fin 128))
    (rows : Fin L → Finset (Fin b)) {α : ℝ}
    (hα0 : 0 < α) (hα1 : α ≤ 1 / 10000) (hq : density rows = (4 / 5) * α) :
    (ConcreteRoutedEncoder.experimentLaw L b R).expect
        (fun ω => (1 - (8 / 5) * α) ^ weight e rows ω) ≤
      (((b : ℝ) + 1) ^ activeRows rows * ((L : ℝ) + 1) ^ b) *
        (2048 * (1 - 96 * α) ^ R) :=
  sparse_moment_bound hL hb e rows hα0 hα1 hq

#print axioms Spin.Structured.Routing.uniform_permutation_shuffleLaw
#print axioms Spin.Structured.Routing.permutedLaw_eq_shuffledLaw
#print axioms Spin.Structured.Routing.shuffleLaw_dominates_rowBernoulli
#print axioms Spin.Structured.Routing.shuffleLaw_prob_card
#print axioms Spin.Structured.Routing.shuffledLaw_fiber_uniform
#print axioms Spin.Structured.Routing.shuffledLaw_prob_card
#print axioms Spin.Structured.ConcreteRoute.transpose_bernoulli
#print axioms Spin.Structured.ConcreteRoute.law_dominates
#print axioms Spin.Structured.ConcreteRoute.permutation_law_eq
#print axioms Spin.Structured.ConcreteRoute.permutation_dominates
#print axioms Spin.Structured.ConcreteRoute.permutation_expect_le
#print axioms Spin.Structured.ConcreteRoute.permutation_expect_le_of_weight
#print axioms Spin.Structured.ConcreteRoute.routeEval_totalWeight
#print axioms Spin.Structured.ConcreteEncoder.inputMoment_nonneg
#print axioms Spin.Structured.ConcreteEncoder.transvection_expect
#print axioms Spin.Structured.ConcreteEncoder.nextState_expect
#print axioms Spin.Structured.ConcreteEncoder.averagedMoment_eq_iterate
#print axioms Spin.Structured.ConcreteEncoder.encoder_moment_bound
#print axioms Spin.Structured.ConcreteEncoder.sparse_encoder_moment
#print axioms Spin.Structured.ConcreteEncoder.sparse_encoder_experiment
#print axioms Spin.Structured.ConcreteRoute.reshape_iid
#print axioms Spin.Structured.ConcreteRoute.regionMajor_position
#print axioms Spin.Structured.ConcreteRoutedEncoder.moment_eq
#print axioms Spin.Structured.ConcreteRoutedEncoder.moment_bound
#print axioms Spin.Structured.ConcreteRoutedEncoder.sparse_moment_bound
#print axioms Spin.Structured.ConcreteRoutedEncoder.stream_moment_bound
