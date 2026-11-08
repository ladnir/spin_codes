import SpinCodes.UniformMixer.Generated.FullCover.ScalableDenseFullCover
import SpinCodes.UniformMixer.ScalableDistance
import SpinCodes.UniformMixer.ScalableComplexity
import SpinCodes.Structured.ScheduleUpper

/-!
# Checked native-length Scalable distance, rate, and circuit bounds

The generated rational cover discharges the final numerical hypothesis.
The event uses one shared sampled outer, one sampled regional route, and
the actual continuous GL state recursion for all messages.

Build the generated cover using `scripts/build_scalable_dense_packed_cover.ps1`
before compiling this module. Its generated inputs are proposed externally;
all leaf inequalities and coverage are checked by Lean.
The imported scalar-gate circuits compute the actual encoder and its
transpose with at most `272 * Nsched m` gates for every fixed setup.
-/

noncomputable section
attribute [local instance] Classical.propDecidable
namespace Spin.UniformMixer

open Spin.Structured Spin.Structured.ConcreteNativeFamily Filter

theorem scalable_failure_tendsto :
    Tendsto (fun m => (scalableSetupLaw m).prob (scalableBad m))
      atTop (nhds 0) := by
  exact scalable_failure_tendsto_of_dense_rates
    ScalableDenseFullCover.rates (by norm_num) (by norm_num)

#print axioms scalable_failure_tendsto

/-- The paper's native-length minimum-distance conclusion under the original law. -/
theorem scalable_minimum_distance_failure_tendsto :
    Tendsto (fun m => (scalableSetupLaw m).prob (fun seed =>
      scalableMinimumDistance m seed ≤ ⌊(0.11 : ℝ) * Nsched m⌋₊))
      atTop (nhds 0) := by
  simpa only [threshold_floor] using
    scalable_minimum_distance_failure_tendsto_of_dense_rates
      ScalableDenseFullCover.rates (by norm_num) (by norm_num)

#print axioms scalable_minimum_distance_failure_tendsto
#print axioms scalableRealizedCode_rate
#print axioms scalablePrimitiveCircuit_eval
#print axioms scalablePrimitiveCircuit_gates_le
#print axioms scalablePrimitiveTranspose_pairing
#print axioms scalablePrimitiveTranspose_gates_le
#print axioms bsched_logb_error_eventually

end Spin.UniformMixer
