import SpinCodes.UniformMixer.ScalableOuterCircuit
import SpinCodes.UniformMixer.ScalableRouteCircuit
import SpinCodes.UniformMixer.ScalableInnerCircuit
import SpinCodes.UniformMixer.ScalableLinearCode
import SpinCodes.LinearCircuitPrimitive

/-!
# Linear-size circuits for the actual Scalable encoder and its transpose

All sampled coefficients and permutation wires are fixed setup data.
The circuit follows the real outer, route, and continuous inner. The bound
is deliberately coarse and uniform over every setup and native length.
-/
noncomputable section
namespace Spin.UniformMixer

open Spin.Structured Spin.Structured.ConcreteOuter
  Spin.Structured.ConcreteBinaryEncoding Spin.LinearCircuit

def scalableCircuit (m : ℕ) (seed : ScalableSetupSeed m) :
    Circuit (Fin (Nsched m / 2))
      (Fin (bsched m * scalableStepsPerRegion m) × Fin 64) :=
  .comp (scalableInnerCircuit (scalableFlatMaps m seed.2.2))
    (.comp (scalableRouteCircuit m seed.2.1) (nativeOuterCircuit m seed.1))

theorem scalableCircuit_eval (m : ℕ) (seed : ScalableSetupSeed m)
    (x : Fin (Nsched m / 2) → Bit) :
    (scalableCircuit m seed).eval x =
      fun p => scalableCodeword m seed x p.1 p.2 := by
  have hout : nativeOuterVectorLinearMap m seed.1 x =
      fun p => bitEquiv (decide (p.2 ∈
        nativeRows m seed.1 ((wordEquiv (Nsched m / 2)).symm x) p.1)) := by
    funext p
    rw [nativeOuterVectorLinearMap_rows]
    simp only [rowsToWord, Equiv.symm_apply_apply]
  change (scalableInnerCircuit (scalableFlatMaps m seed.2.2)).eval
    ((scalableRouteCircuit m seed.2.1).eval ((nativeOuterCircuit m seed.1).eval x)) = _
  rw [nativeOuterCircuit_eval, hout, scalableRouteCircuit_eval, scalableInnerCircuit_eval]
  rfl

theorem scalableCircuit_size (m : ℕ) (seed : ScalableSetupSeed m) :
    (scalableCircuit m seed).size =
      14792 * (bsched m * scalableStepsPerRegion m) + 40 * Nsched m := by
  simp only [scalableCircuit, Circuit.size, scalableInnerCircuit_size,
    scalableRouteCircuit_size, nativeOuterCircuit_size, Nat.zero_add]

theorem scalableCircuit_size_le (m : ℕ) (seed : ScalableSetupSeed m) :
    (scalableCircuit m seed).size ≤ 272 * Nsched m := by
  rw [scalableCircuit_size]
  have hn := scalable_output_length m
  omega

theorem scalableTransposeCircuit_size_le (m : ℕ) (seed : ScalableSetupSeed m) :
    (scalableCircuit m seed).transpose.size ≤ 272 * Nsched m := by
  rw [Circuit.size_transpose]
  exact scalableCircuit_size_le m seed

/-- The reverse circuit is the adjoint of the actual complete encoder. -/
theorem scalableTransposeCircuit_pairing (m : ℕ) (seed : ScalableSetupSeed m)
    (x : Fin (Nsched m / 2) → Bit)
    (y : (Fin (bsched m * scalableStepsPerRegion m) × Fin 64) → Bit) :
    ((scalableCircuit m seed).transpose.eval y) ⬝ᵥ x =
      y ⬝ᵥ (fun p => scalableCodeword m seed x p.1 p.2) := by
  rw [← scalableCircuit_eval]
  exact Circuit.transpose_pairing _ x y

/-- Explicit coefficient-AND/XOR circuit, with the macros expanded. -/
def scalablePrimitiveCircuit (m : ℕ) (seed : ScalableSetupSeed m) :=
  (scalableCircuit m seed).lower

def scalablePrimitiveTranspose (m : ℕ) (seed : ScalableSetupSeed m) :=
  (scalableCircuit m seed).transpose.lower

theorem scalablePrimitiveCircuit_eval (m : ℕ) (seed : ScalableSetupSeed m)
    (x : Fin (Nsched m / 2) → Bit) :
    (scalablePrimitiveCircuit m seed).eval x =
      fun p => scalableCodeword m seed x p.1 p.2 := by
  rw [scalablePrimitiveCircuit, Circuit.lower_eval, scalableCircuit_eval]

theorem scalablePrimitiveCircuit_gates_le (m : ℕ) (seed : ScalableSetupSeed m) :
    (scalablePrimitiveCircuit m seed).gates ≤ 272 * Nsched m := by
  rw [scalablePrimitiveCircuit, Circuit.lower_gates]
  exact scalableCircuit_size_le m seed

theorem scalablePrimitiveTranspose_gates_le (m : ℕ) (seed : ScalableSetupSeed m) :
    (scalablePrimitiveTranspose m seed).gates ≤ 272 * Nsched m := by
  rw [scalablePrimitiveTranspose, Circuit.lower_gates]
  exact scalableTransposeCircuit_size_le m seed

theorem scalablePrimitiveTranspose_pairing (m : ℕ) (seed : ScalableSetupSeed m)
    (x : Fin (Nsched m / 2) → Bit)
    (y : (Fin (bsched m * scalableStepsPerRegion m) × Fin 64) → Bit) :
    ((scalablePrimitiveTranspose m seed).eval y) ⬝ᵥ x =
      y ⬝ᵥ (fun p => scalableCodeword m seed x p.1 p.2) := by
  rw [scalablePrimitiveTranspose, Circuit.lower_eval]
  exact scalableTransposeCircuit_pairing m seed x y

#print axioms scalableCircuit_eval
#print axioms scalableCircuit_size_le
#print axioms scalableTransposeCircuit_pairing
#print axioms scalableTransposeCircuit_size_le
#print axioms scalablePrimitiveCircuit_eval
#print axioms scalablePrimitiveCircuit_gates_le
#print axioms scalablePrimitiveTranspose_pairing
#print axioms scalablePrimitiveTranspose_gates_le

end Spin.UniformMixer
