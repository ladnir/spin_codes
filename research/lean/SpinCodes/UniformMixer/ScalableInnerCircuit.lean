import SpinCodes.LinearCircuitScan
import SpinCodes.UniformMixer.ScalableCodeProperties

/-! A chronological circuit for the actual sampled 22-state, 64-bit inner.
The state is carried through the whole scan, including region boundaries. -/
noncomputable section
set_option maxHeartbeats 2000000
namespace Spin.UniformMixer
open Spin.LinearCircuit

theorem stepFeedback_smul (c : Bit) (x : Input64) :
    stepFeedback (c • x) = c • stepFeedback x := by
  funext i
  simp only [stepFeedback, Pi.smul_apply, smul_eq_mul, Finset.mul_sum]
  apply Finset.sum_congr rfl
  intro j _
  exact mul_assoc c (x j) (coordCoeff j i)
/-- One actual local update, with output coordinates followed by the next state. -/
def scalableStepLinearMap (g : StateMap) :
    (Sum (Fin 22) (Fin 64) → Bit) →ₗ[Bit]
      (Sum (Fin 64) (Fin 22) → Bit) where
  toFun v := Sum.elim
    (fun j => coordMap j (fun i => v (.inl i)) + v (.inr j))
    (fun i => (g (fun i => v (.inl i))) i +
      stepFeedback (fun j => v (.inr j)) i)
  map_add' v w := by
    funext z
    cases z with
    | inl j =>
      change coordMap j ((fun i => v (.inl i)) + (fun i => w (.inl i))) +
        (v (.inr j) + w (.inr j)) =
        (coordMap j (fun i => v (.inl i)) + v (.inr j)) +
          (coordMap j (fun i => w (.inl i)) + w (.inr j))
      rw [(coordMap j).map_add]
      exact add_add_add_comm _ _ _ _
    | inr i =>
      change (g ((fun i => v (.inl i)) + (fun i => w (.inl i)))) i +
        stepFeedback ((fun j => v (.inr j)) + (fun j => w (.inr j))) i =
        ((g (fun i => v (.inl i))) i + stepFeedback (fun j => v (.inr j)) i) +
          ((g (fun i => w (.inl i))) i + stepFeedback (fun j => w (.inr j)) i)
      rw [g.map_add, stepFeedback_add]
      exact add_add_add_comm _ _ _ _
  map_smul' c v := by
    funext z
    cases z with
    | inl j =>
      change coordMap j (c • (fun i => v (.inl i))) + c • v (.inr j) =
        c • (coordMap j (fun i => v (.inl i)) + v (.inr j))
      rw [(coordMap j).map_smul]
      exact (smul_add c _ _).symm
    | inr i =>
      change (g (c • (fun i => v (.inl i)))) i +
        stepFeedback (c • (fun j => v (.inr j))) i =
        c • ((g (fun i => v (.inl i))) i + stepFeedback (fun j => v (.inr j)) i)
      rw [g.map_smul, stepFeedback_smul]
      exact (smul_add c _ _).symm
def scalableStepCircuit (g : StateMap) :
    Circuit (Sum (Fin 22) (Fin 64)) (Sum (Fin 64) (Fin 22)) :=
  Circuit.ofLinearMap (scalableStepLinearMap g)

theorem scalableStepCircuit_eval (g : StateMap) (q : State22) (x : Input64) :
    (scalableStepCircuit g).eval (Sum.elim q x) =
      Sum.elim (fun j => coordMap j q + x j) (g q + stepFeedback x) := by
  rw [scalableStepCircuit, Circuit.eval_ofLinearMap]
  rfl

theorem scalableStepCircuit_size (g : StateMap) :
    (scalableStepCircuit g).size = 14792 := by
  simp [scalableStepCircuit, Circuit.size_ofLinearMap, Fintype.card_sum]

def scalableInnerCircuit {R : ℕ} (gs : Fin R → StateMap) :
    Circuit (Fin R × Fin 64) (Fin R × Fin 64) :=
  zeroScan (fun i => scalableStepCircuit (gs i))

/-- The scan retains exactly the output of the continuous sampled recursion. -/
theorem scalableInner_runScan {R : ℕ} (gs : Fin R → StateMap)
    (xs : Fin R → Input64) (q : State22) :
    (runScan (fun i => scalableStepCircuit (gs i)) xs q).1 = sampledOutput xs gs q := by
  induction R generalizing q with
  | zero => rfl
  | succ R ih =>
    simp only [runScan, scalableStepCircuit_eval, Sum.elim_inr, Sum.elim_inl,
      sampledOutput]
    funext i
    refine Fin.cases rfl (fun k => ?_) i
    exact congrFun (ih (Fin.tail gs) (Fin.tail xs)
      ((gs 0) q + stepFeedback (xs 0))) k

theorem scalableInnerCircuit_eval {R : ℕ} (gs : Fin R → StateMap)
    (xs : Fin R → Input64) :
    (scalableInnerCircuit gs).eval (fun ix => xs ix.1 ix.2) =
      fun iy => sampledOutput xs gs 0 iy.1 iy.2 := by
  rw [scalableInnerCircuit, zeroScan_eval, scalableInner_runScan]

theorem scalableInnerCircuit_size {R : ℕ} (gs : Fin R → StateMap) :
    (scalableInnerCircuit gs).size = 14792 * R := by
  simp [scalableInnerCircuit, zeroScan_size, scalableStepCircuit_size,
    Nat.mul_comm]

#print axioms scalableInnerCircuit_eval
#print axioms scalableInnerCircuit_size
end Spin.UniformMixer
